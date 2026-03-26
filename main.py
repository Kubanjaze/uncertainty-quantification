import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score, mean_absolute_error
RDLogger.DisableLog("rdApp.*")

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}


def load_compounds(path):
    df = pd.read_csv(path)
    records, n_bad = [], 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            n_bad += 1
            continue
        try:
            pic50 = float(row["pic50"])
        except (KeyError, ValueError):
            continue
        if np.isnan(pic50):
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useChirality=True)
        fam = str(row["compound_name"]).split("_")[0]
        records.append({
            "compound_name": str(row["compound_name"]),
            "family": fam if fam in FAMILY_COLORS else "other",
            "pic50": pic50,
            "fp": list(fp),
        })
    print(f"  {len(records)} valid ({n_bad} skipped)")
    return pd.DataFrame(records)


def rf_tree_variance(model, X):
    """Uncertainty = variance across individual tree predictions."""
    tree_preds = np.array([tree.predict(X) for tree in model.estimators_])
    return tree_preds.mean(axis=0), tree_preds.std(axis=0)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--n-estimators", type=int, default=200)
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)
    X = np.array(df["fp"].tolist(), dtype=float)
    y = df["pic50"].values
    n = len(y)

    # LOO-CV: get uncertainty estimate for each compound
    loo = LeaveOneOut()
    y_pred_loo = np.zeros(n)
    y_std_loo = np.zeros(n)

    print(f"Running LOO-CV (n={n}, {args.n_estimators} trees)...")
    for train_idx, test_idx in loo.split(X):
        rf = RandomForestRegressor(n_estimators=args.n_estimators, random_state=42, n_jobs=-1)
        rf.fit(X[train_idx], y[train_idx])
        mu, sigma = rf_tree_variance(rf, X[test_idx])
        y_pred_loo[test_idx] = mu
        y_std_loo[test_idx] = sigma

    r2 = r2_score(y, y_pred_loo)
    mae = mean_absolute_error(y, y_pred_loo)
    residuals = np.abs(y - y_pred_loo)
    print(f"  LOO-CV R²={r2:.3f}, MAE={mae:.3f}")
    print(f"  Uncertainty: mean sigma={y_std_loo.mean():.3f}, range=[{y_std_loo.min():.3f}, {y_std_loo.max():.3f}]")

    # Correlation: does high uncertainty predict high error?
    corr = np.corrcoef(y_std_loo, residuals)[0, 1]
    print(f"  Pearson r (sigma vs |error|) = {corr:.3f}")

    # Calibration: within k*sigma?
    for k in [1.0, 1.5, 2.0]:
        fraction = (residuals <= k * y_std_loo).mean()
        print(f"  Fraction within {k}*sigma: {fraction:.3f}")

    # Save results
    df["y_pred_loo"] = y_pred_loo
    df["uncertainty_sigma"] = y_std_loo
    df["abs_error"] = residuals
    df["within_1sigma"] = (residuals <= y_std_loo).astype(int)
    df.drop(columns=["fp"]).to_csv(os.path.join(args.output_dir, "uncertainty_results.csv"), index=False)
    print(f"Saved: {args.output_dir}/uncertainty_results.csv")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: actual vs predicted, bubble size = uncertainty
    sizes = 30 + 200 * (y_std_loo / y_std_loo.max())
    for fam, color in FAMILY_COLORS.items():
        mask = df["family"] == fam
        if mask.sum() == 0:
            continue
        axes[0].scatter(y[mask], y_pred_loo[mask], c=color, s=sizes[mask],
                        label=fam, alpha=0.75, edgecolors="white", lw=0.4)
    lo, hi = y.min() - 0.2, y.max() + 0.2
    axes[0].plot([lo, hi], [lo, hi], "k--", lw=1)
    axes[0].set_xlabel("Actual pIC50", fontsize=10)
    axes[0].set_ylabel("Predicted pIC50 (LOO-CV)", fontsize=10)
    axes[0].set_title(f"Predictions (bubble size = uncertainty)\nR²={r2:.3f}, MAE={mae:.3f}", fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=7)
    axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)

    # Panel 2: uncertainty vs absolute error
    axes[1].scatter(y_std_loo, residuals, c=[FAMILY_COLORS.get(f, "#808080") for f in df["family"]],
                    s=55, alpha=0.8, edgecolors="white", lw=0.4)
    m, b = np.polyfit(y_std_loo, residuals, 1)
    x_fit = np.linspace(y_std_loo.min(), y_std_loo.max(), 50)
    axes[1].plot(x_fit, m * x_fit + b, "k--", lw=1.2, label=f"r={corr:.3f}")
    axes[1].set_xlabel("RF Uncertainty (tree std)", fontsize=10)
    axes[1].set_ylabel("|Actual - Predicted|", fontsize=10)
    axes[1].set_title("Uncertainty vs Error\n(Is uncertainty useful?)", fontsize=10, fontweight="bold")
    axes[1].legend(fontsize=9)
    axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)

    # Panel 3: uncertainty per family (boxplot)
    families = [f for f in FAMILY_COLORS if f in df["family"].values]
    data_list = [y_std_loo[df["family"].values == f] for f in families]
    bp = axes[2].boxplot(data_list, patch_artist=True)
    for patch, fam in zip(bp["boxes"], families):
        patch.set_facecolor(FAMILY_COLORS[fam])
        patch.set_alpha(0.8)
    axes[2].set_xticks(range(1, len(families) + 1))
    axes[2].set_xticklabels(families, fontsize=10)
    axes[2].set_ylabel("RF Uncertainty (tree std)", fontsize=10)
    axes[2].set_title("Uncertainty by Scaffold Family", fontsize=10, fontweight="bold")
    axes[2].spines["top"].set_visible(False); axes[2].spines["right"].set_visible(False)

    plt.suptitle("Uncertainty Quantification: RF Tree Variance (LOO-CV)", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "uncertainty_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {args.output_dir}/uncertainty_plot.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
