# Phase 54 — Uncertainty Quantification Demo

**Version:** 1.1 | **Tier:** Standard | **Date:** 2026-03-26

## Goal
Use RF tree variance as a practical uncertainty estimate for pIC50 predictions.
Evaluate: (1) calibration (fraction within k*sigma), (2) whether sigma predicts actual error.

CLI: `python main.py --input data/compounds.csv`

Outputs: uncertainty_results.csv, uncertainty_plot.png

## Logic
- Compute ECFP4 fingerprints → RF (200 trees, LOO-CV)
- Per test point: uncertainty sigma = std across all 200 tree predictions
- Calibration: fraction of compounds where |y_true - y_hat| <= k * sigma (k=1, 1.5, 2)
- Pearson r between sigma and |error| to assess uncertainty utility
- Plot: actual vs predicted (bubble size = sigma) + sigma vs |error| scatter + sigma by family boxplot

## Results

| Metric | Value |
|---|---|
| LOO-CV R² | 0.729 |
| LOO-CV MAE | 0.268 |
| Mean sigma | 0.359 |
| Sigma range | [0.173, 0.590] |
| Pearson r (sigma vs |error|) | **-0.115** |
| Within 1.0*sigma | 71.1% |
| Within 1.5*sigma | 84.4% |
| Within 2.0*sigma | 88.9% |

## Key Concepts
- RF tree variance as uncertainty estimate (std across 200 tree predictions)
- Calibration analysis: fraction of predictions within k*sigma (k=1, 1.5, 2)
- Pearson correlation between predicted uncertainty (sigma) and actual error (|y - y_hat|)
- LOO-CV for unbiased prediction + uncertainty estimation

## Verification Checklist
- [x] LOO-CV R² = 0.729, MAE = 0.268 (consistent with Phase 37)
- [x] Pearson r (sigma vs |error|) = -0.115 (poor uncertainty-error correlation)
- [x] 71.1% within 1*sigma (close to Gaussian 68%)
- [x] uncertainty_results.csv and uncertainty_plot.png saved to output/

## Risks
- RF tree variance reflects ensemble disagreement, not true epistemic uncertainty
- Negative correlation (r=-0.115) means RF sigma is slightly anti-predictive of actual errors on this dataset

## Key Insights
- RF tree variance is poorly correlated with actual error (r=-0.115) — doesn't reliably flag uncertain predictions
- Despite poor correlation, calibration is reasonable: 71% within 1σ (Gaussian would give 68%)
- RF uncertainty tends to be highest for out-of-distribution samples (sparse feature space), not necessarily wrong predictions
- For production UQ on small datasets, conformal prediction (Phase 48) is more reliable than RF variance
- Combining: use conformal intervals for coverage guarantees + RF sigma to flag potential OOD compounds

## Deviations from Plan
- None; plan implemented as specified.
