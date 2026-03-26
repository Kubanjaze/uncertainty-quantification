# uncertainty-quantification — Phase 54

RF tree variance as uncertainty quantification for pIC50 predictions. LOO-CV with 200 trees; uncertainty sigma per compound. Evaluates calibration (fraction within k*sigma) and whether high sigma predicts high error.

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv
```

## Outputs

- `output/uncertainty_results.csv` — per-compound predictions, sigma, abs error
- `output/uncertainty_plot.png` — actual vs predicted (bubble=uncertainty) + sigma vs error + uncertainty by family
