# ADAM visualization scripts (diversity/activity focused)

These scripts generate **tile heatmaps (no smoothing)** with:

- functional groups on Y axis
- diversity/activity metrics on X axis (`H_*`, `N_*`, `h`, `n`)

## 1) From an existing `run_adam()` result

```bash
Rscript scripts/visualization/plot_entropy_activity_heatmap.R \
  --input outputs/models/rnaseq_airway_partial.rds \
  --output outputs/models/rnaseq_airway_entropy_activity_heatmap.png \
  --comparison sampleA,sampleB \
  --top 40
```

Arguments:

- `--input`: path to `.rds` from `run_adam()` (`SummarizedExperiment`)
- `--output`: output image path (`.png`)
- `--comparison`: comparison ID to plot (if omitted, first comparison is used)
- `--top`: number of functional groups to display

## 2) End-to-end demo (Aedes data from ADAM package)

```bash
Rscript scripts/visualization/demo_entropy_activity_heatmap.R
```

Outputs:

- `outputs/models/demo_aedes_partial_for_heatmap.rds`
- `outputs/models/demo_aedes_entropy_activity_heatmap.png`

## Notes

- This is intentionally **not** a volcano-style plot.
- The view is tailored to ADAM's core signals (diversity/activity).
- Values are z-scored by metric column for visual comparability across scales.
