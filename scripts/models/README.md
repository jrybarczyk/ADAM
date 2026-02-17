# ADAM analysis model templates (Bioconductor datasets)

This folder provides three runnable templates that use Bioconductor datasets and
execute ADAM through `run_adam()`.

## Templates

- `run_rnaseq_model.R`: bulk RNA-seq model using `airway::airway`
- `run_microarray_model.R`: microarray model using `ALL::ALL`
- `run_singlecell_model.R`: single-cell model using `scRNAseq::ZeiselBrainData`

## What each template does

1. Installs/loads required Bioconductor packages if missing
2. Loads a reference dataset from Bioconductor
3. Builds a data-driven `DBSpecies` annotation table (mode `AnalysisDomain = "own"`)
4. Runs `ADAM::run_adam()` in `mode = "partial"`
5. Saves outputs under `outputs/models/`

## Run examples

From the repository root:

```bash
Rscript scripts/models/run_rnaseq_model.R
Rscript scripts/models/run_microarray_model.R
Rscript scripts/models/run_singlecell_model.R
```

## Outputs

Each template writes:

- an `.rds` object with full `SummarizedExperiment` output
- a `.csv` table derived from `rowData(result)`

Files are created in `outputs/models/`.

## Notes

- These templates use **data-driven modules** as functional groups for `DBSpecies`.
  They are intended as reproducible starting points and can be replaced by
  curated annotation mappings when available.
- To run complete ADAM (bootstrap significance), change `mode = "partial"` to
  `mode = "complete"` and set complete-mode arguments (e.g. `BootstrapNumber`).
