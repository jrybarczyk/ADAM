![Bioconductor availability](https://bioconductor.org/shields/availability/release/ADAM.svg)
![Bioconductor downloads](https://bioconductor.org/shields/downloads/release/ADAM.svg)
![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/ADAM.svg)
![Bioconductor build status](https://bioconductor.org/shields/build/release/bioc/ADAM.svg)
![Bioconductor dependencies](https://bioconductor.org/shields/dependencies/release/ADAM.svg)

# ADAM: Activity and Diversity Analysis Module

**Authors:** Andre L. Molan, Giordano B. S. Seco, Agnes A. S. Takeda, Jose L. Rybarczyk-Filho  
**Maintainer:** Jose L. Rybarczyk-Filho

## What ADAM does

ADAM performs gene set enrichment analysis by quantifying two complementary
signals for each functional group (GFAG):

- **Diversity** (`H`, relative `h`)
- **Activity** (`N`, relative `n`)

Comparisons are always built as **control vs experiment** pairs.

The package now provides a high-level wrapper, `run_adam()`, that supports both
legacy expression inputs and Bioconductor containers.

## Installation

`ADAM` is available on Bioconductor.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ADAM")
```

Development version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("jrybarczyk/ADAM")
```

## Main API: `run_adam()`

`run_adam()` is the recommended entry point.

It accepts:

- `data.frame` or tab-delimited file path (legacy)
- `SummarizedExperiment`
- `SingleCellExperiment`
- `ExpressionSet`

It returns a `SummarizedExperiment` with:

- assay `metrics` (numeric result matrix)
- `rowData` with GFAG annotations
- `metadata()` containing `legacy_result`, `gene_function_map`, `comparison_id`,
  and `mode`

### Input requirements

- Expression matrix/table must have gene identifiers and at least two samples.
- ADAM comparisons use the form `"control_sample,experiment_sample"`.
- You must provide valid `DBSpecies`, `AnalysisDomain`, and `GeneIdentifier`
  compatible with your annotation source.

For the bundled Aedes example, use:

- `DBSpecies = KeggPathwaysAedes`
- `AnalysisDomain = "own"`
- `GeneIdentifier = "geneStableID"`

## Quick start (legacy table input)

```r
library(ADAM)
data(ExpressionAedes)
data(KeggPathwaysAedes)

res_partial <- run_adam(
  x = ExpressionAedes,
  mode = "partial",
  comparison_id = "control1,experiment1",
  MinGene = 3L,
  MaxGene = 20L,
  DBSpecies = KeggPathwaysAedes,
  AnalysisDomain = "own",
  GeneIdentifier = "geneStableID"
)

res_partial
SummarizedExperiment::assay(res_partial, "metrics")[1:5, ]
S4Vectors::metadata(res_partial)$comparison_id
```

## Quick start (`SummarizedExperiment` with auto-comparisons)

```r
library(ADAM)
library(SummarizedExperiment)
library(S4Vectors)

data(ExpressionAedes)
data(KeggPathwaysAedes)

expr <- as.matrix(ExpressionAedes[, -1])
rownames(expr) <- as.character(ExpressionAedes$gene)

col_data <- S4Vectors::DataFrame(
  condition = c("control", "experiment", "control", "experiment"),
  row.names = colnames(expr)
)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = expr),
  colData = col_data
)

res_se <- run_adam(
  x = se,
  mode = "partial",
  assay_name = "counts",
  group_col = "condition",
  contrast = c("control", "experiment"),
  DBSpecies = KeggPathwaysAedes,
  AnalysisDomain = "own",
  GeneIdentifier = "geneStableID"
)

res_se
head(S4Vectors::metadata(res_se)$gene_function_map)
```

## Complete mode (`GFAGAnalysis` pipeline)

Use `mode = "complete"` to enable bootstrap-based significance (plus optional
Wilcoxon/Fisher tests).

```r
res_complete <- run_adam(
  x = ExpressionAedes,
  mode = "complete",
  comparison_id = "control1,experiment1",
  MinGene = 3L,
  MaxGene = 20L,
  SeedNumber = 1049,
  BootstrapNumber = 50L,
  PCorrection = 0.05,
  PCorrectionMethod = "fdr",
  WilcoxonTest = FALSE,
  FisherTest = FALSE,
  DBSpecies = KeggPathwaysAedes,
  AnalysisDomain = "own",
  GeneIdentifier = "geneStableID"
)

res_complete
SummarizedExperiment::assay(res_complete, "metrics")[1:5, ]
```

## Optional pseudo-bulk aggregation

`run_adam()` can aggregate samples before analysis using metadata columns.

```r
col_data2 <- S4Vectors::DataFrame(
  sample_id = c("S1", "S1", "S2", "S2"),
  cluster = c("A", "B", "A", "B"),
  condition = c("control", "control", "experiment", "experiment"),
  row.names = colnames(expr)
)

se2 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = expr),
  colData = col_data2
)

res_agg <- run_adam(
  x = se2,
  mode = "partial",
  assay_name = "counts",
  group_col = "condition",
  contrast = c("control", "experiment"),
  aggregate_by = c("sample_id", "cluster", "condition"),
  aggregate_fun = "sum",
  DBSpecies = KeggPathwaysAedes,
  AnalysisDomain = "own",
  GeneIdentifier = "geneStableID"
)
```

## Legacy APIs

These functions remain available and are still valid:

- `ADAnalysis()` for partial analysis
- `GFAGAnalysis()` for complete analysis

They return the legacy list structure (`gene-function map`, `comparison tables`).
`run_adam()` wraps these APIs and converts output to `SummarizedExperiment`.

## Documentation and support

- Function reference: `?run_adam`, `?ADAnalysis`, `?GFAGAnalysis`
- Vignettes: `vignettes/ADAM.Rmd`, `vignettes/BiocContainers.Rmd`
- Issues: <https://github.com/jrybarczyk/ADAM/issues>
