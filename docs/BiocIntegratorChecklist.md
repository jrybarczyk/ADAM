# ADAM Bioconductor Integration Checklist

Branch target: `devel`

This checklist tracks the transition to a Bioconductor-first workflow, with
clear priorities and an execution log.

## P0 (Critical)

- [~] Make `run_adam()` the primary interface across docs, examples, and vignettes.
- [ ] Keep Bioconductor container support as first-class input/output:
  `SummarizedExperiment`, `SingleCellExperiment`, `ExpressionSet`.
- [ ] Ensure outputs are consistently returned as `SummarizedExperiment`
  with stable fields in `assays`, `rowData`, and `metadata`.
- [ ] Keep `DESCRIPTION` Bioc-compliant:
  avoid duplicate author metadata, minimize hard dependencies,
  move optional packages to `Suggests`.
- [ ] Keep CI/BiocCheck green on every push to `devel`.
- [ ] Avoid generated artifacts in VCS (`.so`, `.o`, benchmark/check outputs,
  rendered files).

## P1 (High)

- [ ] Reduce non-Bioc dependency surface when possible
  (`dplyr`, `stringr`), preferring base R + Bioc core APIs where practical.
- [ ] Add a Bioconductor end-to-end vignette centered on `run_adam()`.
- [ ] Expand tests with real S4/Bioc object scenarios and edge cases.
- [ ] Add compatibility checks with common Bioc workflows
  (`DESeq2`, `edgeR`, `limma`, `scater`) at integration level.
- [ ] Standardize metadata naming conventions for easier interop
  (stable `rowData`/`metadata` contract).

## P2 (Medium)

- [ ] Add support guidance for large-scale backends
  (`DelayedArray`, `HDF5Array`).
- [ ] Evaluate `ExperimentHub`/`AnnotationHub` integration for annotation loading.
- [ ] Add migration notes from legacy API (`ADAnalysis`, `GFAGAnalysis`) to
  `run_adam()`.
- [ ] Add benchmark/report templates focused on Bioc-native workflows.

## Execution Log

- [x] 2026-02-17: Created `devel` and `adam-legacy` branches.
- [x] 2026-02-17: Added Bioconductor dataset model templates under `scripts/models/`.
- [x] 2026-02-17: Updated `vignettes/ADAM.Rmd` to make `run_adam()` the default workflow.
- [ ] Next: align man page examples and remaining legacy-first references with `run_adam()`-first messaging.

## Working Rules (devel)

- Keep feature work on `devel` only (no direct merge required yet).
- Use small, reviewable commits tied to one checklist item when possible.
- Run CI/BiocCheck after each P0/P1 change set.
