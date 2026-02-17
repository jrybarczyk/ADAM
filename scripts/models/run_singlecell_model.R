#!/usr/bin/env Rscript

source("scripts/models/common_model_utils.R")

ensure_packages(
    c("ADAM", "scRNAseq", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors"),
    bioc = TRUE
)

library(ADAM)
library(scRNAseq)
library(SingleCellExperiment)
library(SummarizedExperiment)

sce <- scRNAseq::ZeiselBrainData()

assay_name <- if ("counts" %in% names(SummarizedExperiment::assays(sce))) {
    "counts"
} else {
    names(SummarizedExperiment::assays(sce))[1]
}

candidate_group_cols <- c("level1class", "level2class", "cell_type1", "cell_type", "label.main", "cluster")
available <- candidate_group_cols[candidate_group_cols %in% colnames(SummarizedExperiment::colData(sce))]
if (length(available) == 0L) {
    stop("Could not find a supported cell-group column in colData(sce).")
}

group_col <- available[1]
selected_groups <- pick_top_two_groups(SummarizedExperiment::colData(sce)[[group_col]])

n_cells_per_group <- 200L
idx_control <- which(as.character(SummarizedExperiment::colData(sce)[[group_col]]) == selected_groups[1])[1:n_cells_per_group]
idx_experiment <- which(as.character(SummarizedExperiment::colData(sce)[[group_col]]) == selected_groups[2])[1:n_cells_per_group]
idx <- c(idx_control, idx_experiment)
idx <- idx[!is.na(idx)]

if (length(idx) < 40L) {
    stop("Not enough cells to build the single-cell template.")
}

sce_sub <- sce[, idx]
cell_group <- as.character(SummarizedExperiment::colData(sce_sub)[[group_col]])

pseudo_id <- character(length(cell_group))
for (grp in unique(cell_group)) {
    grp_idx <- which(cell_group == grp)
    pseudo_id[grp_idx] <- paste0(grp, "_rep", ((seq_along(grp_idx) - 1L) %% 4L) + 1L)
}

SummarizedExperiment::colData(sce_sub)$adam_group <- cell_group
SummarizedExperiment::colData(sce_sub)$pseudo_id <- pseudo_id

expr <- SummarizedExperiment::assay(sce_sub, assay_name)
keep_genes <- rowSums(expr) > 0
sce_sub <- sce_sub[keep_genes, ]
expr <- expr[keep_genes, , drop = FALSE]

own_db <- build_data_driven_db(expr, n_groups = 80L, prefix = "SCRNA")

result <- ADAM::run_adam(
    x = sce_sub,
    mode = "partial",
    assay_name = assay_name,
    group_col = "adam_group",
    contrast = selected_groups,
    aggregate_by = c("pseudo_id", "adam_group"),
    aggregate_fun = "sum",
    MinGene = 10L,
    MaxGene = 500L,
    DBSpecies = own_db,
    AnalysisDomain = "own",
    GeneIdentifier = "geneid"
)

dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file = "outputs/models/scrna_zeisel_partial.rds")
utils::write.csv(
    as.data.frame(SummarizedExperiment::rowData(result)),
    file = "outputs/models/scrna_zeisel_partial_table.csv",
    row.names = FALSE
)

message("Single-cell model completed: outputs/models/scrna_zeisel_partial.rds")
