#!/usr/bin/env Rscript

source("scripts/models/common_model_utils.R")

ensure_packages(c("ADAM", "ALL", "Biobase"), bioc = TRUE)

library(ADAM)
library(ALL)
library(Biobase)

data("ALL", package = "ALL")
eset <- ALL::ALL
pheno <- Biobase::pData(eset)

candidate_group_cols <- c("mol.biol", "BT")
available <- candidate_group_cols[candidate_group_cols %in% colnames(pheno)]
if (length(available) == 0L) {
    stop("Could not find group columns 'mol.biol' or 'BT' in ALL phenoData.")
}

group_col <- available[1]
selected_groups <- pick_top_two_groups(pheno[[group_col]])

n_per_group <- 4L
idx_control <- which(as.character(pheno[[group_col]]) == selected_groups[1])[1:n_per_group]
idx_experiment <- which(as.character(pheno[[group_col]]) == selected_groups[2])[1:n_per_group]
idx <- c(idx_control, idx_experiment)
idx <- idx[!is.na(idx)]

if (length(idx) < 4L) {
    stop("Not enough samples to build the microarray contrast template.")
}

eset_sub <- eset[, idx]
pheno_sub <- Biobase::pData(eset_sub)
pheno_sub$adam_group <- as.character(pheno_sub[[group_col]])
Biobase::pData(eset_sub) <- pheno_sub

expr <- Biobase::exprs(eset_sub)
own_db <- build_data_driven_db(expr, n_groups = 60L, prefix = "MICRO")

result <- ADAM::run_adam(
    x = eset_sub,
    mode = "partial",
    group_col = "adam_group",
    contrast = selected_groups,
    MinGene = 15L,
    MaxGene = 600L,
    DBSpecies = own_db,
    AnalysisDomain = "own",
    GeneIdentifier = "probe"
)

dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file = "outputs/models/microarray_all_partial.rds")
utils::write.csv(
    as.data.frame(SummarizedExperiment::rowData(result)),
    file = "outputs/models/microarray_all_partial_table.csv",
    row.names = FALSE
)

message("Microarray model completed: outputs/models/microarray_all_partial.rds")
