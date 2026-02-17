#!/usr/bin/env Rscript

source("scripts/models/common_model_utils.R")

ensure_packages(c("ADAM", "airway", "SummarizedExperiment", "S4Vectors"), bioc = TRUE)

library(ADAM)
library(airway)
library(SummarizedExperiment)

se <- airway::airway
counts <- SummarizedExperiment::assay(se, "counts")

keep_genes <- rowSums(counts) > 10
se <- se[keep_genes, ]
counts <- counts[keep_genes, , drop = FALSE]

sample_group <- as.character(SummarizedExperiment::colData(se)$dex)
if (!all(c("untrt", "trt") %in% unique(sample_group))) {
    stop("Expected airway dex groups 'untrt' and 'trt'.")
}

comparison_id <- make_pairwise_comparisons(
    control_samples = colnames(se)[sample_group == "untrt"],
    experiment_samples = colnames(se)[sample_group == "trt"],
    max_pairs = 4L
)

own_db <- build_data_driven_db(counts, n_groups = 50L, prefix = "AIRWAY")

result <- ADAM::run_adam(
    x = se,
    mode = "partial",
    assay_name = "counts",
    comparison_id = comparison_id,
    MinGene = 10L,
    MaxGene = 400L,
    DBSpecies = own_db,
    AnalysisDomain = "own",
    GeneIdentifier = "ensembl"
)

dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file = "outputs/models/rnaseq_airway_partial.rds")
utils::write.csv(
    as.data.frame(SummarizedExperiment::rowData(result)),
    file = "outputs/models/rnaseq_airway_partial_table.csv",
    row.names = FALSE
)

message("RNA-seq model completed: outputs/models/rnaseq_airway_partial.rds")
