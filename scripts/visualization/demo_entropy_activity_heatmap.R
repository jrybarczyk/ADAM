#!/usr/bin/env Rscript

# End-to-end demo: run ADAM (partial) and render entropy/activity heatmap.

if (!requireNamespace("ADAM", quietly = TRUE)) {
    stop("ADAM package is required. Install/load package first.")
}

suppressPackageStartupMessages(library(ADAM))

data("ExpressionAedes", package = "ADAM")
data("KeggPathwaysAedes", package = "ADAM")

res <- suppressMessages(ADAM::run_adam(
    x = ExpressionAedes,
    mode = "partial",
    comparison_id = c("control1,experiment1", "control2,experiment2"),
    MinGene = 3L,
    MaxGene = 20L,
    DBSpecies = KeggPathwaysAedes,
    AnalysisDomain = "own",
    GeneIdentifier = "geneStableID"
))

dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(res, "outputs/models/demo_aedes_partial_for_heatmap.rds")

cmd <- c(
    "scripts/visualization/plot_entropy_activity_heatmap.R",
    "--input", "outputs/models/demo_aedes_partial_for_heatmap.rds",
    "--output", "outputs/models/demo_aedes_entropy_activity_heatmap.png",
    "--comparison", "control1,experiment1",
    "--top", "40"
)

status <- system2(command = "Rscript", args = cmd)
if (status != 0) {
    stop("Heatmap script failed.")
}

message("Demo completed: outputs/models/demo_aedes_entropy_activity_heatmap.png")
