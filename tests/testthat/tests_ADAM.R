### Tests
test_that("data size from ADAnalysis is ok", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)
    expect_equal(length(as.data.frame(suppressMessages(
            ADAnalysis(ComparisonID = c("control1,experiment1"),
            ExpressionData = ExpressionAedes, MinGene = 3,
            DBSpecies = KeggPathwaysAedes, MaxGene = 1000,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"))[[2]])$ID),
            length(unique(KeggPathwaysAedes$pathwayID)))
})


test_that("data size from GFAGAnalysis is ok", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)
    expect_equal(length(as.data.frame(suppressMessages(
            GFAGAnalysis(ComparisonID  = c("control1,experiment1"),
            ExpressionData = ExpressionAedes, MinGene = 3, MaxGene = 20, 
            SeedNumber = 1049, BootstrapNumber = 1000, PCorrection = 0.05,
            DBSpecies = KeggPathwaysAedes, PCorrectionMethod = "fdr", 
            WilcoxonTest = FALSE, FisherTest = FALSE,
            AnalysisDomain = "own",GeneIdentifier="geneStableID"))[[2]])$ID),
            length(unique(KeggPathwaysAedes$pathwayID)))
})

test_that("run_adam works with SingleCellExperiment input", {
    skip_if_not_installed("SingleCellExperiment")

    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    expr <- as.matrix(ExpressionAedes[, -1])
    rownames(expr) <- as.character(ExpressionAedes$gene)

    col_data <- S4Vectors::DataFrame(
        condition = c("control", "experiment", "control", "experiment"),
        sample_id = colnames(expr)
    )
    rownames(col_data) <- colnames(expr)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = expr),
        colData = col_data
    )

    result <- suppressMessages(run_adam(
        x = sce,
        mode = "partial",
        assay_name = "counts",
        group_col = "condition",
        contrast = c("control", "experiment"),
        DBSpecies = KeggPathwaysAedes,
        AnalysisDomain = "own",
        GeneIdentifier = "geneStableID"
    ))

    expect_true(methods::is(result, "SummarizedExperiment"))
    expect_true("metrics" %in% names(SummarizedExperiment::assays(result)))
    expect_identical(S4Vectors::metadata(result)$mode, "partial")
})

test_that("run_adam works with SummarizedExperiment group_col and contrast", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

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

    result <- suppressMessages(run_adam(
        x = se,
        mode = "partial",
        assay_name = "counts",
        group_col = "condition",
        contrast = c("control", "experiment"),
        DBSpecies = KeggPathwaysAedes,
        AnalysisDomain = "own",
        GeneIdentifier = "geneStableID"
    ))

    expect_true(methods::is(result, "SummarizedExperiment"))
    expect_true(nrow(result) > 0)
    expect_identical(S4Vectors::metadata(result)$mode, "partial")
})

test_that("run_adam supports pseudo-bulk aggregation with aggregate_by", {
    skip_if_not_installed("SingleCellExperiment")

    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    expr <- as.matrix(ExpressionAedes[, -1])
    rownames(expr) <- as.character(ExpressionAedes$gene)

    # Duplicate cells to emulate a small single-cell-like layout.
    expr <- cbind(expr, expr)
    colnames(expr) <- paste0("cell_", seq_len(ncol(expr)))

    sample_id <- rep(c("S1", "S2", "S3", "S4"), 2)
    cluster <- rep(c("A", "A", "B", "B"), 2)
    condition <- ifelse(sample_id %in% c("S1", "S3"), "control", "experiment")

    col_data <- S4Vectors::DataFrame(
        sample_id = sample_id,
        cluster = cluster,
        condition = condition,
        row.names = colnames(expr)
    )

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = expr),
        colData = col_data
    )

    result <- suppressMessages(run_adam(
        x = sce,
        mode = "partial",
        assay_name = "counts",
        group_col = "condition",
        contrast = c("control", "experiment"),
        aggregate_by = c("sample_id", "cluster", "condition"),
        aggregate_fun = "sum",
        DBSpecies = KeggPathwaysAedes,
        AnalysisDomain = "own",
        GeneIdentifier = "geneStableID"
    ))

    expect_true(methods::is(result, "SummarizedExperiment"))
    expect_true(nrow(result) > 0)
    expect_identical(S4Vectors::metadata(result)$mode, "partial")
})
