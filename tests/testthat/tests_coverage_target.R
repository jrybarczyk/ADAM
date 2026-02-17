test_that("checkExpressionData handles SummarizedExperiment inputs", {
    mat <- matrix(
        c(1, 2, 3, 4, 5, 6, 7, 8, 9),
        nrow = 3,
        dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3"))
    )
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))

    checked <- ADAM:::.checkExpressionData(se)
    expect_true(is.data.frame(checked))
    expect_identical(colnames(checked)[1], "gene")
    expect_identical(as.character(checked$gene), rownames(mat))
    expect_identical(as.numeric(checked$s1), as.numeric(mat[, "s1"]))
})

test_that("run_adam and helpers cover metadata/comparison branches", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    expr <- as.matrix(ExpressionAedes[, -1])
    rownames(expr) <- as.character(ExpressionAedes$gene)

    col_data <- S4Vectors::DataFrame(
        condition = c("control", "experiment", "other", "other"),
        batch = c("b1", "b1", "b2", "b2"),
        row.names = colnames(expr)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = expr),
        colData = col_data
    )

    expect_error(
        ADAM::run_adam(
            x = ExpressionAedes,
            mode = "partial",
            DBSpecies = KeggPathwaysAedes,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"
        ),
        "provide comparison_id explicitly"
    )

    expect_error(
        ADAM::run_adam(
            x = 123,
            mode = "partial",
            DBSpecies = KeggPathwaysAedes,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"
        ),
        "Unsupported input type"
    )

    expect_error(
        ADAM:::.adam_build_comparisons(
            sample_names = colnames(expr),
            sample_data = as.data.frame(col_data),
            group_col = "missing",
            contrast = NULL,
            reference = NULL
        ),
        "Provide a valid group_col"
    )

    expect_error(
        ADAM:::.adam_build_comparisons(
            sample_names = colnames(expr),
            sample_data = as.data.frame(col_data),
            group_col = "condition",
            contrast = NULL,
            reference = "not_here"
        ),
        "reference group"
    )

    expect_error(
        ADAM:::.adam_build_comparisons(
            sample_names = colnames(expr),
            sample_data = as.data.frame(col_data),
            group_col = "condition",
            contrast = c("control"),
            reference = NULL
        ),
        "contrast must have length 2"
    )

    expect_error(
        ADAM:::.adam_build_comparisons(
            sample_names = colnames(expr),
            sample_data = as.data.frame(col_data),
            group_col = "condition",
            contrast = c("absent", "control"),
            reference = NULL
        ),
        "No samples found"
    )

    two_group_col_data <- S4Vectors::DataFrame(
        condition = c("control", "experiment", "control", "experiment"),
        row.names = colnames(expr)
    )
    comps_ref <- ADAM:::.adam_build_comparisons(
        sample_names = colnames(expr),
        sample_data = as.data.frame(two_group_col_data),
        group_col = "condition",
        contrast = NULL,
        reference = "control"
    )
    expect_true(length(comps_ref) > 0)

    agg_mean <- ADAM:::.adam_aggregate_expression(
        expr_mat = expr,
        sample_data = data.frame(
            sample = c("s1", "s2", "s1", "s2"),
            condition = c("control", "experiment", "control", "experiment"),
            stringsAsFactors = FALSE
        ),
        aggregate_by = c("sample", "condition"),
        aggregate_fun = "mean"
    )
    expect_true(is.matrix(agg_mean$expr_mat))
    expect_true(is.data.frame(agg_mean$sample_data))

    expect_error(
        ADAM:::.adam_aggregate_expression(
            expr_mat = expr,
            sample_data = as.data.frame(two_group_col_data),
            aggregate_by = c("missing"),
            aggregate_fun = "sum"
        ),
        "aggregate_by fields"
    )

    rd <- S4Vectors::DataFrame(gid = as.character(ExpressionAedes$gene))
    se_with_rowdata <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = expr),
        rowData = rd,
        colData = two_group_col_data
    )
    expect_error(
        ADAM::run_adam(
            x = se_with_rowdata,
            mode = "partial",
            assay_name = "counts",
            group_col = "condition",
            contrast = c("control", "experiment"),
            gene_id_source = "missing_col",
            DBSpecies = KeggPathwaysAedes,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"
        ),
        "gene_id_source"
    )

    se_no_gene_ids <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = expr),
        rowData = S4Vectors::DataFrame(row.names = rownames(expr)),
        colData = two_group_col_data
    )
    rownames(se_no_gene_ids) <- NULL
    expect_error(
        ADAM::run_adam(
            x = se_no_gene_ids,
            mode = "partial",
            assay_name = "counts",
            group_col = "condition",
            contrast = c("control", "experiment"),
            DBSpecies = KeggPathwaysAedes,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"
        ),
        "Gene IDs are missing"
    )
})

test_that("run_adam complete mode and ExpressionSet path work", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    result_complete <- suppressMessages(ADAM::run_adam(
        x = ExpressionAedes,
        mode = "complete",
        comparison_id = c("control1,experiment1"),
        MinGene = 3,
        MaxGene = 20,
        SeedNumber = 1049,
        BootstrapNumber = 10,
        PCorrection = 0.05,
        DBSpecies = KeggPathwaysAedes,
        PCorrectionMethod = "fdr",
        WilcoxonTest = FALSE,
        FisherTest = FALSE,
        AnalysisDomain = "own",
        GeneIdentifier = "geneStableID"
    ))
    expect_true(methods::is(result_complete, "SummarizedExperiment"))
    expect_identical(S4Vectors::metadata(result_complete)$mode, "complete")

    skip_if_not_installed("Biobase")
    eset <- Biobase::ExpressionSet(
        assayData = as.matrix(ExpressionAedes[, -1])
    )
    Biobase::featureNames(eset) <- as.character(ExpressionAedes$gene)
    Biobase::pData(eset)$condition <- c("control", "experiment", "control", "experiment")
    rownames(Biobase::pData(eset)) <- Biobase::sampleNames(eset)

    result_eset <- suppressMessages(ADAM::run_adam(
        x = eset,
        mode = "partial",
        group_col = "condition",
        contrast = c("control", "experiment"),
        DBSpecies = KeggPathwaysAedes,
        AnalysisDomain = "own",
        GeneIdentifier = "geneStableID"
    ))
    expect_true(methods::is(result_eset, "SummarizedExperiment"))
})

test_that("conversion helper handles heterogeneous legacy tables", {
    empty_legacy <- list(
        data.frame(gene = character(), GroupID = character(), stringsAsFactors = FALSE),
        list(data.frame(ID = character(), Description = character(), stringsAsFactors = FALSE))
    )
    se_empty <- ADAM:::.adam_to_summarized_experiment(
        legacy_result = empty_legacy,
        comparison_id = "c1,e1",
        mode = "partial"
    )
    expect_true(methods::is(se_empty, "SummarizedExperiment"))
    expect_identical(nrow(se_empty), 0L)

    legacy <- list(
        data.frame(gene = c("g1", "g2"), GroupID = c("id1", "id2")),
        list(
            data.frame(ID = "id1", Description = "d1", n = 0.5),
            data.frame(ID = "id2", extra = "x", n = 0.3)
        )
    )
    se <- ADAM:::.adam_to_summarized_experiment(
        legacy_result = legacy,
        comparison_id = c("c1,e1", "c2,e2"),
        mode = "partial"
    )
    expect_true(methods::is(se, "SummarizedExperiment"))
    expect_true(nrow(se) == 2)
})

test_that("ECGMainData defaults and bootstrap finalization branches are exercised", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    obj <- suppressMessages(ADAM:::ECGMainData(
        ComparisonID = c("control1,experiment1"),
        ExpressionData = ExpressionAedes,
        MinGene = NULL,
        MaxGene = NULL,
        SeedNumber = NULL,
        BootstrapNumber = NULL,
        PCorrection = NULL,
        DBSpecies = KeggPathwaysAedes,
        PCorrectionMethod = NULL,
        WilcoxonTest = NULL,
        FisherTest = NULL,
        AnalysisDomain = "own",
        GeneIdentifier = "geneStableID",
        completeTest = TRUE
    ))
    expect_identical(obj@MinGene, 3L)
    expect_identical(obj@MaxGene, 2000L)
    expect_identical(obj@SeedNumber, 1049)
    expect_identical(obj@BootstrapNumber, 1000L)

    custom_obj <- methods::new(
        Class = "ECGMainData",
        ComparisonID = c("c1,e1"),
        ExpressionData = data.frame(gene = c("g1", "g2"), c1 = c(1, 2), e1 = c(2, 3)),
        MinGene = 1L,
        MaxGene = 10L,
        SeedNumber = 1,
        BootstrapNumber = 10L,
        PCorrection = 0.5,
        DBSpeciesFunctionsSample = list(),
        DBSpeciesFunctionsRaw = list("id1 <==> desc1" = 2),
        PCorrectionMethod = "fdr",
        WilcoxonTest = TRUE,
        FisherTest = TRUE,
        AnalysisDomain = "own",
        GeneIdentifier = "gene"
    )

    result_boot <- data.frame(
        ID = "id1",
        Description = "desc1",
        Sample_Number_Genes = 2,
        H_c1 = 0.4,
        H_e1 = 0.6,
        N_c1 = 2,
        N_e1 = 3,
        h = 0.6,
        n = 0.6,
        pValue_h = 0.05,
        pValue_n = 0.08,
        stringsAsFactors = FALSE
    )
    groups <- list(data.frame(gene = c("g1", "g2"), control = c(1, 2), experiment = c(2, 3)))
    comparison_df <- data.frame(
        gene = c("g1", "g2", "g3"),
        control = c(1, 2, 1),
        experiment = c(2, 3, 0),
        stringsAsFactors = FALSE
    )

    finalized <- suppressMessages(ADAM:::.adam_finalize_bootstrap(
        ResultBootstrap = result_boot,
        ECGObject = custom_obj,
        GroupComparisons = groups,
        Comparison = comparison_df,
        compID = c("c1", "e1")
    ))
    expect_true("Wilcox_pvalue" %in% colnames(finalized))
    expect_true("Fisher_pvalue" %in% colnames(finalized))

    zero_stats <- ADAM:::SampleStatistics(
        GroupComparisons = data.frame(gene = c("g1", "g2"), control = c(0, 0), experiment = c(0, 0)),
        Control = "c1",
        Experiment = "e1"
    )
    expect_identical(zero_stats$n, 0)
    expect_identical(ADAM:::Diversity(data.frame(a = c(0), b = c(0)), 1), 0)

    expect_identical(
        ADAM:::WCAnalysis(data.frame(gene = c("g1", "g2"), c1 = c(1, 1), e1 = c(1, 1))),
        1
    )
    expect_identical(
        ADAM:::WCAnalysis(data.frame(gene = "g1", c1 = 1, e1 = 1)),
        1
    )

    fisher_degenerate <- ADAM:::FSHAnalysis(
        pValuesGroups = data.frame(
            gene = c("g1", "g2"),
            control = c(0, 0),
            experiment = c(0, 0)
        ),
        Comparison = data.frame(
            gene = c("g1", "g2", "g3"),
            control = c(0, 0, 0),
            experiment = c(0, 0, 0)
        )
    )
    expect_identical(fisher_degenerate, 1)
})

test_that("comparison and raw gene lookup use exact ID/sample matching", {
    expr <- data.frame(
        gene = c("g1", "g2"),
        control1 = c(1, 2),
        control10 = c(100, 200),
        experiment1 = c(3, 4),
        stringsAsFactors = FALSE
    )

    obj <- methods::new(
        Class = "ECGMainData",
        ComparisonID = c("control1,experiment1"),
        ExpressionData = expr,
        MinGene = 1L,
        MaxGene = 10L,
        SeedNumber = 1,
        BootstrapNumber = 10L,
        PCorrection = 0.5,
        DBSpeciesFunctionsSample = list("id1 <==> desc1" = c("g1", "g2")),
        DBSpeciesFunctionsRaw = list("id1" = 2, "id10" = 999),
        PCorrectionMethod = "fdr",
        WilcoxonTest = FALSE,
        FisherTest = FALSE,
        AnalysisDomain = "own",
        GeneIdentifier = "gene"
    )

    parsed <- ADAM:::.adam_build_comparison("control1,experiment1", obj)
    expect_identical(parsed$Comparison$control, c(1, 2))
    expect_identical(parsed$Comparison$experiment, c(3, 4))

    gs <- list("id1 <==> desc1" = data.frame(
        Sample_Number_Genes = 2,
        H_control1 = 0.2,
        H_experiment1 = 0.8,
        N_control1 = 2,
        N_experiment1 = 3,
        h = 0.8,
        n = 0.6
    ))
    partial <- ADAM:::.adam_partial_analysis_result(gs, obj, c("control1", "experiment1"))
    expect_identical(partial$Raw_Number_Genes, 2)
})

test_that("bootstrap complete mode is reproducible with fixed seed", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    run_once <- function() {
        suppressMessages(ADAM::GFAGAnalysis(
            ComparisonID = c("control1,experiment1"),
            ExpressionData = ExpressionAedes,
            MinGene = 3,
            MaxGene = 20,
            SeedNumber = 1049,
            BootstrapNumber = 50,
            PCorrection = 0.05,
            DBSpecies = KeggPathwaysAedes,
            PCorrectionMethod = "fdr",
            WilcoxonTest = FALSE,
            FisherTest = FALSE,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"
        ))[[2]][[1]]
    }

    r1 <- run_once()
    r2 <- run_once()

    expect_equal(r1$ID, r2$ID)
    expect_equal(as.numeric(r1$pValue_h), as.numeric(r2$pValue_h), tolerance = 1e-12)
    expect_equal(as.numeric(r1$pValue_n), as.numeric(r2$pValue_n), tolerance = 1e-12)
    expect_equal(as.numeric(r1$qValue_h), as.numeric(r2$qValue_h), tolerance = 1e-12)
    expect_equal(as.numeric(r1$qValue_n), as.numeric(r2$qValue_n), tolerance = 1e-12)
})
