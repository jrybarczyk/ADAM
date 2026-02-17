test_that("checkExpressionData validates core input rules", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)

    checked <- ADAM:::.checkExpressionData(ExpressionAedes)
    expect_true(is.data.frame(checked))
    expect_identical(colnames(checked)[1], "gene")
    expect_true(all(vapply(checked[, -1], is.numeric, logical(1))))

    bad_df <- data.frame(gene = c("a", "b"), sample1 = c(1, 2))
    expect_error(
        ADAM:::.checkExpressionData(bad_df),
        "at list 3 columns"
    )

    expect_error(
        ADAM:::.checkExpressionData("this_file_does_not_exist.tsv"),
        "valid file or a"
    )

    dup_df <- data.frame(
        gene = c("g1", "g1", "g2"),
        c1 = c(1, 2, 3),
        c2 = c(4, 5, 6)
    )
    expect_error(ADAM:::.checkExpressionData(dup_df), "must be unique")

    non_numeric_df <- data.frame(
        gene = c("g1", "g2", "g3"),
        c1 = c("a", "b", "c"),
        c2 = c(1, 2, 3)
    )
    expect_error(ADAM:::.checkExpressionData(non_numeric_df), "must be numeric")

    inf_df <- data.frame(
        gene = c("g1", "g2", "g3"),
        c1 = c(1, Inf, 3),
        c2 = c(4, 5, 6)
    )
    expect_error(ADAM:::.checkExpressionData(inf_df), "must be finite")

    na_df <- data.frame(
        gene = c("g1", "g2", "g3"),
        c1 = c(1, NA, 3),
        c2 = c(4, 5, 6)
    )
    expect_error(ADAM:::.checkExpressionData(na_df), "must be finite")

    nan_df <- data.frame(
        gene = c("g1", "g2", "g3"),
        c1 = c(1, NaN, 3),
        c2 = c(4, 5, 6)
    )
    expect_error(ADAM:::.checkExpressionData(nan_df), "must be finite")

    dup_sample_df <- data.frame(
        gene = c("g1", "g2", "g3"),
        s1 = c(1, 2, 3),
        s1 = c(4, 5, 6),
        check.names = FALSE
    )
    expect_error(ADAM:::.checkExpressionData(dup_sample_df), "Sample names")

    tf <- tempfile(fileext = ".tsv")
    on.exit(unlink(tf), add = TRUE)
    utils::write.table(ExpressionAedes, file = tf, sep = "\t",
        row.names = FALSE, quote = FALSE)
    checked_file <- ADAM:::.checkExpressionData(tf)
    expect_true(is.data.frame(checked_file))
})

test_that("comparison and scalar checks behave as expected", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)

    expect_identical(
        ADAM:::.checkComparisonID("control1,experiment1", ExpressionAedes),
        "control1,experiment1"
    )
    expect_error(
        ADAM:::.checkComparisonID("control1", ExpressionAedes),
        "must have 2 elements"
    )
    expect_error(
        ADAM:::.checkComparisonID("control1,invalid_sample", ExpressionAedes),
        "must correspond to the"
    )
    expect_error(
        ADAM:::.checkComparisonID(data.frame(x = "control1,experiment1"), ExpressionAedes),
        "must be a vector"
    )

    expect_identical(ADAM:::.checkGeneNumbers(3, 10), c(3L, 10L))
    expect_identical(ADAM:::.checkGeneNumbers(3.2, 10.8), c(3L, 10L))
    expect_error(ADAM:::.checkGeneNumbers(11, 2), "greater than the minimum")
    expect_error(ADAM:::.checkGeneNumbers("a", 2), "must be integer")

    expect_identical(ADAM:::.checkSeedNumber(1), 1)
    expect_error(ADAM:::.checkSeedNumber("x"), "must.*numeric|numeric")
    expect_error(ADAM:::.checkSeedNumber(-1), "must be a positive")

    expect_identical(ADAM:::.checkBootstrapNumber(1000), 1000L)
    expect_identical(ADAM:::.checkBootstrapNumber(10.2), 10L)
    expect_error(ADAM:::.checkBootstrapNumber("x"), "must be integer")
    expect_error(ADAM:::.checkBootstrapNumber(0), "must be positive")

    expect_identical(ADAM:::.checkPCorrection(0.05), 0.05)
    expect_error(ADAM:::.checkPCorrection("x"), "must be numeric")
    expect_error(ADAM:::.checkPCorrection(1.5), "between zero and one")

    expect_identical(ADAM:::.checkPCorrectionMethod("bh"), "BH")
    expect_error(ADAM:::.checkPCorrectionMethod("abc"), "not found")

    expect_identical(ADAM:::.checkWilcoxonTest(TRUE), TRUE)
    expect_error(ADAM:::.checkWilcoxonTest("yes"), "must be logical")

    expect_identical(ADAM:::.checkFisherTest(FALSE), FALSE)
    expect_error(ADAM:::.checkFisherTest("no"), "be logical")
})

test_that("species and package helpers work", {
    species <- ADAM:::.CreateDataSpecies()
    expect_true(is.data.frame(species))
    expect_true("org.Hs.eg.db" %in% species$DBSpeciesList)

    expect_true(ADAM:::.checkPackage("stats"))
    expect_error(
        ADAM:::.checkPackage("definitely.not.a.real.package"),
        "not available"
    )

    expect_length(ADAM:::.adam_get_ontology_function("org.At.tair.db"), 4)
    expect_length(ADAM:::.adam_get_ontology_function("org.Unknown.eg.db"), 4)
})

test_that("own-domain path for analysis/domain/species/geneidentifier works", {
    expect_warning(data("ExpressionAedes", package = "ADAM"), NA)
    expect_warning(data("KeggPathwaysAedes", package = "ADAM"), NA)

    result <- ADAM:::.checkAnalysisDomain_DBSpecies_GeneIdentifier(
        AnalysisDomain = "own",
        DBSpecies = KeggPathwaysAedes,
        GeneIdentifier = "geneStableID",
        GeneExpressionData = ExpressionAedes
    )

    expect_identical(result[[1]], "own")
    expect_true(is.list(result[[2]]))
    expect_true(is.list(result[[3]]))
    expect_identical(result[[4]], "genestableid")

    expect_error(
        ADAM:::.checkAnalysisDomain_DBSpecies_GeneIdentifier(
            AnalysisDomain = "invalid",
            DBSpecies = KeggPathwaysAedes,
            GeneIdentifier = "geneStableID",
            GeneExpressionData = ExpressionAedes
        ),
        "does not exist"
    )
})

test_that("database and conversion helpers cover success and failure branches", {
    expr <- data.frame(
        gene = c("g1", "g2", "g3"),
        s1 = c(1, 2, 3),
        s2 = c(3, 2, 1)
    )
    db_ok <- data.frame(
        gene = c("g1", "g2"),
        ID = c("id1", "id1"),
        Description = c("d1", "d1")
    )

    checked <- ADAM:::.checkDB(db_ok, expr)
    expect_true(is.list(checked))
    expect_length(checked, 2)

    expect_error(ADAM:::.checkDB(db_ok[, c("gene", "ID")], expr), "must have 3 columns")
    expect_error(ADAM:::.checkDB(db_ok[0, ], expr), "is empty")

    db_nomatch <- data.frame(
        gene = c("x1", "x2"),
        ID = c("id1", "id2"),
        Description = c("d1", "d2")
    )
    expect_error(ADAM:::.checkDB(db_nomatch, expr), "not in the species")

    db_same_cols <- data.frame(
        gene = c("g1", "g2"),
        s1 = c("id1", "id2"),
        s2 = c("d1", "d2")
    )
    expect_error(ADAM:::.checkDB(db_same_cols, expr), "column names")

    expect_identical(
        ADAM:::.CheckTerm(list("a", "b"), "go:1", c("go:1", "go:2")),
        list("a", "b")
    )
    expect_null(ADAM:::.CheckTerm(list("a"), "go:9", c("go:1", "go:2")))

    gene_map <- data.frame(
        entrez = c("e1", "e2"),
        symbol = c("s1", "s2")
    )
    converted <- ADAM:::.GeneIDConversion(
        ListElement = c("e1"),
        ListGenes = gene_map,
        StandardGeneID = "entrez",
        AlternativeGeneID = "symbol"
    )
    expect_identical(converted, "s1")
})

test_that("domain-list helper validates empty and non-empty mappings", {
    expr <- data.frame(
        gene = c("g1", "g2", "g3"),
        s1 = c(1, 2, 3),
        s2 = c(3, 2, 1)
    )
    term_ont <- data.frame(
        ID = c("id1", "id2"),
        Description = c("desc1", "desc2"),
        stringsAsFactors = FALSE
    )

    lists_ok <- list(id1 = c("g1", "g2"), id2 = c("g9"))
    result_ok <- ADAM:::.adam_build_domain_lists(
        DBFunctionsList = lists_ok,
        TermOntologies = term_ont,
        ExpressionData = expr,
        AnalysisDomain = "gobp"
    )
    expect_true(is.list(result_ok))
    expect_length(result_ok$DBFunctionsListSample, 1)

    lists_empty <- list(id2 = c("g9"))
    expect_error(
        ADAM:::.adam_build_domain_lists(
            DBFunctionsList = lists_empty,
            TermOntologies = term_ont,
            ExpressionData = expr,
            AnalysisDomain = "gobp"
        ),
        "not related to any gobp term"
    )
})
