#' @title Run ADAM with Bioconductor Containers
#' @description
#' Wrapper interface to run ADAM using Bioconductor data containers
#' (`SummarizedExperiment`, `SingleCellExperiment`, or `ExpressionSet`) or
#' legacy inputs (`data.frame` / file path). This function builds ADAM
#' comparisons from sample metadata and returns results as a
#' `SummarizedExperiment`.
#' @param x Input expression container. Supported:
#' `SummarizedExperiment`, `SingleCellExperiment`, `ExpressionSet`,
#' `data.frame`, or a tab-delimited file path.
#' @param mode Analysis mode: `"complete"` (GFAGAnalysis) or `"partial"`
#' (ADAnalysis).
#' @param assay_name Assay name for `SummarizedExperiment` /
#' `SingleCellExperiment`. If `NULL`, the first assay is used.
#' @param group_col Column in sample metadata used to define control/experiment
#' groups when `comparison_id` is not provided.
#' @param contrast Character vector of length 2 with group labels in the form
#' `c(control, experiment)`.
#' @param reference Optional control group label used when `contrast` is not
#' provided and exactly two groups are present.
#' @param comparison_id Optional explicit ADAM comparison vector in the format
#' `"sampleA,sampleB"`.
#' @param gene_id_source Optional rowData column to use as gene ID. If `NULL`,
#' row names are used.
#' @param aggregate_by Optional sample metadata columns for pseudo-bulk
#' aggregation (e.g. `c("sample_id", "cluster")`) before analysis.
#' @param aggregate_fun Aggregation function for pseudo-bulk: `"sum"` or
#' `"mean"`.
#' @param MinGene,MaxGene,SeedNumber,BootstrapNumber,PCorrection,
#' DBSpecies,PCorrectionMethod,WilcoxonTest,FisherTest,AnalysisDomain,
#' GeneIdentifier Same arguments used by `GFAGAnalysis` / `ADAnalysis`.
#' @return A `SummarizedExperiment` with one assay named `metrics`, row-level
#' annotations for GFAG outputs, and `metadata()` containing:
#' `legacy_result`, `gene_function_map`, `comparison_id`, and `mode`.
#' @export
run_adam <- function(
    x,
    mode = c("complete", "partial"),
    assay_name = NULL,
    group_col = NULL,
    contrast = NULL,
    reference = NULL,
    comparison_id = NULL,
    gene_id_source = NULL,
    aggregate_by = NULL,
    aggregate_fun = c("sum", "mean"),
    MinGene = 3,
    MaxGene = 2000,
    SeedNumber = 1049,
    BootstrapNumber = 1000,
    PCorrection = 0.05,
    DBSpecies = NULL,
    PCorrectionMethod = "fdr",
    WilcoxonTest = FALSE,
    FisherTest = FALSE,
    AnalysisDomain = NULL,
    GeneIdentifier = NULL
) {
    mode <- match.arg(mode)
    aggregate_fun <- match.arg(aggregate_fun)

    normalized <- .adam_normalize_input(
        x = x,
        assay_name = assay_name,
        group_col = group_col,
        contrast = contrast,
        reference = reference,
        comparison_id = comparison_id,
        gene_id_source = gene_id_source,
        aggregate_by = aggregate_by,
        aggregate_fun = aggregate_fun
    )

    if (mode == "complete") {
        legacy <- GFAGAnalysis(
            ComparisonID = normalized$comparison_id,
            ExpressionData = normalized$expression_data,
            MinGene = MinGene,
            MaxGene = MaxGene,
            SeedNumber = SeedNumber,
            BootstrapNumber = BootstrapNumber,
            PCorrection = PCorrection,
            DBSpecies = DBSpecies,
            PCorrectionMethod = PCorrectionMethod,
            WilcoxonTest = WilcoxonTest,
            FisherTest = FisherTest,
            AnalysisDomain = AnalysisDomain,
            GeneIdentifier = GeneIdentifier
        )
    } else {
        legacy <- ADAnalysis(
            ComparisonID = normalized$comparison_id,
            ExpressionData = normalized$expression_data,
            MinGene = MinGene,
            MaxGene = MaxGene,
            DBSpecies = DBSpecies,
            AnalysisDomain = AnalysisDomain,
            GeneIdentifier = GeneIdentifier
        )
    }

    .adam_to_summarized_experiment(
        legacy_result = legacy,
        comparison_id = normalized$comparison_id,
        mode = mode
    )
}

.adam_normalize_input <- function(
    x,
    assay_name,
    group_col,
    contrast,
    reference,
    comparison_id,
    gene_id_source,
    aggregate_by,
    aggregate_fun
) {
    if (is.data.frame(x) || is.character(x)) {
        if (is.null(comparison_id)) {
            stop("For data.frame/file inputs, provide comparison_id explicitly.")
        }
        return(list(expression_data = x, comparison_id = comparison_id))
    }

    if (methods::is(x, "ExpressionSet")) {
        biobase_pkg <- "Biobase"
        if (!requireNamespace(biobase_pkg, quietly = TRUE)) {
            stop("Biobase is required to use ExpressionSet input.")
        }
        expr_mat <- getExportedValue(biobase_pkg, "exprs")(x)
        sample_data <- getExportedValue(biobase_pkg, "pData")(x)
        gene_ids <- getExportedValue(biobase_pkg, "featureNames")(x)
        sample_names <- colnames(expr_mat)
    } else if (methods::is(x, "SummarizedExperiment") ||
                methods::is(x, "SingleCellExperiment")) {
        if (is.null(assay_name)) {
            assay_name <- names(SummarizedExperiment::assays(x))[1]
        }
        expr_mat <- SummarizedExperiment::assay(x, assay_name)
        sample_data <- as.data.frame(SummarizedExperiment::colData(x))
        sample_names <- colnames(expr_mat)

        if (!is.null(gene_id_source)) {
            if (!gene_id_source %in% colnames(SummarizedExperiment::rowData(x))) {
                stop(sprintf("gene_id_source %s not found in rowData.", gene_id_source))
            }
            gene_ids <- as.character(SummarizedExperiment::rowData(x)[[gene_id_source]])
        } else {
            gene_ids <- rownames(x)
        }
    } else {
        stop("Unsupported input type for run_adam().")
    }

    if (is.null(gene_ids)) {
        stop("Gene IDs are missing. Provide row names or gene_id_source.")
    }

    if (!is.null(aggregate_by)) {
        aggregated <- .adam_aggregate_expression(
            expr_mat = expr_mat,
            sample_data = sample_data,
            aggregate_by = aggregate_by,
            aggregate_fun = aggregate_fun
        )
        expr_mat <- aggregated$expr_mat
        sample_data <- aggregated$sample_data
        sample_names <- colnames(expr_mat)
    }

    expression_data <- as.data.frame(expr_mat, stringsAsFactors = FALSE)
    expression_data <- cbind.data.frame(gene = as.character(gene_ids),
                                        expression_data,
                                        stringsAsFactors = FALSE)

    if (is.null(comparison_id)) {
        comparison_id <- .adam_build_comparisons(
            sample_names = sample_names,
            sample_data = sample_data,
            group_col = group_col,
            contrast = contrast,
            reference = reference
        )
    }

    list(expression_data = expression_data, comparison_id = comparison_id)
}

.adam_build_comparisons <- function(
    sample_names,
    sample_data,
    group_col,
    contrast,
    reference
) {
    if (is.null(group_col) || !group_col %in% colnames(sample_data)) {
        stop("Provide a valid group_col in sample metadata.")
    }

    groups <- as.character(sample_data[[group_col]])
    unique_groups <- unique(groups)

    if (is.null(contrast)) {
        if (!is.null(reference)) {
            if (!reference %in% unique_groups) {
                stop(sprintf("reference group %s not found.", reference))
            }
            other_groups <- setdiff(unique_groups, reference)
            if (length(other_groups) != 1) {
                stop("With reference, metadata must define exactly two groups.")
            }
            contrast <- c(reference, other_groups[1])
        } else {
            if (length(unique_groups) != 2) {
                stop("Provide contrast when metadata contains more than two groups.")
            }
            contrast <- unique_groups
        }
    }

    if (length(contrast) != 2) {
        stop("contrast must have length 2: c(control, experiment).")
    }

    control_samples <- sample_names[groups == contrast[1]]
    experiment_samples <- sample_names[groups == contrast[2]]

    if (length(control_samples) < 1 || length(experiment_samples) < 1) {
        stop("No samples found for one of the contrast groups.")
    }

    as.vector(outer(control_samples, experiment_samples,
                    FUN = function(ctrl, exp) paste(ctrl, exp, sep = ",")))
}

.adam_aggregate_expression <- function(
    expr_mat,
    sample_data,
    aggregate_by,
    aggregate_fun
) {
    if (!all(aggregate_by %in% colnames(sample_data))) {
        stop("All aggregate_by fields must exist in sample metadata.")
    }

    grouping <- apply(sample_data[, aggregate_by, drop = FALSE], 1,
                    function(vals) paste(vals, collapse = "::"))

    if (aggregate_fun == "sum") {
        agg <- t(rowsum(t(expr_mat), group = grouping, reorder = FALSE))
    } else {
        agg_sum <- t(rowsum(t(expr_mat), group = grouping, reorder = FALSE))
        group_sizes <- table(grouping)[colnames(agg_sum)]
        agg <- sweep(agg_sum, 2, group_sizes, "/")
    }

    group_df <- unique(data.frame(grouping = grouping,
                                sample_data[, aggregate_by, drop = FALSE],
                                check.names = FALSE,
                                stringsAsFactors = FALSE))
    rownames(group_df) <- group_df$grouping
    group_df <- group_df[colnames(agg), aggregate_by, drop = FALSE]

    list(expr_mat = agg, sample_data = group_df)
}

.adam_to_summarized_experiment <- function(legacy_result, comparison_id, mode) {
    comparison_tables <- legacy_result[[2]]
    if (!is.list(comparison_tables)) {
        comparison_tables <- list(comparison_tables)
    }

    tables <- lapply(seq_along(comparison_tables), function(i) {
        tbl <- as.data.frame(comparison_tables[[i]], stringsAsFactors = FALSE)
        if (nrow(tbl) == 0L) {
            tbl$Comparison <- character(0)
        } else {
            tbl$Comparison <- rep(comparison_id[i], nrow(tbl))
        }
        tbl
    })

    all_cols <- unique(unlist(lapply(tables, colnames)))
    tables <- lapply(tables, function(tbl) {
        missing_cols <- setdiff(all_cols, colnames(tbl))
        if (length(missing_cols) > 0) {
            for (mc in missing_cols) {
                tbl[[mc]] <- NA
            }
        }
        tbl[, all_cols, drop = FALSE]
    })

    long_result <- do.call(rbind, tables)

    if (nrow(long_result) == 0) {
        return(SummarizedExperiment::SummarizedExperiment())
    }

    metric_cols <- vapply(long_result, is.numeric, logical(1))
    metric_mat <- as.matrix(long_result[, metric_cols, drop = FALSE])
    storage.mode(metric_mat) <- "double"

    row_ids <- make.unique(paste(long_result$ID, long_result$Comparison, sep = "|"))
    rownames(metric_mat) <- row_ids

    annotation <- long_result[, !metric_cols, drop = FALSE]
    row_data <- S4Vectors::DataFrame(annotation, row.names = row_ids,
                                    check.names = FALSE)

    result <- SummarizedExperiment::SummarizedExperiment(
        assays = list(metrics = metric_mat),
        rowData = row_data
    )

    S4Vectors::metadata(result)$comparison_id <- comparison_id
    S4Vectors::metadata(result)$mode <- mode
    S4Vectors::metadata(result)$gene_function_map <- legacy_result[[1]]
    S4Vectors::metadata(result)$legacy_result <- legacy_result

    result
}
