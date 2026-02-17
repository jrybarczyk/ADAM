#!/usr/bin/env Rscript

parse_args <- function() {
    raw <- commandArgs(trailingOnly = TRUE)
    out <- list(
        engine = "local",
        dataset = "airway",
        mode = "partial",
        outdir = "benchmark_results",
        max_genes = "3000",
        min_gene = "3",
        max_gene = "20",
        bootstrap = "200",
        pcut = "0.05",
        pairing = "first"
    )

    i <- 1L
    while (i <= length(raw)) {
        key <- raw[[i]]
        if (!startsWith(key, "--")) {
            i <- i + 1L
            next
        }
        key <- sub("^--", "", key)
        if ((i + 1L) > length(raw)) {
            stop(sprintf("Missing value for --%s", key))
        }
        out[[gsub("-", "_", key)]] <- raw[[i + 1L]]
        i <- i + 2L
    }
    out
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_package <- function(pkg, bioc = FALSE) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        return(invisible(TRUE))
    }
    if (!bioc) {
        stop(sprintf("Package '%s' is required but not installed.", pkg))
    }
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        stop("Package 'BiocManager' is required to install missing Bioconductor packages.")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Failed to install required package '%s'.", pkg))
    }
    invisible(TRUE)
}

install_local_adam <- function(lib, pkg_dir) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    cmd <- file.path(R.home("bin"), "R")
    args <- c("CMD", "INSTALL", sprintf("--library=%s", lib), pkg_dir)
    status <- system2(cmd, args = args)
    if (!identical(status, 0L)) {
        stop("Failed to install local ADAM source.")
    }
}

install_bioc_adam <- function(lib) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    ensure_package("BiocManager")
    ok <- TRUE
    tryCatch(
        BiocManager::install("ADAM", ask = FALSE, update = FALSE, lib = lib),
        error = function(e) {
            ok <<- FALSE
            message(sprintf("Bioc install failed: %s", conditionMessage(e)))
        }
    )
    if (!ok) {
        if (requireNamespace("ADAM", quietly = TRUE)) {
            message("Using preinstalled ADAM from current R library as 'bioc' engine.")
            return(invisible(FALSE))
        }
        stop("Could not install or find preinstalled ADAM for bioc engine.")
    }
    TRUE
}

limit_genes <- function(mat, max_genes) {
    if (nrow(mat) <= max_genes) {
        return(mat)
    }
    vars <- apply(mat, 1, stats::var, na.rm = TRUE)
    keep <- order(vars, decreasing = TRUE)[seq_len(max_genes)]
    mat[keep, , drop = FALSE]
}

build_own_db <- function(genes, term_size = 30L, step = 10L, max_terms = 400L) {
    genes <- unique(as.character(genes))
    if (length(genes) < term_size) {
        term_size <- max(3L, floor(length(genes) / 2L))
    }
    starts <- seq.int(1L, max(1L, length(genes) - term_size + 1L), by = step)
    if (length(starts) > max_terms) {
        starts <- starts[seq_len(max_terms)]
    }

    rows <- lapply(seq_along(starts), function(i) {
        s <- starts[[i]]
        e <- min(length(genes), s + term_size - 1L)
        ids <- genes[s:e]
        term <- sprintf("TERM%04d", i)
        data.frame(
            gene = ids,
            ID = term,
            Description = sprintf("Synthetic term %s", term),
            stringsAsFactors = FALSE
        )
    })

    unique(do.call(rbind, rows))
}

build_comparison_id <- function(sample_names, groups, pairing = c("first", "all")) {
    pairing <- match.arg(pairing)
    groups <- as.character(groups)
    levs <- unique(groups)
    if (length(levs) < 2L) {
        stop("Need at least 2 groups to define comparisons.")
    }
    control <- sample_names[groups == levs[1]]
    experiment <- sample_names[groups == levs[2]]
    if (length(control) == 0L || length(experiment) == 0L) {
        stop("No samples found in one of the two groups.")
    }
    if (pairing == "first") {
        return(sprintf("%s,%s", control[1], experiment[1]))
    }
    as.vector(outer(control, experiment, FUN = function(cn, en) paste(cn, en, sep = ",")))
}

prepare_airway <- function(max_genes, pairing) {
    ensure_package("airway", bioc = TRUE)
    ensure_package("SummarizedExperiment", bioc = TRUE)
    data("airway", package = "airway")
    se <- airway::airway
    mat <- SummarizedExperiment::assay(se)
    groups <- as.character(SummarizedExperiment::colData(se)$dex)
    mat <- limit_genes(mat, max_genes)
    comp <- build_comparison_id(colnames(mat), groups, pairing = pairing)
    expr <- data.frame(gene = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE)
    list(ExpressionData = expr, DBSpecies = build_own_db(rownames(mat)), ComparisonID = comp)
}

prepare_all <- function(max_genes, pairing) {
    ensure_package("ALL", bioc = TRUE)
    ensure_package("Biobase", bioc = TRUE)
    data("ALL", package = "ALL")
    eset <- ALL::ALL
    mat <- Biobase::exprs(eset)
    groups <- as.character(Biobase::pData(eset)$BT)
    keep <- groups %in% c("B", "T")
    mat <- mat[, keep, drop = FALSE]
    groups <- groups[keep]
    # Make values non-negative to keep activity/diversity metrics stable.
    min_val <- min(mat, na.rm = TRUE)
    if (is.finite(min_val) && min_val < 0) {
        mat <- mat - min_val + 1e-06
    }
    mat <- limit_genes(mat, max_genes)
    comp <- build_comparison_id(colnames(mat), groups, pairing = pairing)
    expr <- data.frame(gene = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE)
    list(ExpressionData = expr, DBSpecies = build_own_db(rownames(mat)), ComparisonID = comp)
}

prepare_sce_mock <- function(max_genes, pairing) {
    ensure_package("scuttle", bioc = TRUE)
    ensure_package("SingleCellExperiment", bioc = TRUE)
    sce <- scuttle::mockSCE(ncells = 400, ngenes = max_genes)
    mat <- SummarizedExperiment::assay(sce, "counts")
    mat <- limit_genes(mat, max_genes)
    nc <- ncol(mat)
    cond <- ifelse(seq_len(nc) %% 2L == 0L, "experiment", "control")
    pseudo <- paste0(cond, "_", ceiling(seq_len(nc) / 25L))
    agg <- t(rowsum(t(mat), group = pseudo, reorder = FALSE))
    groups <- ifelse(grepl("^control_", colnames(agg)), "control", "experiment")
    comp <- build_comparison_id(colnames(agg), groups, pairing = pairing)
    expr <- data.frame(gene = rownames(agg), as.data.frame(agg, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE)
    list(ExpressionData = expr, DBSpecies = build_own_db(rownames(agg)), ComparisonID = comp)
}

prepare_aedes <- function(max_genes, pairing) {
    data("ExpressionAedes", package = "ADAM")
    data("KeggPathwaysAedes", package = "ADAM")
    mat <- as.matrix(ExpressionAedes[, -1, drop = FALSE])
    rownames(mat) <- as.character(ExpressionAedes$gene)
    mat <- limit_genes(mat, max_genes)
    groups <- ifelse(grepl("^control", colnames(mat), ignore.case = TRUE), "control", "experiment")
    comp <- build_comparison_id(colnames(mat), groups, pairing = pairing)
    expr <- data.frame(gene = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE)
    list(ExpressionData = expr, DBSpecies = KeggPathwaysAedes, ComparisonID = comp, gene_id = "geneStableID")
}

prepare_dataset <- function(dataset, max_genes, pairing) {
    dataset <- tolower(dataset)
    if (dataset == "airway") {
        return(prepare_airway(max_genes = max_genes, pairing = pairing))
    }
    if (dataset == "all") {
        return(prepare_all(max_genes = max_genes, pairing = pairing))
    }
    if (dataset == "sce_mock") {
        return(prepare_sce_mock(max_genes = max_genes, pairing = pairing))
    }
    if (dataset == "aedes") {
        return(prepare_aedes(max_genes = max_genes, pairing = pairing))
    }
    stop("Invalid --dataset. Use 'airway', 'all', 'sce_mock', or 'aedes'.")
}

run_adam <- function(mode, ds, min_gene, max_gene, bootstrap, pcut) {
    if (mode == "partial") {
        return(ADAM::ADAnalysis(
            ComparisonID = ds$ComparisonID,
            ExpressionData = ds$ExpressionData,
            MinGene = min_gene,
            MaxGene = max_gene,
            DBSpecies = ds$DBSpecies,
            AnalysisDomain = "own",
            GeneIdentifier = ds$gene_id %||% "gene"
        ))
    }
    if (mode == "complete") {
        return(ADAM::GFAGAnalysis(
            ComparisonID = ds$ComparisonID,
            ExpressionData = ds$ExpressionData,
            MinGene = min_gene,
            MaxGene = max_gene,
            SeedNumber = 1049,
            BootstrapNumber = bootstrap,
            PCorrection = pcut,
            DBSpecies = ds$DBSpecies,
            PCorrectionMethod = "fdr",
            WilcoxonTest = FALSE,
            FisherTest = FALSE,
            AnalysisDomain = "own",
            GeneIdentifier = ds$gene_id %||% "gene"
        ))
    }
    stop("Invalid --mode. Use 'partial' or 'complete'.")
}

main <- function() {
    args <- parse_args()
    engine <- tolower(args$engine)
    dataset <- tolower(args$dataset)
    mode <- tolower(args$mode)
    outdir <- normalizePath(args$outdir, mustWork = FALSE)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    max_genes <- as.integer(args$max_genes)
    min_gene <- as.integer(args$min_gene)
    max_gene <- as.integer(args$max_gene)
    bootstrap <- as.integer(args$bootstrap)
    pcut <- as.numeric(args$pcut)
    pairing <- tolower(args$pairing)

    pkg_dir <- normalizePath(".", mustWork = TRUE)
    lib <- file.path(outdir, sprintf("lib_%s", engine))

    use_engine_lib <- TRUE
    if (engine == "local") {
        message("Installing local ADAM...")
        install_local_adam(lib = lib, pkg_dir = pkg_dir)
    } else if (engine == "bioc") {
        message("Installing Bioconductor ADAM...")
        use_engine_lib <- isTRUE(install_bioc_adam(lib = lib))
    } else {
        stop("Invalid --engine. Use 'local' or 'bioc'.")
    }

    if (use_engine_lib) {
        library(ADAM, lib.loc = lib)
    } else {
        library(ADAM)
    }
    ds <- prepare_dataset(dataset = dataset, max_genes = max_genes, pairing = pairing)

    timing <- system.time({
        result <- run_adam(
            mode = mode,
            ds = ds,
            min_gene = min_gene,
            max_gene = max_gene,
            bootstrap = bootstrap,
            pcut = pcut
        )
    })

    result_table <- as.data.frame(result[[2]][[1]], stringsAsFactors = FALSE)
    pkg_version <- as.character(utils::packageVersion("ADAM"))

    summary_df <- data.frame(
        dataset = dataset,
        engine = engine,
        mode = mode,
        adam_version = pkg_version,
        n_terms = nrow(result_table),
        unique_ids = length(unique(result_table$ID)),
        n_comparisons = length(ds$ComparisonID),
        elapsed_sec = unname(timing[["elapsed"]]),
        stringsAsFactors = FALSE
    )

    if ("n" %in% colnames(result_table)) {
        summary_df$mean_n <- mean(as.numeric(result_table$n), na.rm = TRUE)
    }
    if ("h" %in% colnames(result_table)) {
        summary_df$mean_h <- mean(as.numeric(result_table$h), na.rm = TRUE)
    }
    if ("qValue_n" %in% colnames(result_table)) {
        summary_df$signif_qn <- sum(as.numeric(result_table$qValue_n) <= pcut, na.rm = TRUE)
    }

    prefix <- sprintf("%s_%s_%s", dataset, engine, mode)
    saveRDS(result, file = file.path(outdir, sprintf("result_%s.rds", prefix)))
    utils::write.csv(result_table, file = file.path(outdir, sprintf("table_%s.csv", prefix)), row.names = FALSE)
    utils::write.csv(summary_df, file = file.path(outdir, sprintf("summary_%s.csv", prefix)), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), con = file.path(outdir, sprintf("session_%s.txt", prefix)))

    message("Finished.")
    message(sprintf("Table: %s", file.path(outdir, sprintf("table_%s.csv", prefix))))
    message(sprintf("Summary: %s", file.path(outdir, sprintf("summary_%s.csv", prefix))))
}

main()
