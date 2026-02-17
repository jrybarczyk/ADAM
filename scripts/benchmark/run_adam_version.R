#!/usr/bin/env Rscript

parse_args <- function() {
    raw <- commandArgs(trailingOnly = TRUE)
    out <- list(
        engine = "local",
        mode = "partial",
        outdir = "benchmark_results",
        comparison = "control1,experiment1",
        min_gene = "3",
        max_gene = "20",
        bootstrap = "200",
        pcut = "0.05"
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
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
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

main <- function() {
    args <- parse_args()
    engine <- tolower(args$engine)
    mode <- tolower(args$mode)
    outdir <- normalizePath(args$outdir, mustWork = FALSE)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
    data("ExpressionAedes", package = "ADAM")
    data("KeggPathwaysAedes", package = "ADAM")

    comparison <- args$comparison
    min_gene <- as.integer(args$min_gene)
    max_gene <- as.integer(args$max_gene)
    bootstrap <- as.integer(args$bootstrap)
    pcut <- as.numeric(args$pcut)

    timing <- system.time({
        if (mode == "partial") {
            result <- ADAM::ADAnalysis(
                ComparisonID = comparison,
                ExpressionData = ExpressionAedes,
                MinGene = min_gene,
                MaxGene = max_gene,
                DBSpecies = KeggPathwaysAedes,
                AnalysisDomain = "own",
                GeneIdentifier = "geneStableID"
            )
        } else if (mode == "complete") {
            result <- ADAM::GFAGAnalysis(
                ComparisonID = comparison,
                ExpressionData = ExpressionAedes,
                MinGene = min_gene,
                MaxGene = max_gene,
                SeedNumber = 1049,
                BootstrapNumber = bootstrap,
                PCorrection = pcut,
                DBSpecies = KeggPathwaysAedes,
                PCorrectionMethod = "fdr",
                WilcoxonTest = FALSE,
                FisherTest = FALSE,
                AnalysisDomain = "own",
                GeneIdentifier = "geneStableID"
            )
        } else {
            stop("Invalid --mode. Use 'partial' or 'complete'.")
        }
    })

    comparison_table <- as.data.frame(result[[2]][[1]], stringsAsFactors = FALSE)
    pkg_version <- as.character(utils::packageVersion("ADAM"))

    summary_df <- data.frame(
        engine = engine,
        mode = mode,
        adam_version = pkg_version,
        n_terms = nrow(comparison_table),
        unique_ids = length(unique(comparison_table$ID)),
        elapsed_sec = unname(timing[["elapsed"]]),
        stringsAsFactors = FALSE
    )

    if ("n" %in% colnames(comparison_table)) {
        summary_df$mean_n <- mean(as.numeric(comparison_table$n), na.rm = TRUE)
    }
    if ("h" %in% colnames(comparison_table)) {
        summary_df$mean_h <- mean(as.numeric(comparison_table$h), na.rm = TRUE)
    }
    if ("qValue_n" %in% colnames(comparison_table)) {
        summary_df$signif_qn <- sum(as.numeric(comparison_table$qValue_n) <= pcut, na.rm = TRUE)
    }

    rds_file <- file.path(outdir, sprintf("result_%s_%s.rds", engine, mode))
    csv_file <- file.path(outdir, sprintf("table_%s_%s.csv", engine, mode))
    sum_file <- file.path(outdir, sprintf("summary_%s_%s.csv", engine, mode))
    info_file <- file.path(outdir, sprintf("session_%s_%s.txt", engine, mode))

    saveRDS(result, file = rds_file)
    utils::write.csv(comparison_table, file = csv_file, row.names = FALSE)
    utils::write.csv(summary_df, file = sum_file, row.names = FALSE)
    writeLines(capture.output(sessionInfo()), con = info_file)

    message("Finished.")
    message(sprintf("Result table: %s", csv_file))
    message(sprintf("Summary: %s", sum_file))
}

main()
