#!/usr/bin/env Rscript

parse_args <- function() {
    raw <- commandArgs(trailingOnly = TRUE)
    out <- list(
        outdir = "benchmark_results",
        mode = "partial",
        top_n = "20"
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

top_ids <- function(df, top_n) {
    if ("qValue_n" %in% colnames(df)) {
        ord <- order(as.numeric(df$qValue_n), na.last = NA)
        return(as.character(df$ID[ord])[seq_len(min(top_n, length(ord)))])
    }
    if ("n" %in% colnames(df)) {
        ord <- order(as.numeric(df$n), decreasing = TRUE, na.last = NA)
        return(as.character(df$ID[ord])[seq_len(min(top_n, length(ord)))])
    }
    as.character(df$ID[seq_len(min(top_n, nrow(df)))])
}

main <- function() {
    args <- parse_args()
    outdir <- normalizePath(args$outdir, mustWork = TRUE)
    mode <- tolower(args$mode)
    top_n <- as.integer(args$top_n)

    local_file <- file.path(outdir, sprintf("table_local_%s.csv", mode))
    bioc_file <- file.path(outdir, sprintf("table_bioc_%s.csv", mode))

    if (!file.exists(local_file) || !file.exists(bioc_file)) {
        stop("Missing input tables. Run run_adam_version.R for local and bioc first.")
    }

    local <- utils::read.csv(local_file, stringsAsFactors = FALSE)
    bioc <- utils::read.csv(bioc_file, stringsAsFactors = FALSE)

    merged <- merge(local, bioc, by = "ID", suffixes = c("_local", "_bioc"))

    num_local <- names(local)[vapply(local, is.numeric, logical(1))]
    num_bioc <- names(bioc)[vapply(bioc, is.numeric, logical(1))]
    numeric_common <- intersect(num_local, num_bioc)

    metric_rows <- lapply(numeric_common, function(col) {
        x <- as.numeric(merged[[paste0(col, "_local")]])
        y <- as.numeric(merged[[paste0(col, "_bioc")]])
        ok <- is.finite(x) & is.finite(y)
        if (!any(ok)) {
            return(data.frame(
                metric = col,
                spearman = NA_real_,
                mean_abs_diff = NA_real_,
                max_abs_diff = NA_real_,
                stringsAsFactors = FALSE
            ))
        }
        data.frame(
            metric = col,
            spearman = suppressWarnings(cor(x[ok], y[ok], method = "spearman")),
            mean_abs_diff = mean(abs(x[ok] - y[ok]), na.rm = TRUE),
            max_abs_diff = max(abs(x[ok] - y[ok]), na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })
    metric_df <- do.call(rbind, metric_rows)

    top_local <- unique(top_ids(local, top_n))
    top_bioc <- unique(top_ids(bioc, top_n))
    inter <- length(intersect(top_local, top_bioc))
    uni <- length(unique(c(top_local, top_bioc)))
    jaccard <- if (uni == 0) NA_real_ else inter / uni

    overview <- data.frame(
        mode = mode,
        n_local = nrow(local),
        n_bioc = nrow(bioc),
        n_shared = nrow(merged),
        top_n = top_n,
        top_overlap = inter,
        top_jaccard = jaccard,
        stringsAsFactors = FALSE
    )

    utils::write.csv(metric_df, file = file.path(outdir, sprintf("metrics_%s.csv", mode)), row.names = FALSE)
    utils::write.csv(overview, file = file.path(outdir, sprintf("overview_%s.csv", mode)), row.names = FALSE)

    message("Comparison finished.")
    message(sprintf("Overview: %s", file.path(outdir, sprintf("overview_%s.csv", mode))))
    message(sprintf("Metrics: %s", file.path(outdir, sprintf("metrics_%s.csv", mode))))
}

main()
