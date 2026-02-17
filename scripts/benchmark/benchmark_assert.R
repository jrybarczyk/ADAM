#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
    idx <- which(args == flag)
    if (length(idx) == 0L) return(default)
    if (idx[length(idx)] == length(args)) stop(sprintf("Missing value for %s", flag))
    args[idx[length(idx)] + 1L]
}

outdir <- normalizePath(get_arg("--outdir", "benchmark_results"), mustWork = FALSE)
min_jaccard <- as.numeric(get_arg("--min-jaccard", "0.90"))
min_shared_ratio <- as.numeric(get_arg("--min-shared-ratio", "0.90"))
min_spearman <- as.numeric(get_arg("--min-spearman", "0.95"))
max_mean_abs_diff <- as.numeric(get_arg("--max-mean-abs-diff", "0.05"))

overview_files <- Sys.glob(file.path(outdir, "overview_*.csv"))
metrics_files <- Sys.glob(file.path(outdir, "metrics_*.csv"))

if (length(overview_files) == 0L || length(metrics_files) == 0L) {
    stop("No overview/metrics files found. Run benchmark scripts first.")
}

read_all <- function(files) {
    do.call(rbind, lapply(files, function(f) {
        d <- utils::read.csv(f, stringsAsFactors = FALSE)
        d$source_file <- basename(f)
        d
    }))
}

overview <- read_all(overview_files)
metrics <- read_all(metrics_files)
metrics$scenario <- sub("^metrics_(.*)\\.csv$", "\\1", metrics$source_file)
overview$scenario <- sub("^overview_(.*)\\.csv$", "\\1", overview$source_file)

rows <- lapply(unique(overview$scenario), function(sc) {
    o <- overview[overview$scenario == sc, , drop = FALSE]
    m <- metrics[metrics$scenario == sc, , drop = FALSE]

    shared_ratio <- if (nrow(o) > 0L) {
        as.numeric(o$n_shared[1]) / max(as.numeric(o$n_local[1]), as.numeric(o$n_bioc[1]))
    } else {
        NA_real_
    }
    jaccard <- if (nrow(o) > 0L) as.numeric(o$top_jaccard[1]) else NA_real_

    sp <- suppressWarnings(as.numeric(m$spearman))
    sp <- sp[is.finite(sp)]
    spearman_mean <- if (length(sp) > 0L) mean(sp) else NA_real_

    mad <- suppressWarnings(as.numeric(m$mean_abs_diff))
    mad <- mad[is.finite(mad)]
    mean_abs_diff_mean <- if (length(mad) > 0L) mean(mad) else NA_real_

    checks <- c(
        isTRUE(shared_ratio >= min_shared_ratio),
        isTRUE(jaccard >= min_jaccard),
        (is.na(spearman_mean) || spearman_mean >= min_spearman),
        (is.na(mean_abs_diff_mean) || mean_abs_diff_mean <= max_mean_abs_diff)
    )

    status <- if (all(checks)) "PASS" else "FAIL"
    reason <- c()
    if (!checks[1]) reason <- c(reason, sprintf("shared_ratio<%.2f", min_shared_ratio))
    if (!checks[2]) reason <- c(reason, sprintf("top_jaccard<%.2f", min_jaccard))
    if (!checks[3]) reason <- c(reason, sprintf("spearman_mean<%.2f", min_spearman))
    if (!checks[4]) reason <- c(reason, sprintf("mean_abs_diff_mean>%.4f", max_mean_abs_diff))

    data.frame(
        scenario = sc,
        status = status,
        shared_ratio = shared_ratio,
        top_jaccard = jaccard,
        spearman_mean = spearman_mean,
        mean_abs_diff_mean = mean_abs_diff_mean,
        reason = paste(reason, collapse = ";"),
        stringsAsFactors = FALSE
    )
})

result <- do.call(rbind, rows)
csv_out <- file.path(outdir, "benchmark_assertion.csv")
md_out <- file.path(outdir, "benchmark_assertion.md")
utils::write.csv(result, csv_out, row.names = FALSE)

lines <- c(
    "# Benchmark Assertion",
    "",
    sprintf("- min_shared_ratio: %.2f", min_shared_ratio),
    sprintf("- min_jaccard: %.2f", min_jaccard),
    sprintf("- min_spearman: %.2f", min_spearman),
    sprintf("- max_mean_abs_diff: %.4f", max_mean_abs_diff),
    "",
    "| Scenario | Status | Shared Ratio | Top Jaccard | Spearman Mean | Mean Abs Diff Mean | Reason |",
    "|---|---|---:|---:|---:|---:|---|"
)
for (i in seq_len(nrow(result))) {
    r <- result[i, , drop = FALSE]
    lines <- c(lines, sprintf(
        "| %s | %s | %.4f | %.4f | %s | %s | %s |",
        r$scenario, r$status, as.numeric(r$shared_ratio), as.numeric(r$top_jaccard),
        ifelse(is.na(r$spearman_mean), "NA", sprintf("%.4f", as.numeric(r$spearman_mean))),
        ifelse(is.na(r$mean_abs_diff_mean), "NA", sprintf("%.6f", as.numeric(r$mean_abs_diff_mean))),
        ifelse(nchar(r$reason) == 0, "-", r$reason)
    ))
}
writeLines(lines, md_out)

cat(sprintf("Assertion report: %s\n", csv_out))
cat(sprintf("Assertion markdown: %s\n", md_out))

if (any(result$status == "FAIL")) {
    quit(save = "no", status = 1)
}
