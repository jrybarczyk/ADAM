#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1L) args[[1]] else "benchmark_results"
base <- normalizePath(base, mustWork = FALSE)

read_all <- function(files) {
    if (length(files) == 0L) {
        return(data.frame())
    }
    lst <- lapply(files, function(f) {
        d <- utils::read.csv(f, stringsAsFactors = FALSE)
        d$source_file <- basename(f)
        d
    })
    all_cols <- unique(unlist(lapply(lst, colnames)))
    lst <- lapply(lst, function(d) {
        miss <- setdiff(all_cols, colnames(d))
        if (length(miss) > 0L) {
            for (m in miss) d[[m]] <- NA
        }
        d[, all_cols, drop = FALSE]
    })
    do.call(rbind, lst)
}

overview_files <- Sys.glob(file.path(base, "overview_*.csv"))
metrics_files <- Sys.glob(file.path(base, "metrics_*.csv"))
summary_files <- Sys.glob(file.path(base, "summary_*.csv"))

over <- read_all(overview_files)
met <- read_all(metrics_files)
sumdf <- read_all(summary_files)

if (nrow(sumdf) > 0L) {
    sumdf$scenario <- paste(sumdf$dataset, sumdf$mode, sep = "|")
}
if (nrow(over) > 0L) {
    over$scenario <- paste(over$dataset, over$mode, sep = "|")
}
if (nrow(met) > 0L) {
    met$scenario <- sub("^metrics_(.*)\\.csv$", "\\1", met$source_file)
}

utils::write.csv(over, file.path(base, "consolidated_overview.csv"), row.names = FALSE)
utils::write.csv(met, file.path(base, "consolidated_metrics.csv"), row.names = FALSE)
utils::write.csv(sumdf, file.path(base, "consolidated_summary.csv"), row.names = FALSE)

exec <- data.frame()
if (nrow(sumdf) > 0L) {
    keys <- unique(sumdf$scenario)
    rows <- lapply(keys, function(k) {
        x <- sumdf[sumdf$scenario == k, , drop = FALSE]
        l <- x[x$engine == "local", , drop = FALSE]
        b <- x[x$engine == "bioc", , drop = FALSE]
        if (nrow(l) == 0L || nrow(b) == 0L) {
            return(NULL)
        }
        data.frame(
            scenario = k,
            local_version = l$adam_version[1],
            bioc_version = b$adam_version[1],
            local_elapsed_sec = l$elapsed_sec[1],
            bioc_elapsed_sec = b$elapsed_sec[1],
            elapsed_ratio_bioc_vs_local = if (isTRUE(all.equal(l$elapsed_sec[1], 0))) NA_real_ else b$elapsed_sec[1] / l$elapsed_sec[1],
            local_n_terms = l$n_terms[1],
            bioc_n_terms = b$n_terms[1],
            stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) > 0L) {
        exec <- do.call(rbind, rows)
    }
}
utils::write.csv(exec, file.path(base, "consolidated_exec.csv"), row.names = FALSE)

md <- file.path(base, "consolidated_report.md")
lines <- c(
    "# Relatorio Consolidado de Benchmark ADAM",
    "",
    sprintf("- Data de geracao: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("- Cenarios de overview: %d", nrow(over)),
    sprintf("- Linhas de metricas: %d", nrow(met)),
    sprintf("- Linhas de resumo (engine): %d", nrow(sumdf)),
    ""
)

if (nrow(exec) > 0L) {
    lines <- c(
        lines,
        "## Resumo Executivo",
        "",
        "| Cenario | Versao Local | Versao Bioc | Tempo Local (s) | Tempo Bioc (s) | Razao Bioc/Local | Termos Local | Termos Bioc |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"
    )
    for (i in seq_len(nrow(exec))) {
        r <- exec[i, , drop = FALSE]
        lines <- c(lines, sprintf(
            "| %s | %s | %s | %.4f | %.4f | %.4f | %d | %d |",
            r$scenario, r$local_version, r$bioc_version,
            as.numeric(r$local_elapsed_sec), as.numeric(r$bioc_elapsed_sec),
            as.numeric(r$elapsed_ratio_bioc_vs_local),
            as.integer(r$local_n_terms), as.integer(r$bioc_n_terms)
        ))
    }
    lines <- c(lines, "")
}

if (nrow(over) > 0L) {
    lines <- c(
        lines,
        "## Overlap de Resultados (Overview)",
        "",
        "| Dataset | Modo | n_local | n_bioc | n_shared | top_n | top_overlap | top_jaccard |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"
    )
    for (i in seq_len(nrow(over))) {
        r <- over[i, , drop = FALSE]
        lines <- c(lines, sprintf(
            "| %s | %s | %d | %d | %d | %d | %d | %.4f |",
            r$dataset, r$mode, as.integer(r$n_local), as.integer(r$n_bioc),
            as.integer(r$n_shared), as.integer(r$top_n), as.integer(r$top_overlap),
            as.numeric(r$top_jaccard)
        ))
    }
    lines <- c(lines, "")
}

if (nrow(met) > 0L) {
    lines <- c(
        lines,
        "## Metricas Numericas (medias por cenario)",
        "",
        "| Cenario | n_metricas | spearman_medio | mean_abs_diff_medio | max_abs_diff_max |",
        "|---|---:|---:|---:|---:|"
    )
    scenarios <- unique(met$scenario)
    for (sc in scenarios) {
        msub <- met[met$scenario == sc, , drop = FALSE]
        sp_mean <- mean(as.numeric(msub$spearman), na.rm = TRUE)
        mad_mean <- mean(as.numeric(msub$mean_abs_diff), na.rm = TRUE)
        mxd_max <- max(as.numeric(msub$max_abs_diff), na.rm = TRUE)
        lines <- c(lines, sprintf(
            "| %s | %d | %.4f | %.6f | %.6f |",
            sc, nrow(msub), sp_mean, mad_mean, mxd_max
        ))
    }
    lines <- c(lines, "")
}

lines <- c(
    lines,
    "## Arquivos Gerados",
    "",
    "- `benchmark_results/consolidated_report.md`",
    "- `benchmark_results/consolidated_exec.csv`",
    "- `benchmark_results/consolidated_overview.csv`",
    "- `benchmark_results/consolidated_metrics.csv`",
    "- `benchmark_results/consolidated_summary.csv`",
    ""
)

writeLines(lines, md)

cat("Consolidated report created:\n")
cat(file.path(base, "consolidated_report.md"), "\n")
cat(file.path(base, "consolidated_exec.csv"), "\n")
