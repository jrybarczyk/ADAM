ensure_biocmanager <- function() {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
}

ensure_packages <- function(pkgs, bioc = TRUE) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) == 0L) {
        return(invisible(TRUE))
    }

    if (bioc) {
        ensure_biocmanager()
        BiocManager::install(missing, ask = FALSE, update = FALSE)
    } else {
        install.packages(missing, repos = "https://cloud.r-project.org")
    }

    invisible(TRUE)
}

pick_top_two_groups <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) < 2L) {
        stop("Not enough valid group labels to define a contrast.")
    }
    freq <- sort(table(x), decreasing = TRUE)
    if (length(freq) < 2L) {
        stop("At least two groups are required to define a contrast.")
    }
    names(freq)[1:2]
}

make_pairwise_comparisons <- function(control_samples, experiment_samples, max_pairs = 4L) {
    pairs <- as.vector(outer(control_samples, experiment_samples,
                            FUN = function(ctrl, exp) paste(ctrl, exp, sep = ",")))
    pairs[seq_len(min(length(pairs), max_pairs))]
}

build_data_driven_db <- function(expr_mat, n_groups = 40L, seed = 1049L, prefix = "MOD") {
    expr_mat <- as.matrix(expr_mat)
    if (is.null(rownames(expr_mat))) {
        stop("Expression matrix must include row names (gene identifiers).")
    }

    keep <- apply(expr_mat, 1, function(row_vals) {
        all(is.finite(row_vals)) && stats::var(row_vals) > 0
    })
    expr_mat <- expr_mat[keep, , drop = FALSE]

    if (nrow(expr_mat) < 10L) {
        stop("Not enough informative genes after filtering to build modules.")
    }

    scaled <- t(scale(t(log1p(abs(expr_mat)))))
    scaled[is.na(scaled)] <- 0

    centers <- max(2L, min(as.integer(n_groups), floor(nrow(scaled) / 3L)))
    set.seed(seed)
    km <- stats::kmeans(scaled, centers = centers, nstart = 5, iter.max = 100)

    module_id <- sprintf("%s_%03d", prefix, km$cluster)
    data.frame(
        gene = rownames(scaled),
        ID = module_id,
        Description = paste("Data-driven module", module_id),
        stringsAsFactors = FALSE
    )
}
