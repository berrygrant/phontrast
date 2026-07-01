#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || identical(x, "")) y else x
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript analysis/validation/summarize_validation_metrics.R \\",
      "    --metrics analysis/validation/outputs/run/validation_metrics.csv",
      "",
      "Or:",
      "  Rscript analysis/validation/summarize_validation_metrics.R \\",
      "    --validation-dir analysis/validation/outputs/run",
      "",
      "Optional:",
      "  --out-dir path/to/summary_dir",
      "",
      "Writes:",
      "  validation_summary_by_feature_set.csv",
      "  validation_summary_by_category_pair.csv",
      "  validation_metric_mode_counts.csv",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    key <- gsub("-", "_", key)
    if (identical(key, "help")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

read_nonblank_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Metrics file does not exist: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  nonblank <- lines[nzchar(trimws(lines))]
  if (!length(nonblank)) {
    stop(
      "Metrics file has no nonblank lines: ", path, "\n",
      "Check that you are pointing at the completed validation output directory.",
      call. = FALSE
    )
  }
  out <- tryCatch(
    read.csv(text = paste(nonblank, collapse = "\n"), stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    stop("Could not read metrics CSV: ", conditionMessage(out), call. = FALSE)
  }
  if (!nrow(out)) {
    stop(
      "Metrics file contains a header but no metric rows: ", path, "\n",
      "Inspect validation_run_summary.csv and validation_skipped_contrasts.csv for the skip reason.",
      call. = FALSE
    )
  }
  out
}

as_numeric_column <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  suppressWarnings(as.numeric(x))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 7))
}

contrast_ids <- function(df) {
  control <- if ("control_group" %in% names(df)) as.character(df$control_group) else ""
  cat1 <- if ("cat1" %in% names(df)) as.character(df$cat1) else ""
  cat2 <- if ("cat2" %in% names(df)) as.character(df$cat2) else ""
  paste(control, cat1, cat2, sep = "\t")
}

summarize_groups <- function(df, group_cols, metric_cols) {
  missing <- setdiff(group_cols, names(df))
  if (length(missing)) {
    stop("Metrics file is missing required grouping columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  key <- do.call(paste, c(df[group_cols], sep = "\r"))
  pieces <- split(seq_len(nrow(df)), key, drop = TRUE)

  rows <- lapply(pieces, function(idx) {
    base <- df[idx[[1L]], group_cols, drop = FALSE]
    base$n_rows <- length(idx)
    base$n_contrasts <- length(unique(contrast_ids(df[idx, , drop = FALSE])))
    if ("control_group" %in% names(df)) {
      base$n_control_groups <- length(unique(as.character(df$control_group[idx])))
    }
    if ("metric_mode" %in% names(df)) {
      mode <- as.character(df$metric_mode[idx])
      base$n_all_metrics <- sum(mode == "all_metrics", na.rm = TRUE)
      base$n_fallback <- sum(mode != "all_metrics" & !is.na(mode) & nzchar(mode), na.rm = TRUE)
    }

    for (col in metric_cols) {
      values <- as_numeric_column(df[[col]][idx])
      base[[paste0("median_", col)]] <- safe_median(values)
      base[[paste0("mean_", col)]] <- safe_mean(values)
      base[[paste0("q25_", col)]] <- safe_quantile(values, 0.25)
      base[[paste0("q75_", col)]] <- safe_quantile(values, 0.75)
    }
    base
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[do.call(order, out[group_cols]), , drop = FALSE]
}

metric_mode_counts <- function(df) {
  if (!all(c("feature_set", "metric_mode") %in% names(df))) {
    return(data.frame())
  }
  out <- as.data.frame(table(df$feature_set, df$metric_mode), stringsAsFactors = FALSE)
  names(out) <- c("feature_set", "metric_mode", "n")
  out <- out[out$n > 0L, , drop = FALSE]
  rownames(out) <- NULL
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(0L))
  }

  metrics_path <- args$metrics
  if (is.null(metrics_path) && !is.null(args$validation_dir)) {
    metrics_path <- file.path(args$validation_dir, "validation_metrics.csv")
  }
  if (is.null(metrics_path)) {
    usage()
    stop("Missing --metrics or --validation-dir.", call. = FALSE)
  }

  out_dir <- args$out_dir %||% dirname(metrics_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  metrics <- read_nonblank_csv(metrics_path)
  metric_cols <- intersect(
    c(
      "jsd", "js_distance", "percent_overlap", "classifier_auc",
      "classifier_balanced_accuracy", "classifier_accuracy",
      "pillai", "bhatt_dist", "mahalanobis_dist"
    ),
    names(metrics)
  )
  if (!length(metric_cols)) {
    stop("Metrics file has no recognized numeric metric columns.", call. = FALSE)
  }

  by_feature <- summarize_groups(metrics, "feature_set", metric_cols)
  by_pair <- summarize_groups(metrics, c("cat1", "cat2", "feature_set"), metric_cols)
  mode_counts <- metric_mode_counts(metrics)

  by_feature_path <- file.path(out_dir, "validation_summary_by_feature_set.csv")
  by_pair_path <- file.path(out_dir, "validation_summary_by_category_pair.csv")
  mode_counts_path <- file.path(out_dir, "validation_metric_mode_counts.csv")

  write.csv(by_feature, by_feature_path, row.names = FALSE)
  write.csv(by_pair, by_pair_path, row.names = FALSE)
  write.csv(mode_counts, mode_counts_path, row.names = FALSE)

  console_cols <- intersect(
    c(
      "feature_set", "n_rows", "n_contrasts", "n_control_groups",
      "n_all_metrics", "n_fallback", "median_jsd", "median_percent_overlap",
      "median_classifier_auc", "median_classifier_balanced_accuracy"
    ),
    names(by_feature)
  )

  cat("Read metric rows: ", nrow(metrics), "\n", sep = "")
  cat("Feature-set summary:\n")
  print(by_feature[console_cols], row.names = FALSE)
  if (nrow(mode_counts)) {
    cat("Metric-mode counts:\n")
    print(mode_counts, row.names = FALSE)
  }
  cat("Wrote: ", by_feature_path, "\n", sep = "")
  cat("Wrote: ", by_pair_path, "\n", sep = "")
  cat("Wrote: ", mode_counts_path, "\n", sep = "")
  invisible(0L)
}

if (identical(environment(), globalenv())) {
  main()
}
