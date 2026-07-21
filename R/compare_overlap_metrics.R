#' Compute and compare phonological contrast metrics
#'
#' \code{phontrast()} is the package's main entry point. It computes one or more
#' category separation and overlap metrics for a two-category phonological
#' contrast in a single call: Jensen-Shannon divergence and distance, the
#' Pillai-Bartlett trace, Bhattacharyya distance and affinity, Mahalanobis
#' distance, and proportional overlap. Choose the metrics you want with
#' \code{metrics}; the default computes all of them. Results are returned
#' globally or by group, in a wide format (one column per metric, the default)
#' or a tidy long format (one row per metric per comparison). The
#' \code{percent_overlap} values are 0--1 proportions, not 0--100 percentages.
#'
#' Use \code{estimate_jsd()} when Jensen-Shannon divergence is the only outcome
#' of interest, and the lower-level metric helpers when you need direct control
#' over one estimator.
#'
#' Metric directions differ. JSD, Jensen-Shannon distance, Pillai trace,
#' Bhattacharyya distance, and Mahalanobis distance increase as categories
#' become more separated. Percent overlap and Bhattacharyya affinity increase
#' as categories overlap more. Long output includes \code{orientation},
#' \code{separation_value}, and \code{separation_rank} columns so all metrics
#' can be read on a separation-oriented scale.
#'
#' If \code{do_boot = TRUE}, each metric is recomputed on \code{n_boot}
#' nonparametric bootstrap resamples to estimate uncertainty. This can take
#' substantial time because every resample recomputes KDE, MANOVA, and
#' covariance-based metrics. Progress messages are printed by default while
#' bootstrapping is running; set \code{progress = FALSE} to suppress them.
#'
#' @param data Data frame containing category labels and acoustic features.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; column giving the two categories to compare.
#' @param group_col Optional character vector of one or more grouping columns.
#'   If \code{NULL}, metrics are computed globally. Multiple grouping columns
#'   are combined into a labeled \code{group} value such as
#'   \code{"Sex=F | Style=read"}.
#' @param metrics Character vector selecting which contrast metrics to compute.
#'   Any of \code{"jsd"}, \code{"js_distance"}, \code{"pillai"},
#'   \code{"bhattacharyya"}, \code{"mahalanobis"}, and \code{"overlap"}.
#'   Defaults to all of them. \code{"bhattacharyya"} returns both the
#'   Bhattacharyya distance and affinity.
#' @param min_tokens Minimum tokens required globally or per group.
#' @param bw Bandwidth selection method passed to \code{jsd_kde_nd()} and
#'   \code{percent_overlap_kde()}.
#' @param eval_on KDE evaluation points passed to \code{jsd_kde_nd()} and
#'   \code{percent_overlap_kde()}.
#' @param eval_n Optional maximum number of KDE evaluation points passed to
#'   \code{jsd_kde_nd()} and \code{percent_overlap_kde()}.
#' @param eval_seed Optional integer seed for KDE evaluation-point subsampling.
#' @param engine KDE evaluation engine passed to \code{jsd_kde_nd()} and
#'   \code{percent_overlap_kde()}. \code{"fast_diagonal"} is accepted as an
#'   alias for \code{"fast_diag"}.
#' @param chunk_size Chunk size for \code{engine = "fast_diag"}.
#' @param eps Small ridge constant for covariance-based metrics.
#' @param output Output format: \code{"wide"} returns one row per global/group
#'   comparison; \code{"long"} returns one row per metric per comparison.
#' @param do_boot Logical; if \code{TRUE}, compute bootstrap means, standard
#'   deviations, and confidence intervals for each reported metric.
#' @param n_boot Number of bootstrap resamples if \code{do_boot = TRUE}.
#' @param conf_level Confidence level for bootstrap intervals.
#' @param progress Logical; if \code{TRUE}, print progress messages while
#'   bootstrap resamples are running.
#' @param method KDE estimator for the JSD and percent-overlap columns, passed
#'   to \code{jsd_kde_nd()}/\code{percent_overlap_kde()}: \code{"mc"} (default)
#'   for the Monte-Carlo plug-in, or \code{"legacy"} for the pre-1.2.0
#'   self-normalized estimate.
#'
#' @return A data frame containing only the requested \code{metrics}. Wide
#'   output (the default) contains one column per requested metric plus
#'   \code{pillai_p_value} when Pillai is requested; with \code{do_boot = TRUE}
#'   it also includes metric-specific \code{*_mean}, \code{*_sd},
#'   \code{*_ci_lower}, \code{*_ci_upper}, and \code{*_n_boot} columns. Long
#'   output contains \code{metric}, \code{estimate}, \code{orientation},
#'   \code{bounded_0_1}, \code{separation_value}, \code{separation_rank}, and
#'   \code{p_value} (populated for the Pillai row, \code{NA} otherwise) columns;
#'   with \code{do_boot = TRUE} it also includes \code{boot_mean},
#'   \code{boot_sd}, \code{ci_lower}, \code{ci_upper}, \code{n_boot}, and
#'   \code{conf_level}.
#'
#' @examples
#' set.seed(2026)
#' vowels <- data.frame(
#'   speaker = rep(c("s01", "s02"), each = 60),
#'   vowel = rep(rep(c("ih", "eh"), each = 30), 2),
#'   f1 = c(
#'     rnorm(30, 500, 55), rnorm(30, 560, 60),
#'     rnorm(30, 510, 60), rnorm(30, 575, 65)
#'   ),
#'   f2 = c(
#'     rnorm(30, 1980, 150), rnorm(30, 1880, 155),
#'     rnorm(30, 1960, 160), rnorm(30, 1840, 165)
#'   )
#' )
#'
#' # All metrics in one wide comparison table (the default), by speaker.
#' phontrast(
#'   data = vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   group_col = "speaker"
#' )
#'
#' # A single metric in wide format.
#' phontrast(
#'   data = vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   group_col = "speaker",
#'   metrics = "pillai",
#'   output = "wide"
#' )
#'
#' \donttest{
#' # Bootstrapping is useful but slower because every requested metric is
#' # recomputed on every resample. Use a larger n_boot for real analyses.
#' phontrast(
#'   data = vowels,
#'   features = "f1",
#'   category_col = "vowel",
#'   group_col = "speaker",
#'   metrics = c("jsd", "pillai"),
#'   do_boot = TRUE,
#'   n_boot = 5,
#'   progress = FALSE
#' )
#' }
#' @export
phontrast <- function(data,
                      features,
                      category_col,
                      group_col = NULL,
                      metrics = c("jsd", "js_distance", "pillai",
                                  "bhattacharyya", "mahalanobis", "overlap"),
                      min_tokens = 20,
                      bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                      eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                      eval_n = NULL,
                      eval_seed = NULL,
                      engine = c("ks", "fast_diag", "fast_diagonal"),
                      chunk_size = 1000L,
                      eps = 1e-6,
                      output = c("wide", "long"),
                      do_boot = FALSE,
                      n_boot = 1000,
                      conf_level = 0.95,
                      progress = TRUE,
                      method = c("mc", "legacy")) {
  output <- match.arg(output)
  metrics <- .resolve_contrast_metrics(metrics)
  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  engine <- .match_kde_engine(engine)
  method <- match.arg(method)
  .check_positive_count(min_tokens, "min_tokens")
  .check_ridge_eps(eps, "eps")
  if (!is.logical(do_boot) || length(do_boot) != 1L || is.na(do_boot)) {
    stop("`do_boot` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  .check_conf_level(conf_level)
  if (isTRUE(do_boot)) {
    .check_positive_count(n_boot, "n_boot")
  }
  .validate_metric_inputs(data, features, category_col, group_col)

  wide <- .compare_overlap_metrics_point(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    min_tokens = min_tokens,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size,
    eps = eps,
    method = method
  )
  if (!nrow(wide)) {
    .warn_empty_overlap_comparison(
      data = data,
      features = features,
      category_col = category_col,
      group_col = group_col,
      min_tokens = min_tokens
    )
  }

  if (isTRUE(do_boot) && nrow(wide)) {
    boot <- .bootstrap_compare_overlap_metrics(
      data = data,
      point_wide = wide,
      features = features,
      category_col = category_col,
      group_col = group_col,
      min_tokens = min_tokens,
      bw = bw,
      eval_on = eval_on,
      eval_n = eval_n,
      eval_seed = eval_seed,
      engine = engine,
      chunk_size = chunk_size,
      eps = eps,
      n_boot = n_boot,
      conf_level = conf_level,
      progress = progress,
      method = method
    )
    key_cols <- if (is.null(group_col)) c("scope", "n_tokens") else c("scope", "group", "n_tokens")
    wide <- dplyr::left_join(wide, boot, by = key_cols)
    wide <- wide[, c(
      intersect(key_cols, names(wide)),
      intersect(c("n_boot", "conf_level"), names(wide)),
      setdiff(names(wide), c(key_cols, "n_boot", "conf_level"))
    ), drop = FALSE]
  }

  wide <- .select_contrast_columns(wide, metrics)

  if (identical(output, "wide")) {
    return(wide)
  }

  .comparison_long(wide)
}

# ---- metric selection ------------------------------------------------------

.contrast_metric_columns <- function() {
  list(
    jsd           = "jsd",
    js_distance   = "js_distance",
    pillai        = "pillai",
    bhattacharyya = c("bhatt_dist", "bhatt_affinity"),
    mahalanobis   = "mahalanobis_dist",
    overlap       = "percent_overlap"
  )
}

.resolve_contrast_metrics <- function(metrics) {
  choices <- names(.contrast_metric_columns())
  if (is.null(metrics)) {
    return(choices)
  }
  if (!is.character(metrics) || !length(metrics)) {
    stop("`metrics` must be a non-empty character vector.", call. = FALSE)
  }
  metrics <- unique(metrics)
  unknown <- setdiff(metrics, choices)
  if (length(unknown)) {
    stop(
      "Unknown metric(s): ", paste(unknown, collapse = ", "),
      ". Choose from: ", paste(choices, collapse = ", "), ".",
      call. = FALSE
    )
  }
  choices[choices %in% metrics]
}

.select_contrast_columns <- function(wide, metrics) {
  map <- .contrast_metric_columns()
  key_cols <- intersect(
    c("scope", "group", "n_tokens", "n_boot", "conf_level"),
    names(wide)
  )
  metric_cols <- unlist(map[metrics], use.names = FALSE)
  keep <- unlist(lapply(metric_cols, function(col) {
    c(col, paste0(col, c("_n_boot", "_mean", "_sd", "_ci_lower", "_ci_upper")))
  }))
  if ("pillai" %in% metrics) {
    keep <- c(keep, "pillai_p_value")
  }
  keep <- c(key_cols, keep)
  wide[, intersect(names(wide), keep), drop = FALSE]
}

#' Compare phonological category overlap metrics (deprecated)
#'
#' @description
#' `compare_overlap_metrics()` was renamed to [phontrast()] in phontrast 2.0.0
#' (the package formerly released as 'phonJSD'). It remains as a thin wrapper
#' that calls [phontrast()] with `output = "wide"` for backward compatibility
#' and will be removed in a future release. New code should call [phontrast()].
#'
#' @inheritParams phontrast
#' @return See [phontrast()]; wide format by default.
#' @seealso [phontrast()]
#' @keywords internal
#' @export
compare_overlap_metrics <- function(data,
                                    features,
                                    category_col,
                                    group_col = NULL,
                                    min_tokens = 20,
                                    bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                                    eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                                    eval_n = NULL,
                                    eval_seed = NULL,
                                    engine = c("ks", "fast_diag", "fast_diagonal"),
                                    chunk_size = 1000L,
                                    eps = 1e-6,
                                    output = c("wide", "long"),
                                    do_boot = FALSE,
                                    n_boot = 1000,
                                    conf_level = 0.95,
                                    progress = TRUE,
                                    method = c("mc", "legacy")) {
  .Deprecated("phontrast")
  output <- match.arg(output)
  phontrast(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    min_tokens = min_tokens,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size,
    eps = eps,
    output = output,
    do_boot = do_boot,
    n_boot = n_boot,
    conf_level = conf_level,
    progress = progress,
    method = method
  )
}

.compare_overlap_metrics_point <- function(data,
                                           features,
                                           category_col,
                                           group_col = NULL,
                                           min_tokens = 20,
                                           bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                                           eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                                           eval_n = NULL,
                                           eval_seed = NULL,
                                           engine = c("ks", "fast_diag", "fast_diagonal"),
                                           chunk_size = 1000L,
                                           eps = 1e-6,
                                           method = c("mc", "legacy")) {
  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  engine <- .match_kde_engine(engine)
  method <- match.arg(method)

  jsd_out <- estimate_jsd(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    do_boot = FALSE,
    min_tokens = min_tokens,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size,
    method = method
  )
  jsd_wide <- jsd_out[, intersect(c("scope", "group", "n_tokens"), names(jsd_out)), drop = FALSE]
  jsd_wide$jsd <- jsd_out$jsd_point
  jsd_wide$js_distance <- sqrt(jsd_out$jsd_point)

  pillai_out <- estimate_pillai(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    min_tokens = min_tokens
  )
  pillai_wide <- pillai_out[, intersect(c("scope", "group", "n_tokens"), names(pillai_out)), drop = FALSE]
  pillai_wide$pillai <- pillai_out$pillai
  pillai_wide$pillai_p_value <- pillai_out$p_value

  bhatt_out <- estimate_bhatt(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    min_tokens = min_tokens,
    eps = eps
  )
  bhatt_wide <- bhatt_out[, intersect(c("scope", "group", "n_tokens"), names(bhatt_out)), drop = FALSE]
  bhatt_wide$bhatt_dist <- bhatt_out$bhatt_dist
  bhatt_wide$bhatt_affinity <- bhatt_out$bhatt_affinity

  mahal_wide <- .estimate_mahalanobis(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    min_tokens = min_tokens,
    eps = eps
  )

  overlap_out <- estimate_overlap(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    min_tokens = min_tokens,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size,
    method = method
  )
  overlap_wide <- overlap_out[, intersect(c("scope", "group", "n_tokens"), names(overlap_out)), drop = FALSE]
  overlap_wide$percent_overlap <- overlap_out$overlap

  key_cols <- if (is.null(group_col)) c("scope", "n_tokens") else c("scope", "group", "n_tokens")
  pieces <- list(
    pillai_wide,
    bhatt_wide,
    jsd_wide,
    mahal_wide,
    overlap_wide
  )
  if (!is.null(group_col)) {
    pieces <- lapply(pieces, function(piece) {
      if ("group" %in% names(piece)) {
        piece$group <- as.character(piece$group)
      }
      piece
    })
  }
  wide <- Reduce(function(x, y) {
    dplyr::full_join(x, y, by = intersect(key_cols, intersect(names(x), names(y))))
  }, pieces)
  wide <- wide[, c(intersect(key_cols, names(wide)), setdiff(names(wide), key_cols)), drop = FALSE]
  wide
}

.compare_metric_columns <- function() {
  c(
    "pillai",
    "bhatt_dist",
    "bhatt_affinity",
    "jsd",
    "js_distance",
    "mahalanobis_dist",
    "percent_overlap"
  )
}

.warn_empty_overlap_comparison <- function(data,
                                           features,
                                           category_col,
                                           group_col = NULL,
                                           min_tokens = 20) {
  keep_cols <- if (is.null(group_col)) c(category_col, features) else c(group_col, category_col, features)
  df <- tryCatch(.metric_data(data, keep_cols), error = function(e) NULL)
  min_per_category <- .kde_min_category_tokens(length(features))

  if (is.null(df) || !nrow(df)) {
    warning(
      "phontrast() returned no rows after removing missing or non-finite values.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  if (is.null(group_col)) {
    counts <- .observed_category_counts(df[[category_col]])
    warning(
      "phontrast() returned no rows. Observed category counts after filtering were: ",
      paste(names(counts), as.integer(counts), sep = "=", collapse = ", "),
      ". Exactly two observed categories and at least ", min_per_category,
      " observations per category are required for KDE-based metrics.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  groups <- .split_groups(df, group_col)
  group_n <- vapply(groups, nrow, integer(1))
  category_counts <- lapply(groups, function(df_g) .observed_category_counts(df_g[[category_col]]))
  exactly_two <- vapply(category_counts, length, integer(1)) == 2L
  meets_min <- group_n >= min_tokens
  meets_kde <- vapply(
    category_counts,
    function(x) length(x) == 2L && all(x >= min_per_category),
    logical(1)
  )
  min_category_counts <- vapply(
    category_counts,
    function(x) if (length(x)) min(as.integer(x)) else 0L,
    integer(1)
  )

  warning(
    "phontrast() returned no grouped rows. After removing missing/non-finite values, ",
    sum(meets_min & exactly_two), " of ", length(groups),
    " groups had at least min_tokens = ", min_tokens,
    " and exactly two observed categories; ",
    sum(meets_min & meets_kde), " also had at least ", min_per_category,
    " observations per category for KDE metrics. Max group size was ",
    max(group_n), "; largest within-group minimum category count was ",
    max(min_category_counts), ". If you intended a global contrast, omit `group_col`; ",
    "otherwise group at a coarser level or use a feature space with enough observations per category.",
    call. = FALSE
  )
  invisible(NULL)
}

.bootstrap_compare_overlap_metrics <- function(data,
                                               point_wide,
                                               features,
                                               category_col,
                                               group_col = NULL,
                                               min_tokens = 20,
                                               bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                                               eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                                               eval_n = NULL,
                                               eval_seed = NULL,
                                               engine = c("ks", "fast_diag", "fast_diagonal"),
                                               chunk_size = 1000L,
                                               eps = 1e-6,
                                               n_boot = 300,
                                               conf_level = 0.95,
                                               progress = TRUE,
                                               method = "mc") {
  key_cols <- if (is.null(group_col)) c("scope", "n_tokens") else c("scope", "group", "n_tokens")

  if (is.null(group_col)) {
    df <- .metric_data(data, c(category_col, features))
    boot_rows <- list(.bootstrap_one_overlap_source(
      df = df,
      label = "global comparison",
      features = features,
      category_col = category_col,
      min_tokens = min_tokens,
      bw = bw,
      eval_on = eval_on,
      eval_n = eval_n,
      eval_seed = eval_seed,
      engine = engine,
      chunk_size = chunk_size,
      eps = eps,
      n_boot = n_boot,
      conf_level = conf_level,
      progress = progress,
      method = method
    ))
    out <- cbind(point_wide[, key_cols, drop = FALSE], dplyr::bind_rows(boot_rows))
    rownames(out) <- NULL
    return(out)
  }

  df <- .metric_data(data, c(group_col, category_col, features))
  groups <- .split_groups(df, group_col)
  boot_rows <- lapply(seq_len(nrow(point_wide)), function(i) {
    group_id <- as.character(point_wide$group[i])
    df_g <- groups[[group_id]]
    if (is.null(df_g)) {
      return(.empty_boot_overlap_summary(n_boot, conf_level))
    }
    .bootstrap_one_overlap_source(
      df = df_g,
      label = paste0("group `", group_id, "`"),
      features = features,
      category_col = category_col,
      min_tokens = min_tokens,
      bw = bw,
      eval_on = eval_on,
      eval_n = eval_n,
      eval_seed = eval_seed,
      engine = engine,
      chunk_size = chunk_size,
      eps = eps,
      n_boot = n_boot,
      conf_level = conf_level,
      progress = progress,
      method = method
    )
  })
  out <- cbind(point_wide[, key_cols, drop = FALSE], dplyr::bind_rows(boot_rows))
  rownames(out) <- NULL
  out
}

.bootstrap_one_overlap_source <- function(df,
                                          label,
                                          features,
                                          category_col,
                                          min_tokens = 20,
                                          bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                                          eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                                          eval_n = NULL,
                                          eval_seed = NULL,
                                          engine = c("ks", "fast_diag", "fast_diagonal"),
                                          chunk_size = 1000L,
                                          eps = 1e-6,
                                          n_boot = 300,
                                          conf_level = 0.95,
                                          progress = TRUE,
                                          method = "mc") {
  if (isTRUE(progress)) {
    message(
      "Bootstrapping overlap metrics for ", label, " (",
      n_boot, " resamples; conf_level = ", conf_level,
      "). This may take time."
    )
  }

  metric_cols <- .compare_metric_columns()
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(metric_cols))
  colnames(boot_mat) <- metric_cols
  n <- nrow(df)
  progress_every <- max(1L, floor(n_boot / 10))

  for (b in seq_len(n_boot)) {
    if (isTRUE(progress) && (b == 1L || b == n_boot || b %% progress_every == 0L)) {
      message("  ", label, ": bootstrap replicate ", b, " / ", n_boot)
    }

    samp <- df[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
    if (.observed_n_categories(samp[[category_col]]) != 2L) {
      next
    }
    vals <- tryCatch(
      .compare_overlap_metrics_point(
        data = samp,
        features = features,
        category_col = category_col,
        group_col = NULL,
        min_tokens = min_tokens,
        bw = bw,
        eval_on = eval_on,
        eval_n = eval_n,
        eval_seed = eval_seed,
        engine = engine,
        chunk_size = chunk_size,
        eps = eps,
        method = method
      ),
      error = function(e) NULL
    )
    if (is.null(vals) || !nrow(vals)) {
      next
    }
    boot_mat[b, metric_cols] <- as.numeric(vals[1, metric_cols, drop = TRUE])
  }

  .summarize_boot_overlap_metrics(boot_mat, n_boot, conf_level)
}

.summarize_boot_overlap_metrics <- function(boot_mat, n_boot, conf_level) {
  alpha <- 1 - conf_level
  out <- data.frame(
    n_boot = n_boot,
    conf_level = conf_level,
    stringsAsFactors = FALSE
  )

  for (metric in colnames(boot_mat)) {
    vals <- boot_mat[, metric]
    vals <- vals[is.finite(vals)]
    out[[paste0(metric, "_n_boot")]] <- length(vals)
    out[[paste0(metric, "_mean")]] <- if (length(vals)) mean(vals) else NA_real_
    out[[paste0(metric, "_sd")]] <- if (length(vals) > 1L) stats::sd(vals) else NA_real_
    if (length(vals)) {
      qs <- stats::quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
      out[[paste0(metric, "_ci_lower")]] <- qs[1]
      out[[paste0(metric, "_ci_upper")]] <- qs[2]
    } else {
      out[[paste0(metric, "_ci_lower")]] <- NA_real_
      out[[paste0(metric, "_ci_upper")]] <- NA_real_
    }
  }

  out
}

.empty_boot_overlap_summary <- function(n_boot, conf_level) {
  boot_mat <- matrix(
    NA_real_,
    nrow = n_boot,
    ncol = length(.compare_metric_columns()),
    dimnames = list(NULL, .compare_metric_columns())
  )
  .summarize_boot_overlap_metrics(boot_mat, n_boot, conf_level)
}

.mahalanobis_distance <- function(data, features, category_col, eps = 1e-6) {
  .check_columns(data, c(category_col, features))
  data <- .metric_data(data, c(category_col, features))
  .check_numeric_features(data, features)
  .check_ridge_eps(eps, "eps")
  levs <- .two_levels(data[[category_col]], "category_col")
  .check_two_category_sample_size(
    data,
    category_col,
    .kde_min_category_tokens(length(features)),
    "Mahalanobis distance"
  )

  X1 <- as.matrix(data[data[[category_col]] == levs[1], features, drop = FALSE])
  X2 <- as.matrix(data[data[[category_col]] == levs[2], features, drop = FALSE])
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  pooled_cov <- ((n1 - 1) * stats::cov(X1) + (n2 - 1) * stats::cov(X2)) /
    (n1 + n2 - 2)
  pooled_cov <- pooled_cov + diag(eps, ncol(pooled_cov))

  inv_cov <- tryCatch(solve(pooled_cov), error = function(e) NULL)
  if (is.null(inv_cov)) {
    stop(
      "Mahalanobis distance failed: pooled covariance is not positive definite. ",
      "Try increasing `eps` or reducing feature dimensionality.",
      call. = FALSE
    )
  }

  diff <- matrix(colMeans(X2) - colMeans(X1), ncol = 1)
  sqrt(as.numeric(t(diff) %*% inv_cov %*% diff))
}

.estimate_mahalanobis <- function(data,
                                  features,
                                  category_col,
                                  group_col = NULL,
                                  min_tokens = 20,
                                  eps = 1e-6) {
  .check_positive_count(min_tokens, "min_tokens")

  if (is.null(group_col)) {
    df <- .metric_data(data, c(category_col, features))
    n <- nrow(df)
    if (n < min_tokens) {
      stop("Not enough tokens after removing missing values. Got ",
           n, ", need at least ", min_tokens, ".")
    }
    return(data.frame(
      scope = "global",
      n_tokens = n,
      mahalanobis_dist = .mahalanobis_distance(df, features, category_col, eps = eps),
      stringsAsFactors = FALSE
    ))
  }

  group_col <- .check_group_cols(group_col)
  .check_columns(data, c(group_col, category_col, features))
  df <- .metric_data(data, c(group_col, category_col, features))
  out <- lapply(.split_groups(df, group_col), function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens || .observed_n_categories(df_g[[category_col]]) != 2L) {
      return(NULL)
    }
    dist <- tryCatch(
      .mahalanobis_distance(df_g, features, category_col, eps = eps),
      error = function(e) NA_real_
    )
    data.frame(
      scope = "group",
      group = .group_label(df_g, group_col),
      n_tokens = n_tok,
      mahalanobis_dist = dist,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(
      scope = character(),
      group = character(),
      n_tokens = integer(),
      mahalanobis_dist = numeric(),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  .warn_failed_groups(out, "mahalanobis_dist", "estimate_mahalanobis()")
}

.comparison_long <- function(wide) {
  specs <- list(
    list(column = "pillai", metric = "Pillai trace", orientation = "separation",
         bounded = TRUE, transform = identity),
    list(column = "bhatt_dist", metric = "Bhattacharyya distance",
         orientation = "separation", bounded = FALSE, transform = identity),
    list(column = "bhatt_affinity", metric = "Bhattacharyya affinity",
         orientation = "overlap", bounded = TRUE, transform = function(x) 1 - x),
    list(column = "jsd", metric = "Jensen-Shannon divergence",
         orientation = "separation", bounded = TRUE, transform = identity),
    list(column = "js_distance", metric = "Jensen-Shannon distance",
         orientation = "separation", bounded = TRUE, transform = identity),
    list(column = "mahalanobis_dist", metric = "Mahalanobis distance",
         orientation = "separation", bounded = FALSE, transform = identity),
    list(column = "percent_overlap", metric = "Percent overlap",
         orientation = "overlap", bounded = TRUE, transform = function(x) 1 - x)
  )
  key_cols <- intersect(c("scope", "group", "n_tokens"), names(wide))
  rows <- lapply(specs, function(spec) {
    if (!spec$column %in% names(wide)) {
      return(NULL)
    }
    out <- wide[, key_cols, drop = FALSE]
    out$metric <- spec$metric
    out$estimate <- wide[[spec$column]]
    out$orientation <- spec$orientation
    out$bounded_0_1 <- spec$bounded
    out$separation_value <- spec$transform(out$estimate)
    out$p_value <- if (identical(spec$column, "pillai") && "pillai_p_value" %in% names(wide)) {
      wide$pillai_p_value
    } else {
      NA_real_
    }
    boot_cols <- paste0(
      spec$column,
      c("_n_boot", "_mean", "_sd", "_ci_lower", "_ci_upper")
    )
    if (all(boot_cols %in% names(wide))) {
      out$n_boot <- wide[[boot_cols[1]]]
      out$conf_level <- wide$conf_level
      out$boot_mean <- wide[[boot_cols[2]]]
      out$boot_sd <- wide[[boot_cols[3]]]
      out$ci_lower <- wide[[boot_cols[4]]]
      out$ci_upper <- wide[[boot_cols[5]]]
    }
    out
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) {
    return(out)
  }
  out$separation_rank <- stats::ave(
    out$separation_value,
    out$metric,
    FUN = function(x) rank(-x, ties.method = "average", na.last = "keep")
  )
  out
}
