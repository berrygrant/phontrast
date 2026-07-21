#' JSD summary: point estimate and optional bootstrap per group
#'
#' Convenience wrapper that returns both the point-estimate JSD and,
#' optionally, bootstrap-based uncertainty (mean, SD, and CI) for each group.
#'
#' @inheritParams speaker_jsd
#' @param do_boot Logical; if TRUE (default), perform bootstrap via \code{boot_jsd()}.
#' @param n_boot Integer; number of bootstrap resamples per group if
#'   \code{do_boot = TRUE}.
#' @param conf_level Confidence level for bootstrap intervals.
#' @param bw Bandwidth selection method passed to \code{jsd_kde_nd()}.
#' @param eval_on KDE evaluation points passed to \code{jsd_kde_nd()}.
#' @param eval_n Optional maximum number of KDE evaluation points.
#' @param eval_seed Optional integer seed for KDE evaluation-point subsampling.
#' @param engine KDE evaluation engine passed to \code{jsd_kde_nd()}.
#'   \code{"fast_diagonal"} is accepted as an alias for \code{"fast_diag"}.
#' @param chunk_size Chunk size for \code{engine = "fast_diag"}.
#' @param method Estimator passed to \code{jsd_kde_nd()}: \code{"mc"} (default)
#'   or \code{"legacy"} (pre-1.2.0 self-normalized estimate).
#'
#' @return A tibble with one row per group and columns:
#'   \itemize{
#'     \item \code{group} - group ID (e.g., speaker)
#'     \item \code{n_tokens} - number of tokens for that group
#'     \item \code{jsd_point} - single JSD point estimate
#'     \item \code{n_boot}, \code{conf_level}, \code{jsd_mean},
#'       \code{jsd_sd}, \code{ci_lower}, \code{ci_upper},
#'       \code{jsd_low}, \code{jsd_high} - bootstrap summary columns.
#'       These are \code{NA} (or 0 for \code{n_boot}) if
#'       \code{do_boot = FALSE}.
#'   }
#' @export
#' @importFrom dplyr left_join rename
#' @importFrom rlang .data
jsd_summary <- function(data,
                        group_col,
                        category_col,
                        features,
                        do_boot     = TRUE,
                        n_boot      = 1000,
                        min_tokens  = 20,
                        conf_level  = 0.95,
                        bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                        eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                        eval_n = NULL,
                        eval_seed = NULL,
                        engine = c("ks", "fast_diag", "fast_diagonal"),
                        chunk_size = 1000L,
                        method = c("mc", "legacy"),
                        ...) {

  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  engine <- .match_kde_engine(engine)
  method <- match.arg(method)
  .check_conf_level(conf_level)
  if (isTRUE(do_boot)) {
    .check_positive_count(n_boot, "n_boot")
  }
  .check_positive_count(min_tokens, "min_tokens")
  .validate_metric_inputs(data, features, category_col, group_col)

  # Point estimates
  pt <- speaker_jsd(
    data         = data,
    group_col    = group_col,
    category_col = category_col,
    features     = features,
    min_tokens   = min_tokens,
    bw           = bw,
    eval_on      = eval_on,
    eval_n       = eval_n,
    eval_seed    = eval_seed,
    engine       = engine,
    chunk_size   = chunk_size,
    method       = method,
    ...
  ) |>
    dplyr::rename(jsd_point = "jsd")

  if (!do_boot) {
    pt$n_boot <- 0L
    pt$conf_level <- conf_level
    pt$jsd_mean <- NA_real_
    pt$jsd_sd <- NA_real_
    pt$ci_lower <- NA_real_
    pt$ci_upper <- NA_real_
    pt$jsd_low <- NA_real_
    pt$jsd_high <- NA_real_
    return(pt)
  }

  # Bootstrap estimates
  bt <- boot_jsd(
    data         = data,
    group_col    = group_col,
    category_col = category_col,
    features     = features,
    n_boot       = n_boot,
    min_tokens   = min_tokens,
    conf_level   = conf_level,
    bw           = bw,
    eval_on      = eval_on,
    eval_n       = eval_n,
    eval_seed    = eval_seed,
    engine       = engine,
    chunk_size   = chunk_size,
    method       = method,
    ...
  )

  # Join on group + n_tokens (both functions report them)
  out <- dplyr::left_join(
    pt,
    bt,
    by = c("group", "n_tokens")
  )

  out
}
