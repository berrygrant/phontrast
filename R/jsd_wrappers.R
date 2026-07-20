#' Estimate Jensen–Shannon divergence or distance between two categories
#'
#' Use this function when Jensen-Shannon divergence (JSD) or Jensen-Shannon
#' distance is the primary outcome. If you want to compare JSD with Pillai,
#' Bhattacharyya, Mahalanobis, and percent-overlap metrics, start with
#' \code{compare_overlap_metrics()} instead.
#'
#' JSD is bounded from 0 to 1 for two equally weighted distributions. Larger
#' values indicate greater category separation. Jensen-Shannon distance is
#' \code{sqrt(JSD)} and has the same direction of interpretation.
#'
#' @param data Data frame with at least `category_col` and `features`.
#' @param features Character vector of feature column names (e.g., c("F1","F2")).
#' @param category_col Name of the column giving the two-way category factor.
#' @param group_col Optional character vector of one or more grouping columns.
#'   If provided, returns per-group JSD. Multiple grouping columns are combined
#'   into a labeled \code{group} value such as \code{"Sex=F | Style=read"}.
#' @param do_boot Logical; if TRUE, run nonparametric bootstrap.
#' @param n_boot Number of bootstrap resamples.
#' @param min_tokens Minimum total tokens required (globally or per group).
#' @param est_distance Logical; if TRUE, return Jensen–Shannon *distance* (sqrt of divergence).
#' @param conf_level Confidence level for bootstrap interval.
#' @param bw Bandwidth selection method passed to \code{jsd_kde_nd()}.
#' @param eval_on KDE evaluation points passed to \code{jsd_kde_nd()}.
#' @param eval_n Optional maximum number of KDE evaluation points.
#' @param eval_seed Optional integer seed for KDE evaluation-point subsampling.
#' @param engine KDE evaluation engine passed to \code{jsd_kde_nd()}.
#'   \code{"fast_diagonal"} is accepted as an alias for \code{"fast_diag"}.
#' @param chunk_size Chunk size for \code{engine = "fast_diag"}.
#' @param method Estimator passed to \code{jsd_kde_nd()}: \code{"mc"} (default)
#'   for the Monte-Carlo plug-in estimate of the continuous JSD, or
#'   \code{"legacy"} to reproduce the pre-1.2.0 self-normalized sample-point
#'   estimate.
#' @param ... Additional arguments passed to \code{jsd_kde_nd()} (e.g.,
#'   \code{loo}).
#'
#' @return A tibble. Global: one row with columns
#'   scope, n_tokens, n_boot, conf_level, jsd_point, jsd_mean, jsd_sd,
#'   ci_lower, ci_upper, jsd_low, jsd_high.
#'   Grouped: one row per group with columns
#'   scope, group, n_tokens, n_boot, conf_level, jsd_point, jsd_mean,
#'   jsd_sd, ci_lower, ci_upper, jsd_low, jsd_high.
#'   \code{ci_lower} and \code{ci_upper} are the preferred confidence interval
#'   columns; \code{jsd_low} and \code{jsd_high} are retained as legacy aliases.
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
#' # Global JSD with a small bootstrap interval.
#' # Increase n_boot for real analyses.
#' estimate_jsd(
#'   data = vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   do_boot = TRUE,
#'   n_boot = 5,
#'   min_tokens = 20
#' )
#'
#' # Grouped JSD by speaker, again with a small example bootstrap.
#' estimate_jsd(
#'   data = vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   group_col = "speaker",
#'   do_boot = TRUE,
#'   n_boot = 5,
#'   min_tokens = 20
#' )
#' @export
#' @importFrom rlang .data
estimate_jsd <- function(data,
                         features,
                         category_col,
                         group_col    = NULL,
                         do_boot      = FALSE,
                         n_boot       = 1000,
                         min_tokens   = 20,
                         est_distance = FALSE,
                         conf_level   = 0.95,
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
  .check_columns(data, c(category_col, features))

  if (is.null(group_col)) {
    # ---- Global ----
    df <- .metric_data(data, c(category_col, features))
    
    n <- nrow(df)
    if (n < min_tokens) {
      stop("estimate_jsd(): Not enough tokens after removing missing values. Got ",
           n, ", need at least ", min_tokens, ".")
    }
    
    # Ensure exactly two categories
    .two_levels(df[[category_col]], "category_col")
    
    # ---- Point estimate via KDE ----
    jsd_div_point <- jsd_kde_nd(
      data     = df,
      features = features,
      group    = category_col,
      bw        = bw,
      eval_on   = eval_on,
      eval_n    = eval_n,
      eval_seed = eval_seed,
      engine    = engine,
      chunk_size = chunk_size,
      method    = method,
      ...
    )

    jsd_point <- if (est_distance) sqrt(jsd_div_point) else jsd_div_point
    
    # If no bootstrap, just return the point estimate row
    if (!do_boot) {
      return(tibble::tibble(
        scope     = "global",
        n_tokens  = n,
        n_boot    = 0L,
        conf_level = conf_level,
        jsd_point = jsd_point,
        jsd_mean  = NA_real_,
        jsd_sd    = NA_real_,
        ci_lower  = NA_real_,
        ci_upper  = NA_real_,
        jsd_low   = NA_real_,
        jsd_high  = NA_real_
      ))
    }
    
    # ---- Bootstrap case ----
    alpha <- 1 - conf_level
    
    boot_vals <- vapply(seq_len(n_boot), function(i) {
      idx <- sample.int(n, replace = TRUE)
      df_boot <- df[idx, , drop = FALSE]
      
      # Guard: some bootstrap samples may lose one category
      if (.observed_n_categories(df_boot[[category_col]]) != 2L) {
        return(NA_real_)
      }
      
      jsd_div_boot <- tryCatch(
        jsd_kde_nd(
          data     = df_boot,
          features = features,
          group    = category_col,
          bw        = bw,
          eval_on   = eval_on,
          eval_n    = eval_n,
          eval_seed = eval_seed,
          engine    = engine,
          chunk_size = chunk_size,
          method    = method,
          ...
        ),
        error = function(e) NA_real_
      )
      
      if (is.na(jsd_div_boot)) {
        return(NA_real_)
      }
      
      if (est_distance) sqrt(jsd_div_boot) else jsd_div_boot
    }, numeric(1))
    
    boot_vals <- boot_vals[!is.na(boot_vals)]
    if (!length(boot_vals)) {
      stop("estimate_jsd(): All bootstrap samples failed. Try increasing n or reducing n_boot.")
    }
    
    jsd_mean <- mean(boot_vals)
    jsd_sd   <- stats::sd(boot_vals)
    qs       <- stats::quantile(
      boot_vals,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE
    )
    
    return(tibble::tibble(
      scope     = "global",
      n_tokens  = n,
      n_boot    = as.integer(length(boot_vals)),
      conf_level = conf_level,
      jsd_point = jsd_point,
      jsd_mean  = jsd_mean,
      jsd_sd    = jsd_sd,
      ci_lower  = qs[1],
      ci_upper  = qs[2],
      jsd_low   = qs[1],
      jsd_high  = qs[2]
    ))
  }
  
  # ---- Grouped ----
  group_col <- .check_group_cols(group_col)
  .check_columns(data, group_col, "group_col")

  keep_cols <- c(group_col, category_col, features)
  df <- .metric_data(data, keep_cols)
  
  pt <- speaker_jsd(
    data         = df,
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
  
  if (est_distance) {
    pt$jsd_point <- sqrt(pt$jsd_point)
  }
  
  if (!do_boot) {
    out <- pt
    out$scope <- "group"
    out$n_boot <- 0L
    out$conf_level <- conf_level
    out$jsd_mean <- NA_real_
    out$jsd_sd   <- NA_real_
    out$ci_lower <- NA_real_
    out$ci_upper <- NA_real_
    out$jsd_low  <- NA_real_
    out$jsd_high <- NA_real_
    out <- out[, c("scope", "group", "n_tokens", "n_boot",
                   "conf_level", "jsd_point", "jsd_mean", "jsd_sd",
                   "ci_lower", "ci_upper", "jsd_low", "jsd_high")]
    return(out)
  }
  
  bt <- boot_jsd(
    data         = df,
    group_col    = group_col,
    category_col = category_col,
    features     = features,
    n_boot       = n_boot,
    min_tokens   = min_tokens,
    est_distance = est_distance,
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

  out <- dplyr::left_join(pt, bt, by = c("group", "n_tokens"))
  out$scope <- "group"
  out <- out[, c("scope", "group", "n_tokens", "n_boot",
                 "conf_level", "jsd_point", "jsd_mean", "jsd_sd",
                 "ci_lower", "ci_upper", "jsd_low", "jsd_high")]
  out
}
