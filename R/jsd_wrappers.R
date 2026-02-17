#' Estimate Jensen–Shannon divergence or distance between two categories
#'
#' @param data Data frame with at least `category_col` and `features`.
#' @param features Character vector of feature column names (e.g., c("F1","F2")).
#' @param category_col Name of the column giving the two-way category factor.
#' @param group_col Optional grouping column. If provided, returns per-group JSD.
#' @param do_boot Logical; if TRUE, run nonparametric bootstrap.
#' @param n_boot Number of bootstrap resamples.
#' @param min_tokens Minimum total tokens required (globally or per group).
#' @param est_distance Logical; if TRUE, return Jensen–Shannon *distance* (sqrt of divergence).
#' @param conf_level Confidence level for bootstrap interval.
#' @param ... Additional arguments passed to \code{jsd_kde_nd()}.
#'
#' @return A tibble. Global: one row with columns
#'   scope, n_tokens, n_boot, jsd_point, jsd_mean, jsd_sd, jsd_low, jsd_high.
#'   Grouped: one row per group with columns
#'   scope, group, n_tokens, n_boot, jsd_point, jsd_mean, jsd_sd, jsd_low, jsd_high.
#'
#' @export
#' @importFrom rlang .data
estimate_jsd <- function(data,
                         features,
                         category_col,
                         group_col    = NULL,
                         do_boot      = FALSE,
                         n_boot       = 1000,
                         min_tokens   = 5,
                         est_distance = FALSE,
                         conf_level   = 0.95,
                         ...) {
  
  # Basic checks
  if (!all(features %in% names(data))) {
    stop("estimate_jsd(): All `features` must be columns in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("estimate_jsd(): `category_col` must be a column in `data`.")
  }
  
  if (is.null(group_col)) {
    # ---- Global ----
    df <- data[, c(category_col, features), drop = FALSE]
    df <- df[stats::complete.cases(df), , drop = FALSE]
    
    n <- nrow(df)
    if (n < min_tokens) {
      stop("estimate_jsd(): Not enough tokens after removing missing values. Got ",
           n, ", need at least ", min_tokens, ".")
    }
    
    # Ensure exactly two categories
    g <- droplevels(factor(df[[category_col]]))
    if (nlevels(g) != 2L) {
      stop("estimate_jsd(): `category_col` must have exactly 2 levels in the filtered data.")
    }
    
    # ---- Point estimate via KDE ----
    jsd_div_point <- jsd_kde_nd(
      data     = df,
      features = features,
      group    = category_col  # pass column name, not vector
    )
    
    jsd_point <- if (est_distance) sqrt(jsd_div_point) else jsd_div_point
    
    # If no bootstrap, just return the point estimate row
    if (!do_boot) {
      return(tibble::tibble(
        scope     = "global",
        n_tokens  = n,
        n_boot    = 0L,
        jsd_point = jsd_point,
        jsd_mean  = NA_real_,
        jsd_sd    = NA_real_,
        jsd_low   = NA_real_,
        jsd_high  = NA_real_
      ))
    }
    
    # ---- Bootstrap case ----
    alpha <- 1 - conf_level
    
    boot_vals <- replicate(n_boot, {
      idx <- sample.int(n, replace = TRUE)
      df_boot <- df[idx, , drop = FALSE]
      
      # Guard: some bootstrap samples may lose one category
      if (length(unique(df_boot[[category_col]])) < 2L) {
        return(NA_real_)
      }
      
      jsd_div_boot <- tryCatch(
        jsd_kde_nd(
          data     = df_boot,
          features = features,
          group    = category_col
        ),
        error = function(e) NA_real_
      )
      
      if (is.na(jsd_div_boot)) {
        return(NA_real_)
      }
      
      if (est_distance) sqrt(jsd_div_boot) else jsd_div_boot
    })
    
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
      jsd_point = jsd_point,
      jsd_mean  = jsd_mean,
      jsd_sd    = jsd_sd,
      jsd_low   = qs[1],
      jsd_high  = qs[2]
    ))
  }
  
  # ---- Grouped ----
  if (!group_col %in% names(data)) {
    stop("estimate_jsd(): `group_col` must be a column in `data`.")
  }
  
  keep_cols <- c(group_col, category_col, features)
  df <- data[, keep_cols, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  
  pt <- speaker_jsd(
    data         = df,
    group_col    = group_col,
    category_col = category_col,
    features     = features,
    min_tokens   = min_tokens,
    ...
  ) |>
    dplyr::rename(jsd_point = .data$jsd)
  
  if (est_distance) {
    pt$jsd_point <- sqrt(pt$jsd_point)
  }
  
  if (!do_boot) {
    out <- pt
    out$scope <- "group"
    out$n_boot <- 0L
    out$jsd_mean <- NA_real_
    out$jsd_sd   <- NA_real_
    out$jsd_low  <- NA_real_
    out$jsd_high <- NA_real_
    out <- out[, c("scope", "group", "n_tokens", "n_boot",
                   "jsd_point", "jsd_mean", "jsd_sd", "jsd_low", "jsd_high")]
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
    ...
  )
  
  out <- dplyr::left_join(pt, bt, by = c("group", "n_tokens"))
  out$scope <- "group"
  out$n_boot <- as.integer(n_boot)
  out <- out[, c("scope", "group", "n_tokens", "n_boot",
                 "jsd_point", "jsd_mean", "jsd_sd", "jsd_low", "jsd_high")]
  out
}
