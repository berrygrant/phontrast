#' Compare major phonological category overlap metrics
#'
#' Computes the package's main category-separation and overlap metrics in one
#' call: Pillai trace, Bhattacharyya distance and affinity, Jensen-Shannon
#' divergence, Jensen-Shannon distance, Mahalanobis distance, and percent
#' overlap. Results can be returned in a wide format for analysis tables or a
#' long format for plotting and rank-based comparison.
#'
#' This is the recommended entry point when you want to compare more than one
#' overlap metric for the same phonological contrast. Use \code{estimate_jsd()}
#' when JSD is the only outcome of interest, and use lower-level metric helpers
#' only when you need direct control over one estimator.
#'
#' Metric directions differ. JSD, Jensen-Shannon distance, Pillai trace,
#' Bhattacharyya distance, and Mahalanobis distance increase as categories
#' become more separated. Percent overlap and Bhattacharyya affinity increase
#' as categories overlap more. Long output includes \code{orientation},
#' \code{separation_value}, and \code{separation_rank} columns so all metrics
#' can be read on a separation-oriented scale.
#'
#' @param data Data frame containing category labels and acoustic features.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; column giving the two categories to compare.
#' @param group_col Optional string; grouping column. If \code{NULL}, metrics
#'   are computed globally.
#' @param min_tokens Minimum tokens required globally or per group.
#' @param bw Bandwidth selection method passed to \code{jsd_kde_nd()} and
#'   \code{percent_overlap_kde()}.
#' @param eval_on KDE evaluation points passed to \code{jsd_kde_nd()} and
#'   \code{percent_overlap_kde()}.
#' @param eps Small ridge constant for covariance-based metrics.
#' @param output Output format: \code{"wide"} returns one row per global/group
#'   comparison; \code{"long"} returns one row per metric per comparison.
#'
#' @return A data frame. Wide output contains one column per metric. Long output
#'   contains \code{metric}, \code{estimate}, \code{orientation},
#'   \code{bounded_0_1}, \code{separation_value}, and
#'   \code{separation_rank} columns.
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
#' compare_overlap_metrics(
#'   data = vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   group_col = "speaker",
#'   output = "wide"
#' )
#'
#' metrics_long <- compare_overlap_metrics(
#'   data = vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   group_col = "speaker",
#'   output = "long"
#' )
#'
#' metrics_long[, c("group", "metric", "estimate", "orientation",
#'                  "separation_value", "separation_rank")]
#' @export
compare_overlap_metrics <- function(data,
                                    features,
                                    category_col,
                                    group_col = NULL,
                                    min_tokens = 20,
                                    bw = c("Hpi", "Hscv", "Hpi.diag"),
                                    eval_on = c("pooled", "group1", "group2"),
                                    eps = 1e-6,
                                    output = c("wide", "long")) {
  output <- match.arg(output)
  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  .check_positive_count(min_tokens, "min_tokens")
  .check_ridge_eps(eps, "eps")

  jsd_out <- estimate_jsd(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    do_boot = FALSE,
    min_tokens = min_tokens,
    bw = bw,
    eval_on = eval_on
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
    eval_on = eval_on
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
  wide <- Reduce(function(x, y) {
    dplyr::full_join(x, y, by = intersect(key_cols, intersect(names(x), names(y))))
  }, pieces)
  wide <- wide[, c(intersect(key_cols, names(wide)), setdiff(names(wide), key_cols)), drop = FALSE]

  if (identical(output, "wide")) {
    return(wide)
  }

  .comparison_long(wide)
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

  .check_columns(data, c(group_col, category_col, features))
  df <- .metric_data(data, c(group_col, category_col, features))
  out <- lapply(split(df, df[[group_col]]), function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens || length(unique(df_g[[category_col]])) != 2L) {
      return(NULL)
    }
    dist <- tryCatch(
      .mahalanobis_distance(df_g, features, category_col, eps = eps),
      error = function(e) NA_real_
    )
    data.frame(
      scope = "group",
      group = df_g[[group_col]][1],
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
  out
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
    out
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) {
    return(out)
  }
  out$separation_rank <- ave(
    out$separation_value,
    out$metric,
    FUN = function(x) rank(-x, ties.method = "average", na.last = "keep")
  )
  out
}
