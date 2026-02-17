#' Pillai trace for multivariate overlap
#'
#' Computes the Pillai-Bartlett trace from a MANOVA of features ~ category.
#' This is a convenience wrapper for comparison with JSD.
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; column giving categories (>= 2).
#'
#' @return A list with elements `pillai` and `p_value`.
#' @export
#' @importFrom stats manova cov
pillai_overlap <- function(data, features, category_col) {
  if (!category_col %in% names(data)) {
    stop("`category_col` must be in `data`.")
  }
  Y <- as.matrix(data[, features, drop = FALSE])
  cat <- data[[category_col]]
  m <- stats::manova(Y ~ cat)
  s <- summary(m, test = "Pillai")
  list(
    pillai  = s$stats[1, "Pillai"],
    p_value = s$stats[1, "Pr(>F)"]
  )
}

#' Group-level Pillai scores
#'
#' Computes Pillai scores and associated p-values per group (e.g., per speaker).
#'
#' @param data Data frame.
#' @param group_col String; grouping column (e.g., "speaker").
#' @param category_col String; category column (e.g., "vowel").
#' @param features Character vector of numeric feature columns.
#' @param min_tokens Minimum tokens per group.
#'
#' @return A tibble with columns: group, n_tokens, pillai, p_value.
#' @export
#' @importFrom dplyr group_by summarize n n_distinct filter ungroup rename
#' @importFrom rlang .data
speaker_pillai <- function(data,
                           group_col,
                           category_col,
                           features,
                           min_tokens = 20) {

  if (!group_col %in% names(data)) {
    stop("`group_col` must be in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("`category_col` must be in `data`.")
  }

  dplyr::group_by(data, .data[[group_col]]) |>
    dplyr::filter(
      dplyr::n() >= min_tokens,
      dplyr::n_distinct(.data[[category_col]]) >= 2L
    ) |>
    dplyr::summarize(
      n_tokens = dplyr::n(),
      pillai   = pillai_overlap(dplyr::cur_data_all(), features, category_col)$pillai,
      p_value  = pillai_overlap(dplyr::cur_data_all(), features, category_col)$p_value,
      .groups  = "drop"
    ) |>
    dplyr::rename(group = 1)
}

#' Global Pillai trace (point estimate)
#'
#' Computes a single Pillai-Bartlett trace for the full dataset in the
#' specified feature space, returning a one-row data frame that includes
#' the total number of tokens used.
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; column with category labels (>= 2 levels).
#' @param min_tokens Minimum total tokens required after removing missing
#'   values.
#'
#' @return A one-row data frame with columns:
#'   \code{n_tokens}, \code{pillai}, and \code{p_value}.
#' @export
global_pillai <- function(data,
                          features,
                          category_col,
                          min_tokens = 20) {

  if (!category_col %in% names(data)) {
    stop("`category_col` must be in `data`.")
  }
  if (!all(features %in% names(data))) {
    stop("All `features` must be in `data`.")
  }

  keep_cols <- c(category_col, features)

  df <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]),
             keep_cols, drop = FALSE]

  n <- nrow(df)
  if (n < min_tokens) {
    stop("Not enough tokens after removing missing values. Got ",
         n, ", need at least ", min_tokens, ".")
  }

  po <- pillai_overlap(
    data         = df,
    features     = features,
    category_col = category_col
  )

  data.frame(
    n_tokens = n,
    pillai   = po$pillai,
    p_value  = po$p_value,
    stringsAsFactors = FALSE
  )
}

#' Estimate Pillai trace, globally or by group
#'
#' Unified front-end for Pillai-Bartlett trace. If \code{group_col} is
#' \code{NULL}, computes a single global Pillai trace for the full dataset.
#' If \code{group_col} is provided, computes Pillai per group (e.g., per
#' speaker).
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; category column name.
#' @param group_col Optional string; grouping column name. If \code{NULL},
#'   a global Pillai value is returned.
#' @param min_tokens Minimum tokens (globally or per group).
#'
#' @return A data frame with either one global row or one row per group.
#' @export
estimate_pillai <- function(data,
                            features,
                            category_col,
                            group_col  = NULL,
                            min_tokens = 20) {

  if (is.null(group_col)) {
    # global
    gp <- global_pillai(
      data         = data,
      features     = features,
      category_col = category_col,
      min_tokens   = min_tokens
    )
    gp$scope <- "global"
    gp <- gp[, c("scope", setdiff(names(gp), "scope"))]
    return(gp)
  } else {
    # grouped
    sp <- speaker_pillai(
      data         = data,
      group_col    = group_col,
      category_col = category_col,
      features     = features,
      min_tokens   = min_tokens
    )
    sp$scope <- "group"
    sp <- sp[, c("scope", "group", "n_tokens", "pillai", "p_value")]
    return(sp)
  }
}




#' Bhattacharyya distance and affinity under multivariate normality
#'
#' Estimates means and covariances per category and computes
#' Bhattacharyya distance and affinity (exp(-distance)) under the
#' assumption of multivariate normality.
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; column with exactly two categories.
#' @param eps Small ridge constant added to covariance matrices to improve
#'   numerical stability.
#'
#' @return A list with `distance` and `affinity` (exp(-distance)).
#' @export
#' @importFrom stats cov
bhattacharyya_mvnorm <- function(data, features, category_col, eps = 1e-6) {
  if (!category_col %in% names(data)) {
    stop("`category_col` must be in `data`.")
  }

  levs <- unique(data[[category_col]])
  if (length(levs) != 2L) {
    stop("`category_col` must have exactly two levels.")
  }
  X1 <- as.matrix(data[data[[category_col]] == levs[1], features, drop = FALSE])
  X2 <- as.matrix(data[data[[category_col]] == levs[2], features, drop = FALSE])

  mu1 <- colMeans(X1); mu2 <- colMeans(X2)
  S1  <- stats::cov(X1); S2  <- stats::cov(X2)

  if (any(!is.finite(S1)) || any(!is.finite(S2))) {
    stop("Bhattacharyya distance failed: non-finite covariance estimates.")
  }

  # Ridge regularization to improve numerical stability
  S1 <- S1 + diag(eps, ncol(S1))
  S2 <- S2 + diag(eps, ncol(S2))
  S   <- (S1 + S2) / 2

  diff <- matrix(mu2 - mu1, ncol = 1)
  invS <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(invS)) {
    stop("Bhattacharyya distance failed: covariance not positive definite. ",
         "Try increasing `eps` or reducing feature dimensionality.")
  }

  detS  <- tryCatch(det(S),  error = function(e) NA_real_)
  detS1 <- tryCatch(det(S1), error = function(e) NA_real_)
  detS2 <- tryCatch(det(S2), error = function(e) NA_real_)
  if (!is.finite(detS) || !is.finite(detS1) || !is.finite(detS2) ||
      detS <= 0 || detS1 <= 0 || detS2 <= 0) {
    stop("Bhattacharyya distance failed: non-positive determinant. ",
         "Try increasing `eps` or reducing feature dimensionality.")
  }

  term1 <- 0.125 * t(diff) %*% invS %*% diff
  term2 <- 0.5 * log(detS / sqrt(detS1 * detS2))
  d <- as.numeric(term1 + term2)
  list(
    distance = d,
    affinity = exp(-d)
  )
}

#' Bhattacharyya distance by group
#'
#' Computes Bhattacharyya distance and affinity for each group, under
#' a multivariate normal approximation.
#'
#' @param data Data frame.
#' @param group_col String; grouping column (e.g., "speaker").
#' @param category_col String; category column with exactly two levels per group.
#' @param features Character vector of numeric feature columns.
#' @param min_tokens Minimum tokens per group.
#' @param eps Small ridge constant passed to \code{bhattacharyya_mvnorm()}.
#'
#' @return Data frame with columns: group, n_tokens, bhatt_dist, bhatt_affinity.
#' @export
speaker_bhatt <- function(data,
                          group_col,
                          category_col,
                          features,
                          min_tokens = 20,
                          eps = 1e-6) {

  if (!group_col %in% names(data)) {
    stop("`group_col` must be in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("`category_col` must be in `data`.")
  }
  if (!all(features %in% names(data))) {
    stop("All `features` must be in `data`.")
  }

  groups <- split(data, data[[group_col]])
  out_list <- lapply(groups, function(df_g) {
    n_tok <- nrow(df_g)
    levs <- unique(df_g[[category_col]])
    if (n_tok < min_tokens || length(levs) != 2L) {
      return(NULL)
    }
    bh <- bhattacharyya_mvnorm(
      data         = df_g,
      features     = features,
      category_col = category_col,
      eps          = eps
    )
    data.frame(
      group         = df_g[[group_col]][1],
      n_tokens      = n_tok,
      bhatt_dist    = bh$distance,
      bhatt_affinity = bh$affinity,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out_list)
  if (is.null(out)) {
    out <- data.frame(
      group          = character(0),
      n_tokens       = integer(0),
      bhatt_dist     = numeric(0),
      bhatt_affinity = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}

#' Estimate Bhattacharyya distance, globally or by group
#'
#' Unified front-end for Bhattacharyya distance (and affinity) under a
#' multivariate normal approximation.
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; category column name (exactly two levels globally).
#' @param group_col Optional string; grouping column name. If \code{NULL},
#'   a single global Bhattacharyya distance is returned.
#' @param min_tokens Minimum tokens (globally or per group).
#' @param eps Small ridge constant passed to \code{bhattacharyya_mvnorm()}.
#'
#' @return Data frame with either one global row or one row per group.
#' @export
estimate_bhatt <- function(data,
                           features,
                           category_col,
                           group_col  = NULL,
                           min_tokens = 20,
                           eps = 1e-6) {

  if (is.null(group_col)) {
    # global
    keep_cols <- c(category_col, features)
    df <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]),
               keep_cols, drop = FALSE]

    n <- nrow(df)
    if (n < min_tokens) {
      stop("Not enough tokens after removing missing values. Got ",
           n, ", need at least ", min_tokens, ".")
    }

    bh <- bhattacharyya_mvnorm(
      data         = df,
      features     = features,
      category_col = category_col,
      eps          = eps
    )

    out <- data.frame(
      scope          = "global",
      n_tokens       = n,
      bhatt_dist     = bh$distance,
      bhatt_affinity = bh$affinity,
      stringsAsFactors = FALSE
    )
    return(out)
  } else {
    # grouped
    sb <- speaker_bhatt(
      data         = data,
      group_col    = group_col,
      category_col = category_col,
      features     = features,
      min_tokens   = min_tokens,
      eps          = eps
    )
    sb$scope <- "group"
    sb <- sb[, c("scope", "group", "n_tokens", "bhatt_dist", "bhatt_affinity")]
    return(sb)
  }
}
