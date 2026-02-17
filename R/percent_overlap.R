#' Percent overlap between two distributions via KDE
#'
#' Computes the percentage overlap (shared area) between two categories
#' in an n-dimensional acoustic space using multivariate kernel density
#' estimation. Range: 0 = no overlap, 1 = identical.
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; exactly two categories.
#' @param ... Additional arguments passed to KDE bandwidth selection.
#'
#' @return Numeric scalar in \code{[0, 1]}.
#' @export
percent_overlap_kde <- function(data,
                                features,
                                category_col,
                                ...) {

  if (!category_col %in% names(data))
    stop("`category_col` must be a column in `data`.")

  if (!all(features %in% names(data)))
    stop("Not all `features` exist in `data`.")

  levs <- unique(data[[category_col]])
  if (length(levs) != 2L)
    stop("`category_col` must have exactly two categories.")

  d1 <- data[data[[category_col]] == levs[1], features, drop = FALSE]
  d2 <- data[data[[category_col]] == levs[2], features, drop = FALSE]

  X1 <- as.matrix(d1)
  X2 <- as.matrix(d2)
  X_all <- rbind(X1, X2)

  H1 <- ks::Hpi(X1)
  H2 <- ks::Hpi(X2)

  kde1 <- ks::kde(x = X1, H = H1, eval.points = X_all)
  kde2 <- ks::kde(x = X2, H = H2, eval.points = X_all)

  p <- as.numeric(kde1$estimate)
  q <- as.numeric(kde2$estimate)

  # Normalize to discrete probability masses on a shared grid
  p <- p / sum(p)
  q <- q / sum(q)

  # Percent overlap is the shared area on the grid: sum(min(p, q))
  overlap <- sum(pmin(p, q))

  # Bound numerically
  overlap <- max(min(overlap, 1), 0)

  overlap
}

#' Estimate percent overlap globally or by group
#'
#' Unified front-end for KDE-based percent overlap between two categories.
#'
#' @inheritParams estimate_jsd
#'
#' @return A data frame (global = one row; grouped = one per group).
#' @export
estimate_overlap <- function(data,
                             features,
                             category_col,
                             group_col  = NULL,
                             min_tokens = 20) {

  if (is.null(group_col)) {
    # ---- Global ----
    keep_cols <- c(category_col, features)
    df <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]),
               keep_cols, drop = FALSE]

    n <- nrow(df)
    if (n < min_tokens)
      stop("Not enough tokens for global percent overlap.")

    ov <- percent_overlap_kde(
      data         = df,
      features     = features,
      category_col = category_col
    )

    return(data.frame(
      scope        = "global",
      n_tokens     = n,
      overlap      = ov,
      stringsAsFactors = FALSE
    ))
  }

  # ---- Grouped ----
  if (!group_col %in% names(data))
    stop("`group_col` must be a column in `data`.")

  groups <- split(data, data[[group_col]])

  out <- lapply(groups, function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens ||
        length(unique(df_g[[category_col]])) != 2L)
      return(NULL)

    ov <- percent_overlap_kde(
      data         = df_g,
      features     = features,
      category_col = category_col
    )

    data.frame(
      scope    = "group",
      group    = df_g[[group_col]][1],
      n_tokens = n_tok,
      overlap  = ov,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
