#' Percent overlap between two distributions via KDE
#'
#' Computes the percentage overlap (shared area) between two categories
#' in an n-dimensional acoustic space using multivariate kernel density
#' estimation. Range: 0 = no overlap, 1 = identical.
#'
#' @param data Data frame.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; exactly two categories.
#' @param bw Bandwidth selection method. Uses the same options as
#'   \code{jsd_kde_nd()}: \code{"Hpi"}, \code{"Hscv"}, or \code{"Hpi.diag"}.
#' @param eval_on Where to evaluate the KDEs. Uses the same options as
#'   \code{jsd_kde_nd()}: \code{"pooled"}, \code{"group1"}, or \code{"group2"}.
#' @param ... Reserved for future extensions; currently unused.
#'
#' @return Numeric scalar in \code{[0, 1]}.
#' @export
percent_overlap_kde <- function(data,
                                features,
                                category_col,
                                bw = c("Hpi", "Hscv", "Hpi.diag"),
                                eval_on = c("pooled", "group1", "group2"),
                                ...) {

  dens <- .kde_density_pair(
    data = data,
    features = features,
    category_col = category_col,
    bw = bw,
    eval_on = eval_on,
    metric = "percent_overlap_kde()"
  )

  # Normalize to discrete probability masses on a shared grid
  p <- dens$p / sum(dens$p)
  q <- dens$q / sum(dens$q)

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
                             min_tokens = 20,
                             bw = c("Hpi", "Hscv", "Hpi.diag"),
                             eval_on = c("pooled", "group1", "group2"),
                             ...) {

  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  .check_positive_count(min_tokens, "min_tokens")

  if (is.null(group_col)) {
    # ---- Global ----
    keep_cols <- c(category_col, features)
    df <- .metric_data(data, keep_cols)

    n <- nrow(df)
    if (n < min_tokens)
      stop("Not enough tokens for global percent overlap.")

    ov <- percent_overlap_kde(
      data         = df,
      features     = features,
      category_col = category_col,
      bw           = bw,
      eval_on      = eval_on,
      ...
    )

    return(data.frame(
      scope        = "global",
      n_tokens     = n,
      overlap      = ov,
      stringsAsFactors = FALSE
    ))
  }

  # ---- Grouped ----
  .check_columns(data, c(group_col, category_col, features))
  data <- .metric_data(data, c(group_col, category_col, features))

  groups <- .split_groups(data, group_col)

  out <- lapply(groups, function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens ||
        .observed_n_categories(df_g[[category_col]]) != 2L)
      return(NULL)

    ov <- tryCatch(
      percent_overlap_kde(
        data         = df_g,
        features     = features,
        category_col = category_col,
        bw           = bw,
        eval_on      = eval_on,
        ...
      ),
      error = function(e) NA_real_
    )

    data.frame(
      scope    = "group",
      group    = df_g[[group_col]][1],
      n_tokens = n_tok,
      overlap  = ov,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(
      scope = character(),
      group = character(),
      n_tokens = integer(),
      overlap = numeric(),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}
