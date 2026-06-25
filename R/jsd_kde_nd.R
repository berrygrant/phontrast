#' n-dimensional JSD via multivariate kernel density estimation
#'
#' Computes Jensen–Shannon divergence between two categories in an
#' arbitrary n-dimensional acoustic space using multivariate KDE from
#' the \pkg{ks} package.
#'
#' @param data A data frame containing observations from exactly two categories.
#' @param features Character vector of column names giving the acoustic
#'   dimensions (e.g., MFCC1..MFCC13, F1/F2/duration).
#' @param group String: name of the column giving the category labels
#'   (e.g., "vowel", "segment"). Must have exactly two unique values in `data`.
#' @param bw Bandwidth selection method. One of \code{"Hpi"},
#'   \code{"Hscv"}, or \code{"Hpi.diag"}. Passed to \code{ks::Hpi()},
#'   \code{ks::Hscv()}, or \code{ks::Hpi.diag()}.
#' @param eval_on Where to evaluate the KDEs. "pooled" (default) evaluates
#'   on all observations from both categories; "group1" or "group2"
#'   evaluate on the respective group only.
#'
#' @return A single numeric JSD value in bits.
#' @export
#' @importFrom ks Hpi Hscv Hpi.diag kde
#' @importFrom rlang .data
jsd_kde_nd <- function(data,
                       features,
                       group   = "category",
                       bw      = c("Hpi", "Hscv", "Hpi.diag"),
                       eval_on = c("pooled", "group1", "group2")) {

  bw      <- match.arg(bw)
  eval_on <- match.arg(eval_on)

  .check_columns(data, c(group, features))
  data <- .metric_data(data, c(group, features))

  levs <- .two_levels(data[[group]], "group")

  d1 <- data[data[[group]] == levs[1], , drop = FALSE]
  d2 <- data[data[[group]] == levs[2], , drop = FALSE]

  X1    <- as.matrix(d1[, features, drop = FALSE])
  X2    <- as.matrix(d2[, features, drop = FALSE])
  X_all <- as.matrix(data[, features, drop = FALSE])

  eval_pts <- switch(
    eval_on,
    pooled = X_all,
    group1 = X1,
    group2 = X2
  )

  # Bandwidths for each group
  H1 <- switch(
    bw,
    Hpi      = ks::Hpi(X1),
    Hscv     = ks::Hscv(X1),
    Hpi.diag = ks::Hpi.diag(X1)
  )
  H2 <- switch(
    bw,
    Hpi      = ks::Hpi(X2),
    Hscv     = ks::Hscv(X2),
    Hpi.diag = ks::Hpi.diag(X2)
  )

  kde1 <- ks::kde(x = X1, H = H1, eval.points = eval_pts)
  kde2 <- ks::kde(x = X2, H = H2, eval.points = eval_pts)

  p <- as.numeric(kde1$estimate)
  q <- as.numeric(kde2$estimate)

  if (any(!is.finite(p)) || any(!is.finite(q)) || sum(p) <= 0 || sum(q) <= 0) {
    stop("KDE returned invalid density estimates.", call. = FALSE)
  }

  jsd(p, q)
}
