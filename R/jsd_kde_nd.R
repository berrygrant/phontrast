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
#'   \code{ks::Hscv()}, or \code{ks::Hpi.diag()} for multivariate inputs.
#'   For one-dimensional inputs, these map to \code{stats::bw.SJ()},
#'   \code{stats::bw.ucv()}, and \code{stats::bw.nrd0()}, respectively,
#'   with a robust fallback for constant samples.
#' @param eval_on Where to evaluate the KDEs. "pooled" (default) evaluates
#'   on all observations from both categories; "group1" or "group2"
#'   evaluate on the respective group only.
#'
#' @return A single numeric JSD value in bits.
#'
#' @examples
#' set.seed(2026)
#' vowels <- data.frame(
#'   vowel = rep(c("ih", "eh"), each = 40),
#'   f1 = c(rnorm(40, 500, 55), rnorm(40, 565, 60)),
#'   f2 = c(rnorm(40, 1980, 150), rnorm(40, 1870, 155))
#' )
#'
#' # One-dimensional JSD, for example a single formant or duration.
#' jsd_kde_nd(vowels, features = "f1", group = "vowel")
#'
#' # Two-dimensional JSD in F1/F2 space.
#' jsd_kde_nd(vowels, features = c("f1", "f2"), group = "vowel")
#' @export
#' @importFrom ks Hpi Hscv Hpi.diag kde
#' @importFrom rlang .data
jsd_kde_nd <- function(data,
                       features,
                       group   = "category",
                       bw      = c("Hpi", "Hscv", "Hpi.diag"),
                       eval_on = c("pooled", "group1", "group2")) {

  dens <- .kde_density_pair(
    data = data,
    features = features,
    category_col = group,
    bw = bw,
    eval_on = eval_on,
    metric = "jsd_kde_nd()"
  )

  jsd(dens$p, dens$q)
}
