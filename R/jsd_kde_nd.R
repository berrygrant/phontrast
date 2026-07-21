#' n-dimensional JSD via multivariate kernel density estimation
#'
#' Computes Jensen-Shannon divergence between two categories in an
#' arbitrary n-dimensional acoustic space using multivariate KDE. The default
#' engine uses the \pkg{ks} package; a faster diagonal-Gaussian engine is
#' available for diagonal bandwidths.
#'
#' By default (`method = "mc"`) JSD is estimated with a Monte-Carlo plug-in:
#' each category's KDE is evaluated at that category's own observations and the
#' log density ratio against the mixture is averaged. This is a consistent
#' estimator of the continuous JSD in any dimension. `method = "legacy"`
#' reproduces the pre-1.2.0 self-normalized sample-point estimate (a bounded
#' relative separation index rather than the continuous JSD); use it only to
#' reproduce results from phonJSD 1.0.0.
#'
#' @param data A data frame containing observations from exactly two categories.
#' @param features Character vector of column names giving the acoustic
#'   dimensions (e.g., MFCC1..MFCC13, F1/F2/duration).
#' @param group String: name of the column giving the category labels
#'   (e.g., "vowel", "segment"). Must have exactly two unique values in `data`.
#' @param bw Bandwidth selection method. One of \code{"Hpi"},
#'   \code{"Hscv"}, \code{"Hpi.diag"}, or \code{"scott.diag"}. The first
#'   three are passed to \code{ks::Hpi()}, \code{ks::Hscv()}, or
#'   \code{ks::Hpi.diag()} for multivariate inputs. \code{"scott.diag"}
#'   uses a diagonal Scott rule-of-thumb bandwidth matrix.
#'   For one-dimensional inputs, these map to \code{stats::bw.SJ()},
#'   \code{stats::bw.ucv()}, \code{stats::bw.nrd0()}, and Scott's rule,
#'   respectively, with a robust fallback for constant samples.
#' @param eval_on Where to evaluate the KDEs (\code{method = "legacy"} only).
#'   "pooled" (default) evaluates on all observations from both categories;
#'   "group1" or "group2" evaluate on the respective group only.
#'   \code{"pooled_sample"} evaluates on a sampled subset of pooled observations
#'   and requires \code{eval_n}. Ignored when \code{method = "mc"} (which always
#'   evaluates each category at its own observations).
#' @param eval_n Optional positive integer giving the maximum number of
#'   evaluation points to use. If supplied, evaluation points are sampled from
#'   the set chosen by \code{eval_on}.
#' @param eval_seed Optional integer seed used only when \code{eval_n} causes
#'   evaluation-point subsampling. If \code{NULL}, the current R random-number
#'   state is used.
#' @param engine KDE evaluation engine. \code{"ks"} uses \code{ks::kde()}.
#'   \code{"fast_diag"} uses a chunked diagonal-Gaussian evaluator and requires
#'   \code{bw = "scott.diag"} or \code{bw = "Hpi.diag"} for multivariate KDE.
#'   \code{"fast_diagonal"} is accepted as an alias for \code{"fast_diag"}.
#' @param chunk_size Positive integer controlling the number of evaluation
#'   points processed per chunk by \code{engine = "fast_diag"}.
#' @param method Estimator: \code{"mc"} (default) for the Monte-Carlo plug-in
#'   estimate of the continuous JSD, or \code{"legacy"} for the pre-1.2.0
#'   self-normalized sample-point index.
#' @param loo Logical; if \code{TRUE} (default) the Monte-Carlo estimator uses a
#'   leave-one-out correction on each category's self-density to reduce
#'   resubstitution bias. Ignored when \code{method = "legacy"}.
#'
#' @return A single numeric JSD value in bits, bounded in \code{[0, 1]}.
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
#'
#' # Faster high-dimensional path: diagonal Scott bandwidth and sampled
#' # pooled evaluation points.
#' jsd_kde_nd(
#'   vowels,
#'   features = c("f1", "f2"),
#'   group = "vowel",
#'   bw = "scott.diag",
#'   eval_n = 40,
#'   eval_seed = 2026,
#'   engine = "fast_diag"
#' )
#' @export
#' @importFrom ks Hpi Hscv Hpi.diag kde
#' @importFrom rlang .data
jsd_kde_nd <- function(data,
                       features,
                       group   = "category",
                       bw      = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                       eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                       eval_n = NULL,
                       eval_seed = NULL,
                       engine = c("ks", "fast_diag", "fast_diagonal"),
                       chunk_size = 1000L,
                       method = c("mc", "legacy"),
                       loo = TRUE) {

  .validate_metric_inputs(data, features, group)
  method <- match.arg(method)
  .check_bool(loo, "loo")

  if (identical(method, "mc")) {
    mc <- .kde_mc_pair(
      data = data,
      features = features,
      category_col = group,
      bw = bw,
      eval_n = eval_n,
      eval_seed = eval_seed,
      engine = engine,
      chunk_size = chunk_size,
      metric = "jsd_kde_nd()"
    )
    return(.jsd_mc(mc, loo = loo))
  }

  dens <- .kde_density_pair(
    data = data,
    features = features,
    category_col = group,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size,
    metric = "jsd_kde_nd()"
  )

  jsd(dens$p, dens$q)
}
