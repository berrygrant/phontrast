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
#' estimator of the continuous JSD in any dimension.
#'
#' `method = "legacy"` is the self-normalized sample-point estimator: both
#' categories' KDEs are evaluated at one shared set of points (by default the
#' pooled tokens of both categories; see `eval_on`), each density vector is
#' normalized to sum to one over those points, and the discrete JSD of the two
#' vectors is returned. It is a bounded relative separation index rather than a
#' consistent estimate of the continuous JSD, and on overlapping categories it
#' reads lower than `"mc"`. It was phontrast's estimator before 1.2.0, and it
#' is the estimator on which the ranking protocol of `rank_contrasts()` was
#' validated (Berry, under review, Sec. VII A): the protocol's sample-size
#' floors, margin, ceiling, and bandwidth check were calibrated on it, so
#' `recommended_estimator()` returns it and `rank_contrasts()` uses it. Use it
#' to report values against that protocol or to reproduce the study's
#' numbers. No leave-one-out correction is applied under `"legacy"`.
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
#'   estimate of the continuous JSD, or \code{"legacy"} for the self-normalized
#'   sample-point index on which the \code{rank_contrasts()} protocol was
#'   calibrated (see Details). Ignored when \code{density = "mvnorm"}.
#' @param density Density model behind the estimate: \code{"kde"} (default)
#'   estimates each category's density by kernel density estimation;
#'   \code{"mvnorm"} fits one multivariate normal per category and estimates the
#'   continuous JSD between the two Gaussians by Monte-Carlo (no closed form
#'   exists). Under \code{"mvnorm"} the KDE-specific arguments (\code{bw},
#'   \code{engine}, \code{eval_on}, \code{chunk_size}, \code{method},
#'   \code{eval_n}, \code{loo}) do not apply; the Monte-Carlo sample size is set
#'   by \code{mc_n} and \code{eval_seed} makes the draw reproducible.
#' @param mc_n Positive integer; number of Monte-Carlo samples drawn from each
#'   fitted Gaussian when \code{density = "mvnorm"} (default \code{10000}). The
#'   estimator draws \code{mc_n} fresh points from each category's fitted
#'   Gaussian and averages the log density ratio, so it targets the JSD between
#'   the two fitted Gaussians rather than a resubstitution estimate at the
#'   observed points. Larger values reduce Monte-Carlo variance. Ignored when
#'   \code{density = "kde"}.
#' @param loo Logical; if \code{TRUE} (default) the Monte-Carlo estimator uses a
#'   partial leave-one-out correction on each category's self-density to reduce
#'   resubstitution bias. The correction removes a sample-size-scaled fraction
#'   \code{n / (n + 20)} of each point's own kernel: half at 20 tokens per
#'   category (the \code{min_tokens} default), approaching the full
#'   leave-one-out correction as the category grows. Removing only part of the
#'   self-kernel keeps the corrected density strictly positive at isolated
#'   points, so small but real divergences remain small positive values rather
#'   than being floored to exactly 0 (as the full leave-one-out correction did
#'   through phontrast 2.0.2). Ignored when \code{method = "legacy"}, which
#'   applies no leave-one-out correction.
#' @param bw_scale Positive number multiplying the selected kernel bandwidth on
#'   the standard-deviation scale: univariate bandwidths are multiplied by
#'   \code{bw_scale} and bandwidth matrices by \code{bw_scale^2}. The default
#'   \code{1} uses the selected bandwidth unchanged; \code{0.5} and \code{2}
#'   give the halved and doubled bandwidths of the smoothing-sensitivity check
#'   in \code{rank_contrasts()}. Ignored when \code{density = "mvnorm"}.
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
                       density = c("kde", "mvnorm"),
                       mc_n = 10000L,
                       loo = TRUE,
                       bw_scale = 1) {

  .validate_metric_inputs(data, features, group)
  method <- match.arg(method)
  density <- match.arg(density)
  .check_bool(loo, "loo")
  .check_bw_scale(bw_scale)

  if (identical(density, "mvnorm")) {
    mc <- .mvnorm_mc_pair(
      data = data,
      features = features,
      category_col = group,
      mc_n = mc_n,
      eval_seed = eval_seed,
      metric = "jsd_kde_nd()"
    )
    return(.jsd_mc(mc, loo = FALSE))
  }

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
      metric = "jsd_kde_nd()",
      bw_scale = bw_scale
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
    metric = "jsd_kde_nd()",
    bw_scale = bw_scale
  )

  jsd(dens$p, dens$q)
}
