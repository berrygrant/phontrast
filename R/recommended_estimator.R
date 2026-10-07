#' Kernel estimator settings by dimensionality
#'
#' Returns the kernel-density estimator configuration that the simulation
#' study behind \code{rank_contrasts()} used at each dimensionality (Berry,
#' under review, Table III), so that a contrast measured with phontrast is
#' scored with the settings the safe-use envelope was calibrated on.
#'
#' The study computed every \eqn{\sqrt{JSD}} with phontrast's
#' \code{method = "legacy"} estimator: each category's kernel density is
#' evaluated at the pooled tokens of both categories, the two density vectors
#' are self-normalized over those points, and the divergence is read off them.
#' The protocol's floors, its 0.25 margin, its 0.99 ceiling, and its bandwidth
#' check were calibrated on that estimator, so \code{method} is
#' \code{"legacy"} at every dimensionality, and no leave-one-out correction is
#' applied (\code{legacy} has none). The bandwidth check of the study
#' (\code{bracket_bw = "scott.pooled"}) used one diagonal Scott bandwidth
#' selected on the pooled tokens of the pair for both categories. The
#' bandwidth rule, engine, and evaluation subsample switch with the number of
#' features \code{d}:
#'
#' \describe{
#'   \item{\code{d <= 4}}{Plug-in bandwidth (\code{bw = "Hpi"}), the
#'     \code{"ks"} engine, densities evaluated at every pooled token.
#'     Calibrated at \code{d = 2} and \code{4}.}
#'   \item{\code{5 <= d <= 13}}{Diagonal Scott bandwidth
#'     (\code{bw = "scott.diag"}), the \code{"fast_diag"} engine, densities
#'     evaluated at 200 subsampled pooled tokens. Calibrated at \code{d = 8}
#'     and \code{13}; dimensionalities between the calibrated ones take the
#'     settings of the next higher calibrated dimensionality.}
#'   \item{\code{d >= 14}}{The same settings, reproducing the log-space path
#'     the study used at \code{d = 32} and \code{64}. phontrast's
#'     \code{"fast_diag"} engine evaluates kernels in log space, so this tier
#'     runs natively. Above eight dimensions the study gives ordering evidence
#'     only.}
#' }
#'
#' @param d Positive integer; the number of acoustic features (dimensions).
#'
#' @return A list with elements \code{d}, \code{tier} (a label for the row of
#'   Table III applied), \code{calibrated_at} (the dimensionalities the row was
#'   calibrated on), \code{method}, \code{bw}, \code{engine},
#'   \code{eval_on}, \code{eval_n}, \code{loo}, \code{bracket_bw}, and
#'   \code{note}. The \code{method}, \code{bw}, \code{engine},
#'   \code{eval_on}, \code{eval_n}, and \code{loo} elements can be passed
#'   straight to \code{jsd_kde_nd()}, \code{estimate_jsd()}, or
#'   \code{phontrast()}. \code{bracket_bw} is read by \code{rank_contrasts()}
#'   only: \code{"scott.pooled"} for the study's pooled diagonal Scott
#'   bandwidth, or \code{"scott.diag"} for each category's own.
#'
#' @examples
#' recommended_estimator(2)
#' recommended_estimator(13)$engine
#'
#' set.seed(2026)
#' vowels <- data.frame(
#'   vowel = rep(c("ih", "eh"), each = 40),
#'   f1 = c(rnorm(40, 500, 55), rnorm(40, 565, 60)),
#'   f2 = c(rnorm(40, 1980, 150), rnorm(40, 1870, 155))
#' )
#' est <- recommended_estimator(2)
#' jsd_kde_nd(vowels, c("f1", "f2"), "vowel", method = est$method,
#'            bw = est$bw, engine = est$engine, eval_on = est$eval_on)
#' @export
recommended_estimator <- function(d) {
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 1 || d != round(d)) {
    stop("`d` must be a single positive integer (the number of features).", call. = FALSE)
  }
  d <- as.integer(d)

  if (d <= 4L) {
    out <- list(
      tier = "d <= 4",
      calibrated_at = c(2L, 4L),
      method = "legacy",
      bw = "Hpi",
      engine = "ks",
      eval_on = "pooled",
      eval_n = NULL,
      loo = FALSE,
      bracket_bw = "scott.pooled",
      note = paste(
        "Legacy (self-normalized) estimator, plug-in (Hpi) bandwidth, ks engine,",
        "densities evaluated at all pooled tokens, no leave-one-out correction."
      )
    )
  } else if (d <= 13L) {
    out <- list(
      tier = "5 <= d <= 13",
      calibrated_at = c(8L, 13L),
      method = "legacy",
      bw = "scott.diag",
      engine = "fast_diag",
      eval_on = "pooled",
      eval_n = 200L,
      loo = FALSE,
      bracket_bw = "scott.pooled",
      note = paste(
        "Legacy (self-normalized) estimator, diagonal Scott bandwidth, fast_diag",
        "engine, densities evaluated at 200 subsampled pooled tokens, no",
        "leave-one-out correction."
      )
    )
  } else {
    out <- list(
      tier = "d >= 14",
      calibrated_at = c(32L, 64L),
      method = "legacy",
      bw = "scott.diag",
      engine = "fast_diag",
      eval_on = "pooled",
      eval_n = 200L,
      loo = FALSE,
      bracket_bw = "scott.pooled",
      note = paste(
        "Legacy (self-normalized) estimator, diagonal Scott bandwidth, log-space",
        "fast_diag engine, 200 subsampled pooled tokens, no leave-one-out",
        "correction (the simulation code's path at d = 32 and 64). Above eight",
        "dimensions the study gives ordering evidence only."
      )
    )
  }
  c(list(d = d), out)
}
