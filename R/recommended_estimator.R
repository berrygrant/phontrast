#' Kernel estimator settings by dimensionality
#'
#' Returns the kernel-density estimator configuration that the simulation
#' study behind \code{rank_contrasts()} used at each dimensionality (Berry,
#' under review, Table III), so that a contrast measured with phontrast is
#' scored with the settings the safe-use envelope was calibrated on. Three
#' settings switch together with the number of features \code{d}:
#'
#' \describe{
#'   \item{\code{d <= 4}}{Plug-in bandwidth (\code{bw = "Hpi"}), the
#'     \code{"ks"} engine, densities evaluated at every token, and the partial
#'     leave-one-out correction on. Calibrated at \code{d = 2} and \code{4}.}
#'   \item{\code{5 <= d <= 13}}{Diagonal Scott bandwidth
#'     (\code{bw = "scott.diag"}), the \code{"fast_diag"} engine, densities
#'     evaluated at 200 subsampled tokens per category, leave-one-out on.
#'     Calibrated at \code{d = 8} and \code{13}; dimensionalities between the
#'     calibrated ones take the settings of the next higher calibrated
#'     dimensionality.}
#'   \item{\code{d >= 14}}{As above but with no leave-one-out correction,
#'     reproducing the log-space path the study used at \code{d = 32} and
#'     \code{64}. phontrast's \code{"fast_diag"} engine evaluates kernels in
#'     log space, so this tier runs natively. Above eight dimensions the study
#'     gives ordering evidence only.}
#' }
#'
#' @param d Positive integer; the number of acoustic features (dimensions).
#'
#' @return A list with elements \code{d}, \code{tier} (a label for the row of
#'   Table III applied), \code{calibrated_at} (the dimensionalities the row was
#'   calibrated on), \code{bw}, \code{engine}, \code{eval_n}, \code{loo}, and
#'   \code{note}. The \code{bw}, \code{engine}, \code{eval_n}, and \code{loo}
#'   elements can be passed straight to \code{jsd_kde_nd()},
#'   \code{estimate_jsd()}, or \code{phontrast()}.
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
#' jsd_kde_nd(vowels, c("f1", "f2"), "vowel",
#'            bw = est$bw, engine = est$engine, loo = est$loo)
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
      bw = "Hpi",
      engine = "ks",
      eval_n = NULL,
      loo = TRUE,
      note = paste(
        "Plug-in (Hpi) bandwidth, ks engine, densities evaluated at all tokens,",
        "partial leave-one-out correction on."
      )
    )
  } else if (d <= 13L) {
    out <- list(
      tier = "5 <= d <= 13",
      calibrated_at = c(8L, 13L),
      bw = "scott.diag",
      engine = "fast_diag",
      eval_n = 200L,
      loo = TRUE,
      note = paste(
        "Diagonal Scott bandwidth, fast_diag engine, densities evaluated at 200",
        "subsampled tokens per category, partial leave-one-out correction on."
      )
    )
  } else {
    out <- list(
      tier = "d >= 14",
      calibrated_at = c(32L, 64L),
      bw = "scott.diag",
      engine = "fast_diag",
      eval_n = 200L,
      loo = FALSE,
      note = paste(
        "Diagonal Scott bandwidth, log-space fast_diag engine, 200 subsampled",
        "tokens per category, no leave-one-out correction (the simulation code's",
        "path at d = 32 and 64). Above eight dimensions the study gives ordering",
        "evidence only."
      )
    )
  }
  c(list(d = d), out)
}
