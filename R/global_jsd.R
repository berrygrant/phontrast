#' Global JSD with bootstrap confidence interval
#'
#' Computes a single Jensen–Shannon divergence (JSD) value for two categories
#' in an n-dimensional acoustic space, together with bootstrap-based
#' confidence intervals obtained by resampling tokens with replacement.
#'
#' This is the "group-wise" version of JSD: it ignores speakers and treats
#' all tokens as coming from a single population for each category.
#'
#' @param data Data frame containing at least the category column and
#'   the feature columns.
#' @param features Character vector of column names giving the acoustic
#'   dimensions (e.g., c("f1", "f2") or paste0("mfcc", 1:13)).
#' @param category_col String; name of the column giving the two categories
#'   to compare (e.g., "vowel"). Must have exactly two unique values.
#' @param n_boot Integer; number of bootstrap resamples.
#' @param min_tokens Minimum total number of non-missing tokens required.
#' @param est_distance Logical; if TRUE, return Jensen-Shannon distance
#'   (sqrt of divergence) instead of divergence.
#' @param conf_level Confidence level for bootstrap intervals.
#' @param bw Bandwidth selection method passed to \code{jsd_kde_nd()}.
#' @param eval_on KDE evaluation points passed to \code{jsd_kde_nd()}.
#' @param eval_n Optional maximum number of KDE evaluation points.
#' @param eval_seed Optional integer seed for KDE evaluation-point subsampling.
#' @param engine KDE evaluation engine passed to \code{jsd_kde_nd()}.
#'   \code{"fast_diagonal"} is accepted as an alias for \code{"fast_diag"}.
#' @param chunk_size Chunk size for \code{engine = "fast_diag"}.
#' @param method Estimator passed to \code{jsd_kde_nd()}: \code{"mc"} (default)
#'   or \code{"legacy"} (pre-1.1.0 self-normalized estimate).
#' @param ... Additional arguments passed to \code{jsd_kde_nd()}.
#'
#' @return A one-row data frame with columns:
#'   \itemize{
#'     \item \code{n_tokens} – total number of tokens used
#'     \item \code{n_boot} – number of successful bootstrap samples
#'     \item \code{conf_level} – confidence level used for the interval
#'     \item \code{jsd_point} – JSD on the full dataset
#'     \item \code{jsd_mean} – mean JSD across bootstrap samples
#'     \item \code{jsd_sd} – standard deviation of bootstrap JSD
#'     \item \code{ci_lower}, \code{ci_upper} – bootstrap confidence interval
#'     \item \code{jsd_low}, \code{jsd_high} – legacy aliases for
#'       \code{ci_lower} and \code{ci_upper}
#'   }
#' @export
#' @importFrom stats sd quantile
global_boot_jsd <- function(data,
                            features,
                            category_col,
                            n_boot     = 1000,
                            min_tokens = 20,
                            est_distance = FALSE,
                            conf_level = 0.95,
                            bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                            eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                            eval_n = NULL,
                            eval_seed = NULL,
                            engine = c("ks", "fast_diag", "fast_diagonal"),
                            chunk_size = 1000L,
                            method = c("mc", "legacy"),
                            ...) {

  out <- estimate_jsd(
    data = data,
    features = features,
    category_col = category_col,
    group_col = NULL,
    do_boot = TRUE,
    n_boot = n_boot,
    min_tokens = min_tokens,
    est_distance = est_distance,
    conf_level = conf_level,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size,
    method = method,
    ...
  )
  out$scope <- NULL
  tibble::as_tibble(out)
}
