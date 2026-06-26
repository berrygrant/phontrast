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
#' @param ... Additional arguments passed to \code{jsd_kde_nd()}.
#'
#' @return A one-row data frame with columns:
#'   \itemize{
#'     \item \code{n_tokens} – total number of tokens used
#'     \item \code{n_boot} – number of successful bootstrap samples
#'     \item \code{jsd_point} – JSD on the full dataset
#'     \item \code{jsd_mean} – mean JSD across bootstrap samples
#'     \item \code{jsd_sd} – standard deviation of bootstrap JSD
#'     \item \code{jsd_low}, \code{jsd_high} – 95% bootstrap CI (2.5%, 97.5%)
#'   }
#' @export
#' @importFrom stats sd quantile
global_boot_jsd <- function(data,
                            features,
                            category_col,
                            n_boot     = 300,
                            min_tokens = 20,
                            est_distance = FALSE,
                            conf_level = 0.95,
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
    ...
  )
  out$scope <- NULL
  as.data.frame(out)
}
