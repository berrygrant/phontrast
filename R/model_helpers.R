#' Prepare JSD estimates for beta regression / GAMs
#'
#' JSD lives in \eqn{[0,1]}. This helper adds a `jsd_beta` column bounded in
#' \eqn{(0,1)} so it can be used with the Beta family (e.g., in \pkg{mgcv}).
#'
#' @param jsd_df Data frame containing a JSD column.
#' @param jsd_col String: name of the JSD column (default "jsd_mean").
#' @param eps Small constant used to bound JSD away from 0 and 1.
#'
#' @return A modified data frame with an added `jsd_beta` column.
#'
#' @examples
#' jsd_by_speaker <- data.frame(
#'   speaker = paste0("s", 1:8),
#'   age = c(18, 22, 27, 31, 38, 45, 52, 60),
#'   jsd_mean = c(0.02, 0.05, 0.08, 0.13, 0.18, 0.24, 0.31, 0.39)
#' )
#'
#' model_data <- prepare_jsd_beta(jsd_by_speaker)
#' model_data
#'
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   fit <- mgcv::gam(
#'     jsd_beta ~ age,
#'     data = model_data,
#'     family = mgcv::betar(),
#'     method = "REML"
#'   )
#'   stats::predict(fit, type = "response")
#' }
#' @export
prepare_jsd_beta <- function(jsd_df,
                             jsd_col = "jsd_mean",
                             eps     = 1e-6) {

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) ||
      eps <= 0 || eps >= 0.5) {
    stop("`eps` must be a single finite number between 0 and 0.5.", call. = FALSE)
  }
  if (!jsd_col %in% names(jsd_df)) {
    stop("`jsd_col` must be a column in `jsd_df`.")
  }
  jsd <- jsd_df[[jsd_col]]
  if (!is.numeric(jsd)) {
    stop("`jsd_col` must refer to a numeric column.", call. = FALSE)
  }
  bad <- !is.na(jsd) & (!is.finite(jsd) | jsd < 0 | jsd > 1)
  if (any(bad)) {
    stop(
      "`", jsd_col, "` must contain only values in [0, 1] or NA. ",
      "Invalid row(s): ", paste(which(bad), collapse = ", "),
      call. = FALSE
    )
  }
  jsd_beta <- pmin(pmax(jsd, eps), 1 - eps)
  jsd_df$jsd_beta <- jsd_beta
  jsd_df
}
