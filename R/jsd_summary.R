#' JSD summary: point estimate and optional bootstrap per group
#'
#' Convenience wrapper that returns both the point-estimate JSD and,
#' optionally, bootstrap-based uncertainty (mean, SD, and CI) for each group.
#'
#' @inheritParams speaker_jsd
#' @param do_boot Logical; if TRUE (default), perform bootstrap via \code{boot_jsd()}.
#' @param n_boot Integer; number of bootstrap resamples per group if
#'   \code{do_boot = TRUE}.
#' @param conf_level Confidence level for bootstrap intervals.
#'
#' @return A tibble with one row per group and columns:
#'   \itemize{
#'     \item \code{group} – group ID (e.g., speaker)
#'     \item \code{n_tokens} – number of tokens for that group
#'     \item \code{jsd_point} – single JSD point estimate
#'     \item \code{jsd_mean}, \code{jsd_sd}, \code{jsd_low}, \code{jsd_high}
#'       – bootstrap summary (if \code{do_boot = TRUE})
#'   }
#' @export
#' @importFrom dplyr left_join rename
#' @importFrom rlang .data
jsd_summary <- function(data,
                        group_col,
                        category_col,
                        features,
                        do_boot     = TRUE,
                        n_boot      = 300,
                        min_tokens  = 30,
                        conf_level  = 0.95,
                        ...) {

  .check_conf_level(conf_level)
  if (isTRUE(do_boot)) {
    .check_positive_count(n_boot, "n_boot")
  }
  .check_positive_count(min_tokens, "min_tokens")

  # Point estimates
  pt <- speaker_jsd(
    data         = data,
    group_col    = group_col,
    category_col = category_col,
    features     = features,
    min_tokens   = min_tokens,
    ...
  ) |>
    dplyr::rename(jsd_point = "jsd")

  if (!do_boot) {
    return(pt)
  }

  # Bootstrap estimates
  bt <- boot_jsd(
    data         = data,
    group_col    = group_col,
    category_col = category_col,
    features     = features,
    n_boot       = n_boot,
    min_tokens   = min_tokens,
    conf_level   = conf_level,
    ...
  )

  # Join on group + n_tokens (both functions report them)
  out <- dplyr::left_join(
    pt,
    bt,
    by = c("group", "n_tokens")
  )

  out
}
