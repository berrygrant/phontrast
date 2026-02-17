#' Group-level JSD point estimates
#'
#' Computes JSD for each group (e.g., speaker) comparing two categories
#' (e.g., vowels) in an n-dimensional acoustic space.
#'
#' @param data Data frame containing acoustic measurements.
#' @param group_col String: name of column giving the grouping unit
#'   (e.g., "speaker").
#' @param category_col String: name of column giving the category
#'   to compare (e.g., "vowel"). Each group must have exactly two categories.
#' @param features Character vector of column names giving the acoustic space.
#' @param min_tokens Minimum number of tokens per group required to compute
#'   JSD. Groups with fewer tokens are dropped.
#' @param ... Additional arguments passed to \code{jsd_kde_nd()}.
#'
#' @return A tibble with one row per group and columns:
#'   \code{group}, \code{n_tokens}, and \code{jsd}.
#' @export
#' @importFrom dplyr n_distinct
#' @importFrom tibble tibble
speaker_jsd <- function(data,
                        group_col,
                        category_col,
                        features,
                        min_tokens = 20,
                        ...) {

  if (!group_col %in% names(data)) {
    stop("`group_col` must be a column in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("`category_col` must be a column in `data`.")
  }
  if (!all(features %in% names(data))) {
    stop("All `features` must be columns in `data`.")
  }

  groups <- split(data, data[[group_col]])

  out_list <- lapply(groups, function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens ||
        dplyr::n_distinct(df_g[[category_col]]) != 2L) {
      return(NULL)
    }

    jsd_val <- tryCatch(
      jsd_kde_nd(
        data     = df_g,
        features = features,
        group    = category_col,
        ...
      ),
      error = function(e) NA_real_
    )

    tibble::tibble(
      group    = df_g[[group_col]][1],
      n_tokens = n_tok,
      jsd      = jsd_val
    )
  })

  out <- do.call(rbind, out_list)
  if (is.null(out)) {
    out <- tibble::tibble(
      group    = character(0),
      n_tokens = integer(0),
      jsd      = numeric(0)
    )
  }
  out
}

#' Bootstrap JSD for each group
#'
#' Computes bootstrap mean, SD, and confidence interval for JSD within each
#' group (e.g., speaker), using resampling with replacement.
#'
#' @inheritParams speaker_jsd
#' @param n_boot Number of bootstrap resamples per group.
#' @param est_distance Logical; if TRUE, return Jensen–Shannon distance
#'   (sqrt of divergence) instead of divergence.
#'
#' @return A tibble with one row per group and columns:
#'   \code{group}, \code{n_tokens}, \code{jsd_mean}, \code{jsd_sd},
#'   \code{jsd_low}, and \code{jsd_high}.
#' @export
#' @importFrom purrr map
#' @importFrom dplyr bind_rows n_distinct
#' @importFrom stats quantile sd
boot_jsd <- function(data,
                     group_col,
                     category_col,
                     features,
                     n_boot     = 300,
                     min_tokens = 30,
                     est_distance = FALSE,
                     ...) {

  if (!group_col %in% names(data)) {
    stop("`group_col` must be a column in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("`category_col` must be a column in `data`.")
  }

  groups <- split(data, data[[group_col]])

  res_list <- purrr::map(groups, function(df_g) {

    if (nrow(df_g) < min_tokens ||
        dplyr::n_distinct(df_g[[category_col]]) < 2L) {
      return(tibble::tibble(
        group    = df_g[[group_col]][1],
        n_tokens = nrow(df_g),
        jsd_mean = NA_real_,
        jsd_sd   = NA_real_,
        jsd_low  = NA_real_,
        jsd_high = NA_real_
      ))
    }

    jsd_vals <- replicate(
      n_boot,
      {
        samp <- df_g[sample.int(nrow(df_g), size = nrow(df_g), replace = TRUE), ]
        if (dplyr::n_distinct(samp[[category_col]]) < 2L) {
          return(NA_real_)
        }
        jsd_div <- tryCatch(
          jsd_kde_nd(
            data     = samp,
            features = features,
            group    = category_col,
            ...
          ),
          error = function(e) NA_real_
        )
        if (!is.finite(jsd_div)) {
          return(NA_real_)
        }
        if (est_distance) sqrt(jsd_div) else jsd_div
      }
    )

    jsd_vals <- jsd_vals[!is.na(jsd_vals)]
    if (!length(jsd_vals)) {
      return(tibble::tibble(
        group    = df_g[[group_col]][1],
        n_tokens = nrow(df_g),
        jsd_mean = NA_real_,
        jsd_sd   = NA_real_,
        jsd_low  = NA_real_,
        jsd_high = NA_real_
      ))
    }

    tibble::tibble(
      group    = df_g[[group_col]][1],
      n_tokens = nrow(df_g),
      jsd_mean = mean(jsd_vals),
      jsd_sd   = stats::sd(jsd_vals),
      jsd_low  = stats::quantile(jsd_vals, 0.025, names = FALSE),
      jsd_high = stats::quantile(jsd_vals, 0.975, names = FALSE)
    )
  })

  dplyr::bind_rows(res_list)
}
