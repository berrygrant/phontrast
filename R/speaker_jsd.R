#' Group-level JSD point estimates
#'
#' Computes JSD for each group (e.g., speaker) comparing two categories
#' (e.g., vowels) in an n-dimensional acoustic space.
#'
#' @param data Data frame containing acoustic measurements.
#' @param group_col Character vector giving one or more grouping columns
#'   (e.g., \code{"speaker"} or \code{c("Sex", "Style")}).
#' @param category_col String: name of column giving the category
#'   to compare (e.g., "vowel"). Each group must have exactly two categories.
#' @param features Character vector of column names giving the acoustic space.
#' @param min_tokens Minimum number of tokens per group required to compute
#'   JSD. Groups with fewer tokens are dropped.
#' @param bw Bandwidth selection method passed to \code{jsd_kde_nd()}.
#' @param eval_on KDE evaluation points passed to \code{jsd_kde_nd()}.
#' @param eval_n Optional maximum number of KDE evaluation points.
#' @param eval_seed Optional integer seed for KDE evaluation-point subsampling.
#' @param engine KDE evaluation engine passed to \code{jsd_kde_nd()}.
#'   \code{"fast_diagonal"} is accepted as an alias for \code{"fast_diag"}.
#' @param chunk_size Chunk size for \code{engine = "fast_diag"}.
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
                        bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                        eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                        eval_n = NULL,
                        eval_seed = NULL,
                        engine = c("ks", "fast_diag", "fast_diagonal"),
                        chunk_size = 1000L,
                        ...) {

  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  engine <- .match_kde_engine(engine)
  .check_positive_count(min_tokens, "min_tokens")
  group_col <- .check_group_cols(group_col)
  .check_columns(data, c(group_col, category_col, features))
  data <- .metric_data(data, c(group_col, category_col, features))

  groups <- .split_groups(data, group_col)

  out_list <- lapply(groups, function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens ||
        .observed_n_categories(df_g[[category_col]]) != 2L) {
      return(NULL)
    }

    jsd_val <- tryCatch(
      jsd_kde_nd(
        data     = df_g,
        features = features,
        group    = category_col,
        bw       = bw,
        eval_on  = eval_on,
        eval_n   = eval_n,
        eval_seed = eval_seed,
        engine   = engine,
        chunk_size = chunk_size,
        ...
      ),
      error = function(e) NA_real_
    )

    tibble::tibble(
      group    = .group_label(df_g, group_col),
      n_tokens = n_tok,
      jsd      = jsd_val
    )
  })

  out <- dplyr::bind_rows(out_list)
  if (nrow(out)) out else .empty_group_jsd()
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
#' @param conf_level Confidence level for bootstrap intervals.
#'
#' @return A tibble with one row per group and columns:
#'   \code{group}, \code{n_tokens}, \code{n_boot}, \code{conf_level},
#'   \code{jsd_mean}, \code{jsd_sd}, \code{ci_lower}, \code{ci_upper},
#'   \code{jsd_low}, and \code{jsd_high}.
#'   \code{jsd_low} and \code{jsd_high} are retained as legacy aliases for
#'   \code{ci_lower} and \code{ci_upper}.
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
                     conf_level = 0.95,
                     bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                     eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                     eval_n = NULL,
                     eval_seed = NULL,
                     engine = c("ks", "fast_diag", "fast_diagonal"),
                     chunk_size = 1000L,
                     ...) {

  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  engine <- .match_kde_engine(engine)
  .check_positive_count(n_boot, "n_boot")
  .check_positive_count(min_tokens, "min_tokens")
  .check_conf_level(conf_level)
  group_col <- .check_group_cols(group_col)
  .check_columns(data, c(group_col, category_col, features))
  data <- .metric_data(data, c(group_col, category_col, features))
  alpha <- 1 - conf_level

  groups <- .split_groups(data, group_col)

  res_list <- purrr::map(groups, function(df_g) {

    if (nrow(df_g) < min_tokens ||
        .observed_n_categories(df_g[[category_col]]) != 2L) {
      return(tibble::tibble(
        group    = .group_label(df_g, group_col),
        n_tokens = nrow(df_g),
        n_boot   = 0L,
        conf_level = conf_level,
        jsd_mean = NA_real_,
        jsd_sd   = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        jsd_low  = NA_real_,
        jsd_high = NA_real_
      ))
    }

    jsd_vals <- vapply(
      seq_len(n_boot),
      function(i) {
        samp <- df_g[sample.int(nrow(df_g), size = nrow(df_g), replace = TRUE), ]
        if (.observed_n_categories(samp[[category_col]]) != 2L) {
          return(NA_real_)
        }
        jsd_div <- tryCatch(
          jsd_kde_nd(
            data     = samp,
            features = features,
            group    = category_col,
            bw       = bw,
            eval_on  = eval_on,
            eval_n   = eval_n,
            eval_seed = eval_seed,
            engine   = engine,
            chunk_size = chunk_size,
            ...
          ),
          error = function(e) NA_real_
        )
        if (!is.finite(jsd_div)) {
          return(NA_real_)
        }
        if (est_distance) sqrt(jsd_div) else jsd_div
      },
      numeric(1)
    )

    jsd_vals <- jsd_vals[!is.na(jsd_vals)]
    if (!length(jsd_vals)) {
      return(tibble::tibble(
        group    = .group_label(df_g, group_col),
        n_tokens = nrow(df_g),
        n_boot   = 0L,
        conf_level = conf_level,
        jsd_mean = NA_real_,
        jsd_sd   = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        jsd_low  = NA_real_,
        jsd_high = NA_real_
      ))
    }

    qs <- stats::quantile(
      jsd_vals,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE
    )

    tibble::tibble(
      group    = .group_label(df_g, group_col),
      n_tokens = nrow(df_g),
      n_boot   = length(jsd_vals),
      conf_level = conf_level,
      jsd_mean = mean(jsd_vals),
      jsd_sd   = stats::sd(jsd_vals),
      ci_lower = qs[1],
      ci_upper = qs[2],
      jsd_low  = qs[1],
      jsd_high = qs[2]
    )
  })

  dplyr::bind_rows(res_list)
}
