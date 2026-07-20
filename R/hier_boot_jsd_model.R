#' Hierarchical bootstrap for JSD-based models
#'
#' Performs a hierarchical bootstrap: resample groups with replacement,
#' resample tokens within each sampled group, compute JSD per group,
#' fit a model to the bootstrap JSD values, and repeat.
#'
#' This lets you propagate measurement uncertainty in JSD into model
#' parameters (e.g., GAM/LMM coefficients).
#'
#' @param data Data frame with at least: group_col, category_col, features,
#'   and any predictors used in the model.
#' @param group_col String: grouping variable (e.g., "speaker").
#' @param category_col String: category variable with 2 levels (e.g., "vowel").
#' @param features Character vector of acoustic feature columns.
#' @param formula Model formula to pass to `fit_fun` (e.g.,
#'   `jsd_beta ~ s(age) + s(region, bs = "re")`).
#' @param fit_fun A function that takes `(formula, data, ...)` and returns a
#'   fitted model. Defaults to `mgcv::gam` if available, otherwise `stats::lm`.
#' @param n_outer Number of hierarchical bootstrap replicates.
#' @param min_tokens Minimum within-group tokens required.
#' @param eps Small epsilon for bounding JSD in (0, 1) if using Beta family.
#' @param progress Logical; if TRUE, prints progress every 10 replicates.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{boot_id} – bootstrap replicate index
#'     \item \code{term} – model term
#'     \item \code{estimate} – estimate for that term in that replicate
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(2026)
#' speakers <- paste0("s", 1:4)
#' dat <- data.frame(
#'   speaker = rep(speakers, each = 60),
#'   age = rep(c(22, 35, 48, 61), each = 60),
#'   vowel = rep(rep(c("ih", "eh"), each = 30), 4)
#' )
#' dat$f1 <- rnorm(
#'   nrow(dat),
#'   mean = ifelse(dat$vowel == "ih", 500, 560) + dat$age * 0.3,
#'   sd = 55
#' )
#'
#' hier_boot_jsd_model(
#'   data = dat,
#'   group_col = "speaker",
#'   category_col = "vowel",
#'   features = "f1",
#'   formula = jsd_beta ~ age,
#'   fit_fun = stats::lm,
#'   n_outer = 3,
#'   min_tokens = 20,
#'   progress = FALSE
#' )
#' }
#' @export
#' @importFrom dplyr group_by summarize first bind_rows left_join across all_of
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom stats coef
#' @importFrom rlang .data
hier_boot_jsd_model <- function(data,
                                group_col,
                                category_col,
                                features,
                                formula,
                                fit_fun   = NULL,
                                n_outer   = 200,
                                min_tokens = 20,
                                eps       = 1e-6,
                                progress  = TRUE,
                                ...) {

  .check_columns(data, c(group_col, category_col, features))
  .check_positive_count(n_outer, "n_outer")
  .check_positive_count(min_tokens, "min_tokens")
  .validate_metric_inputs(data, features, category_col, group_col)
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a model formula.", call. = FALSE)
  }

  if (is.null(fit_fun)) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
      fit_fun <- mgcv::gam
    } else {
      fit_fun <- stats::lm
    }
  }
  if (!is.function(fit_fun)) {
    stop("`fit_fun` must be a function.", call. = FALSE)
  }

  metric_cols <- c(group_col, category_col, features)
  data <- data[stats::complete.cases(data[, metric_cols, drop = FALSE]), , drop = FALSE]
  .check_numeric_features(data, features)
  finite_feature_rows <- Reduce(
    `&`,
    lapply(data[, features, drop = FALSE], is.finite)
  )
  data <- data[finite_feature_rows, , drop = FALSE]
  groups <- unique(data[[group_col]])
  n_groups <- length(groups)
  if (!n_groups) {
    return(tibble::tibble(
      boot_id = integer(),
      term = character(),
      estimate = numeric()
    ))
  }

  boot_group_col <- ".phonJSD_boot_group"
  original_group_col <- ".phonJSD_original_group"
  while (boot_group_col %in% names(data)) {
    boot_group_col <- paste0(boot_group_col, "_")
  }
  while (original_group_col %in% names(data)) {
    original_group_col <- paste0(original_group_col, "_")
  }

  preds <- data |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarize(
      dplyr::across(
        .cols = !dplyr::all_of(c(category_col, features)),
        .fns  = dplyr::first
      ),
      .groups = "drop"
    )

  boot_results <- vector("list", n_outer)

  for (b in seq_len(n_outer)) {
    if (progress && b %% 10 == 0) {
      message("Bootstrap replicate ", b, " / ", n_outer)
    }

    # Resample groups with replacement. Index into `groups` rather than calling
    # sample(groups, ...) directly: with a single numeric group id, sample() would
    # dispatch to sample.int() and draw from 1:id instead of the id itself.
    boot_group_ids <- groups[sample.int(n_groups, size = n_groups, replace = TRUE)]

    boot_df_list <- Map(function(g, draw_id) {
      df_g <- data[data[[group_col]] == g, , drop = FALSE]
      if (nrow(df_g) < min_tokens ||
          dplyr::n_distinct(df_g[[category_col]]) < 2L) {
        return(NULL)
      }
      # Resample tokens within group
      samp <- df_g[sample.int(nrow(df_g), size = nrow(df_g), replace = TRUE), ]
      samp[[boot_group_col]] <- paste0("boot", b, "_draw", draw_id)
      samp[[original_group_col]] <- g
      samp
    }, boot_group_ids, seq_along(boot_group_ids))

    boot_df <- dplyr::bind_rows(boot_df_list)
    if (nrow(boot_df) == 0L) {
      next
    }

    # Compute JSD per group for this bootstrap sample
    jsd_sum <- jsd_summary(
      data         = boot_df,
      group_col    = boot_group_col,
      category_col = category_col,
      features     = features,
      do_boot      = FALSE,
      min_tokens   = min_tokens
    )
    if (!nrow(jsd_sum)) {
      next
    }

    group_map <- boot_df |>
      dplyr::group_by(.data[[boot_group_col]]) |>
      dplyr::summarize(
        original_group = dplyr::first(.data[[original_group_col]]),
        .groups = "drop"
      )
    names(group_map)[names(group_map) == boot_group_col] <- "group"

    model_df <- jsd_sum |>
      dplyr::left_join(group_map, by = "group") |>
      dplyr::left_join(preds, by = c("original_group" = group_col))
    model_df[[group_col]] <- model_df$original_group
    model_df <- model_df[is.finite(model_df$jsd_point), , drop = FALSE]
    if (!nrow(model_df)) {
      next
    }

    # Prepare jsd_beta using point JSD
    model_df <- prepare_jsd_beta(
      jsd_df = model_df,
      jsd_col = "jsd_point",
      eps = eps
    )

    # Fit model
    fit <- try(fit_fun(formula = formula, data = model_df, ...),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      next
    }

    cf <- stats::coef(fit)
    boot_results[[b]] <- tibble::tibble(
      boot_id  = b,
      term     = names(cf),
      estimate = unname(cf)
    )
  }

  out <- dplyr::bind_rows(boot_results)
  if (nrow(out)) {
    out
  } else {
    tibble::tibble(
      boot_id = integer(),
      term = character(),
      estimate = numeric()
    )
  }
}
