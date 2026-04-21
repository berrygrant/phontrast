## ============================================================
## JSD estimator robustness check
##
## Computes the same contrast-level JSD under several distribution
## estimators so the substantive pattern can be evaluated independently
## of KDE-specific choices.
## ============================================================

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- args_full[grepl("^--file=", args_full)]
  if (length(file_arg) > 0) {
    p <- sub("^--file=", "", file_arg[1])
    return(normalizePath(dirname(p)))
  }
  p <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (!is.na(p)) return(normalizePath(dirname(p)))
  normalizePath(getwd())
}

script_dir <- get_script_dir()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

load_phonjsd <- function(script_dir) {
  pkg_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  use_local <- FALSE
  want_local <- tolower(Sys.getenv("ROBUST_JSD_LOAD_LOCAL", unset = "true")) %in%
    c("1", "true", "yes")
  if (want_local &&
      file.exists(file.path(pkg_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    try({
      devtools::load_all(pkg_root, quiet = TRUE)
      use_local <- TRUE
    }, silent = TRUE)
  }
  if (!use_local) {
    suppressPackageStartupMessages(library(phonJSD))
  }
  use_local
}

using_local_source <- load_phonjsd(script_dir)

parse_int_env <- function(name, default, min_value = 0L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  val <- suppressWarnings(as.integer(raw))
  if (is.na(val) || val < min_value) {
    warning("Invalid ", name, "='", raw, "'. Using ", default, ".")
    return(as.integer(default))
  }
  as.integer(val)
}

parse_num_env <- function(name, default, min_value = -Inf) {
  raw <- Sys.getenv(name, unset = as.character(default))
  val <- suppressWarnings(as.numeric(raw))
  if (is.na(val) || val < min_value) {
    warning("Invalid ", name, "='", raw, "'. Using ", default, ".")
    return(as.numeric(default))
  }
  as.numeric(val)
}

safe_read <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

as_numeric_matrix <- function(df, features) {
  X <- as.matrix(df[, features, drop = FALSE])
  storage.mode(X) <- "double"
  if (is.null(dim(X))) {
    X <- matrix(X, ncol = length(features))
  }
  X
}

normalize_prob <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x <- pmax(x, 0)
  s <- sum(x)
  if (!is.finite(s) || s <= 0) {
    return(rep(NA_real_, length(x)))
  }
  x / s
}

normalize_log_prob <- function(log_x) {
  log_x <- as.numeric(log_x)
  ok <- is.finite(log_x)
  if (!any(ok)) return(rep(NA_real_, length(log_x)))
  m <- max(log_x[ok])
  out <- rep(0, length(log_x))
  out[ok] <- exp(log_x[ok] - m)
  normalize_prob(out)
}

log_sum_exp <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (!any(ok)) return(-Inf)
  m <- max(x[ok])
  m + log(sum(exp(x[ok] - m)))
}

sample_eval_points <- function(X, eval_n) {
  X <- as.matrix(X)
  if (eval_n > 0L && nrow(X) > eval_n) {
    X <- X[sample.int(nrow(X), eval_n), , drop = FALSE]
  }
  X
}

## ============================================================
## KDE estimator
## ============================================================

jsd_kde_1d <- function(x1, x2, n_grid = 512L) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 2L || length(x2) < 2L) {
    stop("Not enough finite values for 1D KDE JSD.")
  }
  rng <- range(c(x1, x2))
  if (!is.finite(rng[1]) || !is.finite(rng[2])) {
    stop("Non-finite 1D range.")
  }
  if (rng[1] == rng[2]) return(0)

  d1 <- stats::density(x1, n = n_grid, from = rng[1], to = rng[2])
  d2 <- stats::density(x2, n = n_grid, from = rng[1], to = rng[2])
  p <- normalize_prob(d1$y)
  q <- normalize_prob(d2$y)
  jsd(p, q)
}

jsd_kde_point <- function(data,
                          features,
                          category_col,
                          bw = "Hpi.diag",
                          eval_on = "pooled") {
  levs <- unique(data[[category_col]])
  if (length(levs) != 2L) stop("Need exactly two categories for KDE JSD.")

  if (length(features) == 1L) {
    x1 <- data[data[[category_col]] == levs[1], features[1], drop = TRUE]
    x2 <- data[data[[category_col]] == levs[2], features[1], drop = TRUE]
    return(jsd_kde_1d(x1, x2))
  }

  jsd_kde_nd(
    data = data,
    features = features,
    group = category_col,
    bw = bw,
    eval_on = eval_on
  )
}

## ============================================================
## Binned estimators
## ============================================================

make_binner <- function(X, bins, type = c("equal_width", "quantile")) {
  type <- match.arg(type)
  X <- as.matrix(X)
  d <- ncol(X)
  breaks <- vector("list", d)
  bins_dim <- integer(d)

  for (j in seq_len(d)) {
    x <- X[, j]
    x <- x[is.finite(x)]
    if (!length(x)) {
      breaks[[j]] <- c(0, 1)
      bins_dim[j] <- 1L
      next
    }

    rng <- range(x)
    if (rng[1] == rng[2]) {
      delta <- max(abs(rng[1]), 1) * 1e-8
      breaks[[j]] <- c(rng[1] - delta, rng[2] + delta)
      bins_dim[j] <- 1L
      next
    }

    if (type == "equal_width") {
      br <- seq(rng[1], rng[2], length.out = bins + 1L)
    } else {
      br <- stats::quantile(
        x,
        probs = seq(0, 1, length.out = bins + 1L),
        names = FALSE,
        type = 8
      )
      br <- unique(as.numeric(br))
      if (length(br) < 2L) {
        delta <- max(abs(rng[1]), 1) * 1e-8
        br <- c(rng[1] - delta, rng[2] + delta)
      }
    }

    br[1] <- br[1] - max(abs(br[1]), 1) * 1e-12
    br[length(br)] <- br[length(br)] + max(abs(br[length(br)]), 1) * 1e-12
    breaks[[j]] <- br
    bins_dim[j] <- length(br) - 1L
  }

  list(breaks = breaks, bins_dim = bins_dim, type = type)
}

apply_binner <- function(X, binner) {
  X <- as.matrix(X)
  codes <- matrix(1L, nrow = nrow(X), ncol = ncol(X))
  for (j in seq_len(ncol(X))) {
    codes[, j] <- as.integer(cut(
      X[, j],
      breaks = binner$breaks[[j]],
      include.lowest = TRUE,
      labels = FALSE
    ))
    codes[is.na(codes[, j]), j] <- 0L
  }

  if (ncol(codes) == 1L) {
    return(as.character(codes[, 1]))
  }
  apply(codes, 1L, paste, collapse = ":")
}

jsd_binned <- function(X1,
                       X2,
                       bins,
                       type = c("equal_width", "quantile"),
                       smooth = 0) {
  type <- match.arg(type)
  X_all <- rbind(X1, X2)
  binner <- make_binner(X_all, bins = bins, type = type)
  keys1 <- apply_binner(X1, binner)
  keys2 <- apply_binner(X2, binner)
  keys <- union(unique(keys1), unique(keys2))
  c1 <- tabulate(match(keys1, keys), nbins = length(keys))
  c2 <- tabulate(match(keys2, keys), nbins = length(keys))

  p <- normalize_prob(c1 + smooth)
  q <- normalize_prob(c2 + smooth)
  if (anyNA(p) || anyNA(q)) stop("Binned estimator produced empty probabilities.")

  list(
    jsd = jsd(p, q),
    bins_per_feature = paste(binner$bins_dim, collapse = ","),
    possible_bins = prod(as.numeric(pmax(binner$bins_dim, 1L))),
    observed_bins = length(keys)
  )
}

## ============================================================
## Gaussian mixture estimator
## ============================================================

jsd_gmm <- function(X1,
                    X2,
                    max_components = 3L,
                    eval_n = 600L) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package `mclust` is required for Gaussian-mixture JSD.")
  }
  suppressPackageStartupMessages(
    require("mclust", quietly = TRUE, character.only = TRUE)
  )

  fit_one <- function(X) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    X <- X[stats::complete.cases(X), , drop = FALSE]
    if (nrow(X) < 2L) {
      stop("Need at least two complete rows to fit a Gaussian mixture.")
    }
    unique_n <- nrow(unique(as.data.frame(X)))
    g <- seq_len(max(1L, min(as.integer(max_components), nrow(X), unique_n)))
    fit <- mclust::Mclust(data = X, G = g, verbose = FALSE)
    if (is.null(fit) || is.null(fit$parameters)) {
      stop("mclust could not fit a Gaussian mixture.")
    }
    fit
  }

  fit1 <- fit_one(X1)
  fit2 <- fit_one(X2)
  X_eval <- sample_eval_points(rbind(X1, X2), eval_n = eval_n)

  logp <- mclust::dens(
    modelName = fit1$modelName,
    data = X_eval,
    parameters = fit1$parameters,
    logarithm = TRUE
  )
  logq <- mclust::dens(
    modelName = fit2$modelName,
    data = X_eval,
    parameters = fit2$parameters,
    logarithm = TRUE
  )

  p <- normalize_log_prob(logp)
  q <- normalize_log_prob(logq)
  if (anyNA(p) || anyNA(q)) stop("mclust produced empty probabilities.")

  list(
    jsd = jsd(p, q),
    engine = "mclust",
    n_eval = nrow(X_eval),
    k_cat1 = fit1$G,
    k_cat2 = fit2$G,
    bic_cat1 = as.numeric(fit1$bic),
    bic_cat2 = as.numeric(fit2$bic),
    model_cat1 = fit1$modelName,
    model_cat2 = fit2$modelName
  )
}

## ============================================================
## Method orchestration
## ============================================================

ok_row <- function(method,
                   estimator_family,
                   res,
                   detail = NA_character_) {
  tibble(
    method = method,
    estimator_family = estimator_family,
    jsd = as.numeric(res$jsd),
    n_eval = if (!is.null(res$n_eval)) as.integer(res$n_eval) else NA_integer_,
    bins_per_feature = if (!is.null(res$bins_per_feature)) res$bins_per_feature else NA_character_,
    possible_bins = if (!is.null(res$possible_bins)) as.numeric(res$possible_bins) else NA_real_,
    observed_bins = if (!is.null(res$observed_bins)) as.integer(res$observed_bins) else NA_integer_,
    gmm_engine = if (!is.null(res$engine)) res$engine else NA_character_,
    gmm_k_cat1 = if (!is.null(res$k_cat1)) as.integer(res$k_cat1) else NA_integer_,
    gmm_k_cat2 = if (!is.null(res$k_cat2)) as.integer(res$k_cat2) else NA_integer_,
    gmm_model_cat1 = if (!is.null(res$model_cat1)) res$model_cat1 else NA_character_,
    gmm_model_cat2 = if (!is.null(res$model_cat2)) res$model_cat2 else NA_character_,
    gmm_bic_cat1 = if (!is.null(res$bic_cat1)) as.numeric(res$bic_cat1) else NA_real_,
    gmm_bic_cat2 = if (!is.null(res$bic_cat2)) as.numeric(res$bic_cat2) else NA_real_,
    detail = detail,
    status = "ok",
    note = NA_character_
  )
}

error_row <- function(method,
                      estimator_family,
                      error,
                      detail = NA_character_) {
  tibble(
    method = method,
    estimator_family = estimator_family,
    jsd = NA_real_,
    n_eval = NA_integer_,
    bins_per_feature = NA_character_,
    possible_bins = NA_real_,
    observed_bins = NA_integer_,
    gmm_engine = NA_character_,
    gmm_k_cat1 = NA_integer_,
    gmm_k_cat2 = NA_integer_,
    gmm_model_cat1 = NA_character_,
    gmm_model_cat2 = NA_character_,
    gmm_bic_cat1 = NA_real_,
    gmm_bic_cat2 = NA_real_,
    detail = detail,
    status = "error",
    note = conditionMessage(error)
  )
}

method_try <- function(method,
                       estimator_family,
                       expr,
                       detail = NA_character_) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    error_row(method, estimator_family, res, detail = detail)
  } else {
    ok_row(method, estimator_family, res, detail = detail)
  }
}

compute_jsd_methods <- function(data,
                                features,
                                category_col,
                                bw,
                                eval_on,
                                hist_bins,
                                empirical_bins,
                                gmm_components,
                                gmm_eval_n,
                                bin_smooth = 0) {
  levs <- unique(data[[category_col]])
  if (length(levs) != 2L) stop("Need exactly two categories.")

  X1 <- as_numeric_matrix(data[data[[category_col]] == levs[1], , drop = FALSE], features)
  X2 <- as_numeric_matrix(data[data[[category_col]] == levs[2], , drop = FALSE], features)

  dplyr::bind_rows(
    method_try(
      method = "kde",
      estimator_family = "kernel_density",
      expr = list(
        jsd = jsd_kde_point(
          data = data,
          features = features,
          category_col = category_col,
          bw = bw,
          eval_on = eval_on
        )
      ),
      detail = paste0("bw=", bw, "; eval_on=", eval_on)
    ),
    method_try(
      method = "histogram_equal_width",
      estimator_family = "binned_histogram",
      expr = jsd_binned(
        X1 = X1,
        X2 = X2,
        bins = hist_bins,
        type = "equal_width",
        smooth = bin_smooth
      ),
      detail = paste0("bins=", hist_bins, "; smooth=", bin_smooth)
    ),
    method_try(
      method = "empirical_quantile_binned",
      estimator_family = "empirical_discrete",
      expr = jsd_binned(
        X1 = X1,
        X2 = X2,
        bins = empirical_bins,
        type = "quantile",
        smooth = bin_smooth
      ),
      detail = paste0("quantile_bins=", empirical_bins, "; smooth=", bin_smooth)
    ),
    method_try(
      method = "gaussian_mixture",
      estimator_family = "gaussian_mixture",
      expr = jsd_gmm(
        X1 = X1,
        X2 = X2,
        max_components = gmm_components,
        eval_n = gmm_eval_n
      ),
      detail = paste0("max_components=", gmm_components, "; covariance=BIC-selected")
    )
  )
}

wide_from_long <- function(results_df) {
  ok <- results_df %>%
    filter(.data$status == "ok") %>%
    select(contrast_id, method, jsd)
  if (!nrow(ok)) return(tibble())

  wide <- stats::reshape(
    as.data.frame(ok),
    idvar = "contrast_id",
    timevar = "method",
    direction = "wide"
  )
  names(wide) <- sub("^jsd\\.", "jsd_", names(wide))
  meta <- results_df %>%
    distinct(
      contrast_id,
      domain,
      cat1,
      cat2,
      expected_sep,
      n_cat1,
      n_cat2,
      feature_set
    )
  dplyr::left_join(meta, tibble::as_tibble(wide), by = "contrast_id")
}

summarize_concordance <- function(results_df) {
  ok <- results_df %>% filter(.data$status == "ok", is.finite(.data$jsd))
  if (!nrow(ok)) return(tibble())

  rows <- lapply(split(ok, ok$method), function(df_m) {
    low <- df_m$jsd[df_m$expected_sep == 0]
    high <- df_m$jsd[df_m$expected_sep == 1]
    rho <- if (length(unique(df_m$expected_sep)) == 2L && nrow(df_m) >= 3L) {
      suppressWarnings(stats::cor(df_m$jsd, df_m$expected_sep, method = "spearman"))
    } else {
      NA_real_
    }
    r <- if (length(unique(df_m$expected_sep)) == 2L && nrow(df_m) >= 3L) {
      suppressWarnings(stats::cor(df_m$jsd, df_m$expected_sep, method = "pearson"))
    } else {
      NA_real_
    }
    all_above <- length(low) > 0L && length(high) > 0L &&
      isTRUE(min(high, na.rm = TRUE) > max(low, na.rm = TRUE))

    tibble(
      method = df_m$method[1],
      estimator_family = df_m$estimator_family[1],
      n_contrasts = nrow(df_m),
      mean_stop_voicing = if (length(low)) mean(low) else NA_real_,
      mean_fricative_place = if (length(high)) mean(high) else NA_real_,
      mean_difference_fricative_minus_stop = if (length(low) && length(high)) {
        mean(high) - mean(low)
      } else {
        NA_real_
      },
      max_stop_voicing = if (length(low)) max(low) else NA_real_,
      min_fricative_place = if (length(high)) min(high) else NA_real_,
      all_fricatives_above_all_stops = all_above,
      spearman_with_expected_sep = rho,
      pearson_with_expected_sep = r,
      pattern_holds = all_above && isTRUE(mean(high) > mean(low))
    )
  })

  dplyr::bind_rows(rows) %>%
    arrange(.data$method)
}

summarize_method_correlations <- function(wide_df) {
  if (!nrow(wide_df)) return(tibble())
  method_cols <- grep("^jsd_", names(wide_df), value = TRUE)
  if (length(method_cols) < 2L) return(tibble())

  pairs <- utils::combn(method_cols, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(cols) {
    x <- wide_df[[cols[1]]]
    y <- wide_df[[cols[2]]]
    ok <- is.finite(x) & is.finite(y)
    tibble(
      method_1 = sub("^jsd_", "", cols[1]),
      method_2 = sub("^jsd_", "", cols[2]),
      n = sum(ok),
      pearson = if (sum(ok) >= 3L) stats::cor(x[ok], y[ok]) else NA_real_,
      spearman = if (sum(ok) >= 3L) {
        suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
      } else {
        NA_real_
      },
      mean_abs_diff = if (sum(ok) > 0L) mean(abs(x[ok] - y[ok])) else NA_real_
    )
  })
  dplyr::bind_rows(rows)
}

write_summary_txt <- function(path, concordance_df, corr_df, skipped_df) {
  if (!nrow(concordance_df)) {
    writeLines("No valid concordance results were produced.", con = path)
    return(invisible(path))
  }

  n_hold <- sum(concordance_df$pattern_holds, na.rm = TRUE)
  n_methods <- nrow(concordance_df)
  lines <- c(
    "JSD estimator robustness summary",
    "",
    paste0(n_hold, " of ", n_methods, " estimators preserve the expected pattern."),
    "Pattern criterion: every fricative-place contrast has higher JSD than every stop-voicing contrast.",
    "",
    "Estimator concordance:"
  )

  for (i in seq_len(nrow(concordance_df))) {
    row <- concordance_df[i, ]
    lines <- c(
      lines,
      paste0(
        "- ", row$method,
        ": pattern_holds=", row$pattern_holds,
        "; mean_stop=", signif(row$mean_stop_voicing, 4),
        "; mean_fricative=", signif(row$mean_fricative_place, 4),
        "; spearman_expected=", signif(row$spearman_with_expected_sep, 4)
      )
    )
  }

  if (nrow(corr_df)) {
    lines <- c(lines, "", "Pairwise method correlations:")
    for (i in seq_len(nrow(corr_df))) {
      row <- corr_df[i, ]
      lines <- c(
        lines,
        paste0(
          "- ", row$method_1, " vs ", row$method_2,
          ": n=", row$n,
          "; spearman=", signif(row$spearman, 4),
          "; mean_abs_diff=", signif(row$mean_abs_diff, 4)
        )
      )
    }
  }

  if (nrow(skipped_df)) {
    lines <- c(lines, "", paste0("Skipped contrasts: ", nrow(skipped_df)))
  }

  writeLines(lines, con = path)
  invisible(path)
}

maybe_write_plot <- function(path, results_df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 is not installed; skipping robustness plot.")
    return(invisible(FALSE))
  }

  plot_df <- results_df %>%
    filter(.data$status == "ok", is.finite(.data$jsd)) %>%
    mutate(
      contrast = paste(.data$cat1, .data$cat2, sep = "-"),
      expected_class = ifelse(.data$expected_sep == 1, "expected high", "expected low")
    )

  if (!nrow(plot_df)) return(invisible(FALSE))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = contrast, y = jsd, color = expected_class)
  ) +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(
      x = NULL,
      y = "Jensen-Shannon divergence",
      color = NULL,
      title = "JSD estimator robustness by contrast"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  ggplot2::ggsave(path, p, width = 10, height = 6, dpi = 300)
  invisible(TRUE)
}

## ============================================================
## Configuration
## ============================================================

set.seed(parse_int_env("ROBUST_JSD_SEED", 20260421L, min_value = 1L))

full_alignment_path <- file.path(script_dir, "sbcae_classical_consonant_metrics_full_alignments.csv")
regular_path <- file.path(script_dir, "sbcae_classical_consonant_metrics.csv")
default_data_path <- if (file.exists(full_alignment_path)) full_alignment_path else regular_path

data_path <- Sys.getenv("ROBUST_JSD_DATA_PATH", unset = default_data_path)
out_dir <- Sys.getenv(
  "ROBUST_JSD_OUT_DIR",
  unset = file.path(script_dir, "estimator_robustness_outputs")
)

bw <- Sys.getenv("ROBUST_JSD_BW", unset = "Hpi.diag")
eval_on <- Sys.getenv("ROBUST_JSD_EVAL_ON", unset = "pooled")
min_per_category <- parse_int_env("ROBUST_JSD_MIN_PER_CATEGORY", 30L, min_value = 2L)
max_per_category <- parse_int_env("ROBUST_JSD_MAX_PER_CATEGORY", 300L, min_value = 10L)
hist_bins_1d <- parse_int_env("ROBUST_JSD_HIST_BINS_1D", 32L, min_value = 2L)
hist_bins_nd <- parse_int_env("ROBUST_JSD_HIST_BINS_ND", 4L, min_value = 2L)
emp_bins_1d <- parse_int_env("ROBUST_JSD_EMP_BINS_1D", 32L, min_value = 2L)
emp_bins_nd <- parse_int_env("ROBUST_JSD_EMP_BINS_ND", 4L, min_value = 2L)
bin_smooth <- parse_num_env("ROBUST_JSD_BIN_SMOOTH", 0, min_value = 0)
gmm_components <- parse_int_env("ROBUST_JSD_GMM_COMPONENTS", 3L, min_value = 1L)
gmm_eval_n <- parse_int_env("ROBUST_JSD_GMM_EVAL_N", 600L, min_value = 50L)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ============================================================
## Data and contrasts
## ============================================================

df <- safe_read(data_path)
if (!all(c("segment", "speaker") %in% names(df))) {
  stop("Input must include at least columns: segment, speaker")
}

for (nm in setdiff(names(df), c("file", "speaker", "segment", "raw_label"))) {
  if (is.character(df[[nm]])) {
    suppressWarnings(df[[nm]] <- as.numeric(df[[nm]]))
  }
}

fric_features <- c(
  "cog_hz",
  "spec_sd_hz",
  "spec_skew",
  "spec_kurt",
  "peak_hz",
  "spec_slope_db_per_khz",
  "band_ratio_hi_lo_db",
  "intensity_db"
)

contrast_specs <- list(
  list(contrast_id = "stop_voicing_b_p", domain = "stop_voicing", cat1 = "b", cat2 = "p", features = c("vot_ms")),
  list(contrast_id = "stop_voicing_d_t", domain = "stop_voicing", cat1 = "d", cat2 = "t", features = c("vot_ms")),
  list(contrast_id = "stop_voicing_g_k", domain = "stop_voicing", cat1 = "g", cat2 = "k", features = c("vot_ms")),
  list(contrast_id = "fricative_place_s_sh", domain = "fricative_place", cat1 = "s", cat2 = "sh", features = fric_features),
  list(contrast_id = "fricative_place_f_th", domain = "fricative_place", cat1 = "f", cat2 = "th", features = fric_features),
  list(contrast_id = "fricative_place_z_zh", domain = "fricative_place", cat1 = "z", cat2 = "zh", features = fric_features)
)

result_rows <- list()
skipped_rows <- list()

for (spec in contrast_specs) {
  feat <- spec$features
  missing_feat <- setdiff(feat, names(df))
  if (length(missing_feat)) {
    skipped_rows[[length(skipped_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = paste0("Missing features: ", paste(missing_feat, collapse = ", "))
    )
    next
  }

  d <- df %>%
    filter(.data$segment %in% c(spec$cat1, spec$cat2)) %>%
    mutate(category = as.character(.data$segment)) %>%
    select(all_of(c("category", feat))) %>%
    filter(stats::complete.cases(.)) %>%
    group_by(.data$category) %>%
    group_modify(~ {
      if (nrow(.x) > max_per_category) {
        slice_sample(.x, n = max_per_category)
      } else {
        .x
      }
    }) %>%
    ungroup()

  n1 <- sum(d$category == spec$cat1)
  n2 <- sum(d$category == spec$cat2)
  if (n1 < min_per_category || n2 < min_per_category) {
    skipped_rows[[length(skipped_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = paste0("Too few tokens after filtering: ", spec$cat1, "=", n1, ", ", spec$cat2, "=", n2)
    )
    next
  }

  bins_hist <- if (length(feat) == 1L) hist_bins_1d else hist_bins_nd
  bins_emp <- if (length(feat) == 1L) emp_bins_1d else emp_bins_nd

  method_df <- tryCatch(
    compute_jsd_methods(
      data = d,
      features = feat,
      category_col = "category",
      bw = bw,
      eval_on = eval_on,
      hist_bins = bins_hist,
      empirical_bins = bins_emp,
      gmm_components = gmm_components,
      gmm_eval_n = gmm_eval_n,
      bin_smooth = bin_smooth
    ),
    error = function(e) {
      error_row(
        method = "all_methods",
        estimator_family = "all_methods",
        error = e
      )
    }
  )

  result_rows[[length(result_rows) + 1L]] <- method_df %>%
    mutate(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      expected_sep = ifelse(spec$domain == "fricative_place", 1, 0),
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      n_cat1 = n1,
      n_cat2 = n2,
      feature_set = paste(feat, collapse = ","),
      .before = 1
    )
}

results_df <- if (length(result_rows)) bind_rows(result_rows) else tibble()
skipped_df <- if (length(skipped_rows)) bind_rows(skipped_rows) else tibble()
wide_df <- wide_from_long(results_df)
concordance_df <- summarize_concordance(results_df)
method_corr_df <- summarize_method_correlations(wide_df)

results_path <- file.path(out_dir, "robust_jsd_by_method.csv")
wide_path <- file.path(out_dir, "robust_jsd_wide.csv")
concordance_path <- file.path(out_dir, "robust_jsd_method_concordance.csv")
corr_path <- file.path(out_dir, "robust_jsd_method_correlations.csv")
skipped_path <- file.path(out_dir, "robust_jsd_skipped_contrasts.csv")
summary_path <- file.path(out_dir, "robust_jsd_summary.txt")
plot_path <- file.path(out_dir, "robust_jsd_by_method.png")

write.csv(results_df, results_path, row.names = FALSE)
write.csv(wide_df, wide_path, row.names = FALSE)
write.csv(concordance_df, concordance_path, row.names = FALSE)
write.csv(method_corr_df, corr_path, row.names = FALSE)
write.csv(skipped_df, skipped_path, row.names = FALSE)
write_summary_txt(summary_path, concordance_df, method_corr_df, skipped_df)
maybe_write_plot(plot_path, results_df)

message("Finished JSD estimator robustness check.")
message("Using local source: ", using_local_source)
message("Input data: ", normalizePath(data_path, mustWork = FALSE))
message("By-method results: ", results_path)
message("Wide results: ", wide_path)
message("Method concordance: ", concordance_path)
message("Method correlations: ", corr_path)
message("Skipped contrasts: ", skipped_path)
message("Summary: ", summary_path)
