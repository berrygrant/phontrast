## ============================================================
## Validate JSD on classical consonant metrics
## - stop voicing via VOT
## - fricative place via spectral moments + common cues
## ============================================================

set.seed(20260224)

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
  library(tidyr)
  library(tibble)
})

load_phonjsd <- function(script_dir) {
  pkg_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  use_local <- FALSE
  want_local <- tolower(Sys.getenv("CLASSICAL_LOAD_LOCAL", unset = "false")) %in% c("1", "true", "yes")
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

has_arg <- function(fn, arg) {
  arg %in% names(formals(fn))
}

parse_int_env <- function(name, default, min_value = 0L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  val <- suppressWarnings(as.integer(raw))
  if (is.na(val) || val < min_value) {
    warning("Invalid ", name, "='", raw, "'. Using ", default, ".")
    return(as.integer(default))
  }
  val
}

safe_read <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

cohen_d <- function(x1, x2) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  n1 <- length(x1)
  n2 <- length(x2)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s1 <- stats::sd(x1)
  s2 <- stats::sd(x2)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x2) - mean(x1)) / sp
}

jsd_point_compat <- function(data, features, category_col, bw, eval_on) {
  if (length(features) == 1L) {
    feature <- features[1]
    levs <- unique(data[[category_col]])
    if (length(levs) != 2L) stop("Need exactly two categories for 1D JSD.")
    x1 <- data[data[[category_col]] == levs[1], feature, drop = TRUE]
    x2 <- data[data[[category_col]] == levs[2], feature, drop = TRUE]
    x1 <- x1[is.finite(x1)]
    x2 <- x2[is.finite(x2)]
    if (length(x1) < 2L || length(x2) < 2L) stop("Not enough finite values for 1D JSD.")

    rng <- range(c(x1, x2))
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
      return(0)
    }

    n_grid <- 512L
    d1 <- stats::density(x1, n = n_grid, from = rng[1], to = rng[2])
    d2 <- stats::density(x2, n = n_grid, from = rng[1], to = rng[2])
    p <- pmax(d1$y, 0)
    q <- pmax(d2$y, 0)
    p <- p / sum(p)
    q <- q / sum(q)
    return(jsd(p, q))
  }

  args <- list(
    data = data,
    features = features,
    group = category_col
  )
  if (has_arg(jsd_kde_nd, "bw")) args$bw <- bw
  if (has_arg(jsd_kde_nd, "eval_on")) args$eval_on <- eval_on
  do.call(jsd_kde_nd, args)
}

jsd_summary_boot <- function(data,
                             features,
                             category_col,
                             n_boot = 0L,
                             bw = "Hpi.diag",
                             eval_on = "pooled") {
  pt <- jsd_point_compat(data, features, category_col, bw, eval_on)
  if (n_boot <= 0L) {
    return(list(
      jsd_point = pt,
      jsd_mean = NA_real_,
      jsd_sd = NA_real_,
      jsd_low = NA_real_,
      jsd_high = NA_real_,
      n_boot = 0L
    ))
  }

  n <- nrow(data)
  vals <- replicate(n_boot, {
    idx <- sample.int(n, n, replace = TRUE)
    d <- data[idx, , drop = FALSE]
    if (dplyr::n_distinct(d[[category_col]]) != 2L) return(NA_real_)
    tryCatch(
      jsd_point_compat(d, features, category_col, bw, eval_on),
      error = function(e) NA_real_
    )
  })
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(list(
      jsd_point = pt,
      jsd_mean = NA_real_,
      jsd_sd = NA_real_,
      jsd_low = NA_real_,
      jsd_high = NA_real_,
      n_boot = 0L
    ))
  }
  qs <- stats::quantile(vals, c(0.025, 0.975), names = FALSE)
  list(
    jsd_point = pt,
    jsd_mean = mean(vals),
    jsd_sd = stats::sd(vals),
    jsd_low = qs[1],
    jsd_high = qs[2],
    n_boot = as.integer(length(vals))
  )
}

data_path <- Sys.getenv(
  "CLASSICAL_DATA_PATH",
  unset = file.path(script_dir, "sbcae_classical_consonant_metrics.csv")
)
out_dir <- Sys.getenv(
  "CLASSICAL_OUT_DIR",
  unset = file.path(script_dir, "classical_outputs")
)
bw <- Sys.getenv("CLASSICAL_BW", unset = "Hpi.diag")
eval_on <- Sys.getenv("CLASSICAL_EVAL_ON", unset = "pooled")
n_boot <- parse_int_env("CLASSICAL_N_BOOT", 0L, min_value = 0L)
min_per_category <- parse_int_env("CLASSICAL_MIN_PER_CATEGORY", 30L, min_value = 2L)
max_per_category <- parse_int_env("CLASSICAL_MAX_PER_CATEGORY", 800L, min_value = 10L)
low_jsd_threshold <- suppressWarnings(as.numeric(Sys.getenv("CLASSICAL_LOW_JSD_THRESHOLD", unset = "0.20")))
if (!is.finite(low_jsd_threshold) || low_jsd_threshold <= 0) low_jsd_threshold <- 0.20

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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
  list(
    contrast_id = "stop_voicing_b_p",
    domain = "stop_voicing",
    cat1 = "b",
    cat2 = "p",
    features = c("vot_ms"),
    expected_jsd = "low"
  ),
  list(
    contrast_id = "stop_voicing_d_t",
    domain = "stop_voicing",
    cat1 = "d",
    cat2 = "t",
    features = c("vot_ms"),
    expected_jsd = "low"
  ),
  list(
    contrast_id = "stop_voicing_g_k",
    domain = "stop_voicing",
    cat1 = "g",
    cat2 = "k",
    features = c("vot_ms"),
    expected_jsd = "low"
  ),
  list(
    contrast_id = "fricative_place_s_sh",
    domain = "fricative_place",
    cat1 = "s",
    cat2 = "sh",
    features = fric_features,
    expected_jsd = "high"
  ),
  list(
    contrast_id = "fricative_place_f_th",
    domain = "fricative_place",
    cat1 = "f",
    cat2 = "th",
    features = fric_features,
    expected_jsd = "high"
  ),
  list(
    contrast_id = "fricative_place_z_zh",
    domain = "fricative_place",
    cat1 = "z",
    cat2 = "zh",
    features = fric_features,
    expected_jsd = "high"
  )
)

contrast_rows <- list()
feature_rows <- list()
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
    select(all_of(c("speaker", "category", feat))) %>%
    filter(stats::complete.cases(.))

  if (max_per_category > 0L) {
    d <- d %>%
      group_by(.data$category) %>%
      group_modify(~ {
        if (nrow(.x) > max_per_category) {
          slice_sample(.x, n = max_per_category)
        } else {
          .x
        }
      }) %>%
      ungroup()
  }

  n1 <- sum(d$category == spec$cat1)
  n2 <- sum(d$category == spec$cat2)

  if (n1 < min_per_category || n2 < min_per_category) {
    skipped_rows[[length(skipped_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = paste0(
        "Too few tokens after filtering: ",
        spec$cat1, "=", n1, ", ",
        spec$cat2, "=", n2
      )
    )
    next
  }

  jsd_res <- tryCatch(
    jsd_summary_boot(
      data = d,
      features = feat,
      category_col = "category",
      n_boot = n_boot,
      bw = bw,
      eval_on = eval_on
    ),
    error = function(e) e
  )

  if (inherits(jsd_res, "error")) {
    contrast_rows[[length(contrast_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      feature_set = paste(feat, collapse = ","),
      n_tokens_total = nrow(d),
      n_cat1 = n1,
      n_cat2 = n2,
      n_boot = 0L,
      jsd_point = NA_real_,
      jsd_mean = NA_real_,
      jsd_sd = NA_real_,
      jsd_low = NA_real_,
      jsd_high = NA_real_,
      expected_jsd = spec$expected_jsd,
      low_jsd_threshold = low_jsd_threshold,
      jsd_is_low = NA,
      pillai = NA_real_,
      status = "error",
      note = conditionMessage(jsd_res)
    )
    next
  }

  pillai_val <- NA_real_
  if (length(feat) > 1L) {
    pillai_val <- tryCatch(
      estimate_pillai(
        data = d,
        features = feat,
        category_col = "category"
      )$pillai[1],
      error = function(e) NA_real_
    )
  }

  contrast_rows[[length(contrast_rows) + 1L]] <- tibble(
    contrast_id = spec$contrast_id,
    domain = spec$domain,
    cat1 = spec$cat1,
    cat2 = spec$cat2,
    feature_set = paste(feat, collapse = ","),
    n_tokens_total = nrow(d),
    n_cat1 = n1,
    n_cat2 = n2,
    n_boot = as.integer(jsd_res$n_boot),
    jsd_point = jsd_res$jsd_point,
    jsd_mean = jsd_res$jsd_mean,
    jsd_sd = jsd_res$jsd_sd,
    jsd_low = jsd_res$jsd_low,
    jsd_high = jsd_res$jsd_high,
    expected_jsd = spec$expected_jsd,
    low_jsd_threshold = low_jsd_threshold,
    jsd_is_low = isTRUE(jsd_res$jsd_point <= low_jsd_threshold),
    pillai = pillai_val,
    status = "ok",
    note = NA_character_
  )

  for (f in feat) {
    x1 <- d %>% filter(.data$category == spec$cat1) %>% pull(.data[[f]])
    x2 <- d %>% filter(.data$category == spec$cat2) %>% pull(.data[[f]])
    feature_rows[[length(feature_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      feature = f,
      mean_cat1 = mean(x1),
      mean_cat2 = mean(x2),
      sd_cat1 = stats::sd(x1),
      sd_cat2 = stats::sd(x2),
      mean_diff_cat2_minus_cat1 = mean(x2) - mean(x1),
      cohen_d_cat2_vs_cat1 = cohen_d(x1, x2),
      n_cat1 = length(x1),
      n_cat2 = length(x2)
    )
  }
}

contrast_df <- if (length(contrast_rows)) bind_rows(contrast_rows) else tibble()
feature_df <- if (length(feature_rows)) bind_rows(feature_rows) else tibble()
skipped_df <- if (length(skipped_rows)) bind_rows(skipped_rows) else tibble()

contrast_path <- file.path(out_dir, "classical_contrast_jsd_summary.csv")
feature_path <- file.path(out_dir, "classical_feature_descriptives.csv")
skipped_path <- file.path(out_dir, "classical_skipped_contrasts.csv")

write.csv(contrast_df, contrast_path, row.names = FALSE)
write.csv(feature_df, feature_path, row.names = FALSE)
write.csv(skipped_df, skipped_path, row.names = FALSE)

message("Finished classical boundary validation.")
message("Using local source: ", using_local_source)
message("Contrast summary: ", contrast_path)
message("Feature descriptives: ", feature_path)
message("Skipped contrasts: ", skipped_path)
