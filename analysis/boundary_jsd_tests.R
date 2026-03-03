## ============================================================
## Boundary JSD tests: consonants + tones
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
  if (!is.na(p)) {
    return(normalizePath(dirname(p)))
  }
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
  want_local <- tolower(Sys.getenv("BOUNDARY_LOAD_LOCAL", unset = "false")) %in% c("1", "true", "yes")

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

has_dots <- function(fn) {
  "..." %in% names(formals(fn))
}

has_arg <- function(fn, arg) {
  arg %in% names(formals(fn))
}

estimate_jsd_global_compat <- function(data,
                                       features,
                                       category_col,
                                       n_boot,
                                       min_tokens,
                                       bw,
                                       eval_on) {
  args <- list(
    data = data,
    features = features,
    category_col = category_col,
    do_boot = n_boot > 0L,
    n_boot = n_boot,
    min_tokens = min_tokens
  )

  if (has_dots(estimate_jsd) || has_arg(estimate_jsd, "bw")) {
    args$bw <- bw
  }
  if (has_dots(estimate_jsd) || has_arg(estimate_jsd, "eval_on")) {
    args$eval_on <- eval_on
  }

  do.call(estimate_jsd, args)
}

estimate_jsd_group_compat <- function(data,
                                      features,
                                      category_col,
                                      group_col,
                                      n_boot,
                                      min_tokens,
                                      bw,
                                      eval_on) {
  if (!group_col %in% names(data)) {
    stop("Group column not found: ", group_col)
  }

  jsd_args_base <- list(features = features, group = category_col)
  if (has_arg(jsd_kde_nd, "bw")) {
    jsd_args_base$bw <- bw
  }
  if (has_arg(jsd_kde_nd, "eval_on")) {
    jsd_args_base$eval_on <- eval_on
  }

  groups <- split(data, data[[group_col]])
  rows <- lapply(groups, function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens || dplyr::n_distinct(df_g[[category_col]]) != 2L) {
      return(NULL)
    }

    args_pt <- c(list(data = df_g), jsd_args_base)
    jsd_point <- tryCatch(
      do.call(jsd_kde_nd, args_pt),
      error = function(e) NA_real_
    )

    if (!is.finite(jsd_point)) {
      return(NULL)
    }

    if (n_boot > 0L) {
      vals <- replicate(n_boot, {
        samp <- df_g[sample.int(n_tok, n_tok, replace = TRUE), , drop = FALSE]
        if (dplyr::n_distinct(samp[[category_col]]) != 2L) {
          return(NA_real_)
        }
        args_bt <- c(list(data = samp), jsd_args_base)
        tryCatch(do.call(jsd_kde_nd, args_bt), error = function(e) NA_real_)
      })
      vals <- vals[is.finite(vals)]

      if (length(vals)) {
        jsd_mean <- mean(vals)
        jsd_sd <- stats::sd(vals)
        qs <- stats::quantile(vals, probs = c(0.025, 0.975), names = FALSE)
        jsd_low <- qs[1]
        jsd_high <- qs[2]
      } else {
        jsd_mean <- NA_real_
        jsd_sd <- NA_real_
        jsd_low <- NA_real_
        jsd_high <- NA_real_
      }
    } else {
      jsd_mean <- NA_real_
      jsd_sd <- NA_real_
      jsd_low <- NA_real_
      jsd_high <- NA_real_
    }

    tibble(
      group = as.character(df_g[[group_col]][1]),
      n_tokens = n_tok,
      n_boot = as.integer(n_boot),
      jsd_point = jsd_point,
      jsd_mean = jsd_mean,
      jsd_sd = jsd_sd,
      jsd_low = jsd_low,
      jsd_high = jsd_high
    )
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(tibble(
      group = character(0),
      n_tokens = integer(0),
      n_boot = integer(0),
      jsd_point = numeric(0),
      jsd_mean = numeric(0),
      jsd_sd = numeric(0),
      jsd_low = numeric(0),
      jsd_high = numeric(0)
    ))
  }

  bind_rows(rows)
}

parse_csv_string <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

parse_int_env <- function(name, default) {
  raw <- Sys.getenv(name, unset = as.character(default))
  val <- suppressWarnings(as.integer(raw))
  min_value <- if (grepl("^BOUNDARY_N_BOOT_", name)) 0L else 1L
  if (is.na(val) || val < min_value) {
    warning("Invalid value for ", name, "='", raw, "'. Using default ", default, ".")
    return(as.integer(default))
  }
  val
}

normalize_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x)] <- ""
  gsub("[-[:space:]_.]+", "", x)
}

read_input_data <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

default_feature_cols <- paste0("mfcc", 1:13)

data_path <- Sys.getenv("BOUNDARY_DATA_PATH", unset = "")
if (!nzchar(data_path)) {
  stop("Set BOUNDARY_DATA_PATH to a CSV/RDS with segment/tone labels and acoustic features.")
}

group_col <- Sys.getenv("BOUNDARY_GROUP_COL", unset = "speaker")
segment_col <- Sys.getenv("BOUNDARY_SEGMENT_COL", unset = "segment")
tone_col <- Sys.getenv("BOUNDARY_TONE_COL", unset = "tone")

feature_cols_env <- parse_csv_string(Sys.getenv("BOUNDARY_FEATURES", unset = ""))
n_boot_global <- parse_int_env("BOUNDARY_N_BOOT_GLOBAL", 300L)
n_boot_group <- parse_int_env("BOUNDARY_N_BOOT_GROUP", 150L)
min_tokens_global <- parse_int_env("BOUNDARY_MIN_TOKENS_GLOBAL", 80L)
min_tokens_group <- parse_int_env("BOUNDARY_MIN_TOKENS_GROUP", 30L)
min_per_category <- parse_int_env("BOUNDARY_MIN_PER_CATEGORY", 15L)

bw <- Sys.getenv("BOUNDARY_BW", unset = "Hpi.diag")
eval_on <- Sys.getenv("BOUNDARY_EVAL_ON", unset = "pooled")

out_dir <- Sys.getenv("BOUNDARY_OUT_DIR", unset = "")
if (!nzchar(out_dir)) {
  out_dir <- file.path(script_dir, "boundary_outputs")
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

raw_df <- read_input_data(data_path)

if (length(feature_cols_env)) {
  feature_cols <- feature_cols_env
} else if (all(default_feature_cols %in% names(raw_df))) {
  feature_cols <- default_feature_cols
} else {
  numeric_cols <- names(raw_df)[vapply(raw_df, is.numeric, logical(1))]
  auto_exclude <- unique(c(
    group_col, segment_col, tone_col,
    "start", "end", "t_start", "t_end"
  ))
  feature_cols <- setdiff(numeric_cols, auto_exclude)
  if (!length(feature_cols)) {
    stop("Could not determine feature columns. Set BOUNDARY_FEATURES explicitly.")
  }
  message("Auto-detected feature columns: ", paste(feature_cols, collapse = ", "))
}

missing_features <- setdiff(feature_cols, names(raw_df))
if (length(missing_features)) {
  stop("Missing feature columns: ", paste(missing_features, collapse = ", "))
}

contrast_specs <- list(
  list(
    contrast_id = "stop_voicing_b_p",
    domain = "consonant",
    source_type = "segment",
    cat1 = "b",
    cat2 = "p",
    cat1_alias = c("b", "bcl"),
    cat2_alias = c("p", "pcl")
  ),
  list(
    contrast_id = "stop_voicing_d_t",
    domain = "consonant",
    source_type = "segment",
    cat1 = "d",
    cat2 = "t",
    cat1_alias = c("d", "dcl"),
    cat2_alias = c("t", "tcl")
  ),
  list(
    contrast_id = "stop_voicing_g_k",
    domain = "consonant",
    source_type = "segment",
    cat1 = "g",
    cat2 = "k",
    cat1_alias = c("g", "gcl"),
    cat2_alias = c("k", "kcl")
  ),
  list(
    contrast_id = "fricative_s_sh",
    domain = "consonant",
    source_type = "segment",
    cat1 = "s",
    cat2 = "sh",
    cat1_alias = c("s"),
    cat2_alias = c("sh")
  ),
  list(
    contrast_id = "fricative_z_zh",
    domain = "consonant",
    source_type = "segment",
    cat1 = "z",
    cat2 = "zh",
    cat1_alias = c("z"),
    cat2_alias = c("zh")
  ),
  list(
    contrast_id = "fricative_f_th",
    domain = "consonant",
    source_type = "segment",
    cat1 = "f",
    cat2 = "th",
    cat1_alias = c("f"),
    cat2_alias = c("th")
  ),
  list(
    contrast_id = "tone_high_low",
    domain = "tone",
    source_type = "tone",
    cat1 = "high",
    cat2 = "low",
    cat1_alias = c("high", "h", "55", "44", "t1", "1"),
    cat2_alias = c("low", "l", "11", "22", "t3", "3")
  ),
  list(
    contrast_id = "tone_rising_falling",
    domain = "tone",
    source_type = "tone",
    cat1 = "rising",
    cat2 = "falling",
    cat1_alias = c("rising", "r", "35", "24", "t2", "2"),
    cat2_alias = c("falling", "f", "51", "53", "t4", "4")
  )
)

empty_global <- tibble(
  contrast_id = character(0),
  domain = character(0),
  source_col = character(0),
  cat1 = character(0),
  cat2 = character(0),
  n_tokens_total = integer(0),
  n_cat1 = integer(0),
  n_cat2 = integer(0),
  n_boot = integer(0),
  jsd_point = numeric(0),
  jsd_mean = numeric(0),
  jsd_sd = numeric(0),
  jsd_low = numeric(0),
  jsd_high = numeric(0),
  status = character(0),
  note = character(0)
)

empty_group <- tibble(
  contrast_id = character(0),
  domain = character(0),
  source_col = character(0),
  cat1 = character(0),
  cat2 = character(0),
  group = character(0),
  n_tokens = integer(0),
  n_boot = integer(0),
  jsd_point = numeric(0),
  jsd_mean = numeric(0),
  jsd_sd = numeric(0),
  jsd_low = numeric(0),
  jsd_high = numeric(0),
  status = character(0),
  note = character(0)
)

empty_skipped <- tibble(
  contrast_id = character(0),
  domain = character(0),
  source_col = character(0),
  cat1 = character(0),
  cat2 = character(0),
  status = character(0),
  note = character(0)
)

global_rows <- list()
group_rows <- list()
skipped_rows <- list()

for (spec in contrast_specs) {
  source_col <- if (spec$source_type == "segment") segment_col else tone_col

  if (!source_col %in% names(raw_df)) {
    skipped_rows[[length(skipped_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      source_col = source_col,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = paste0("Missing source column '", source_col, "'.")
    )
    next
  }

  keep_cols <- unique(c(source_col, group_col, feature_cols))
  keep_cols <- keep_cols[keep_cols %in% names(raw_df)]
  df <- raw_df[, keep_cols, drop = FALSE]
  df <- df[stats::complete.cases(df[, c(source_col, feature_cols), drop = FALSE]), , drop = FALSE]

  lbl <- normalize_label(df[[source_col]])
  is_cat1 <- lbl %in% normalize_label(spec$cat1_alias)
  is_cat2 <- lbl %in% normalize_label(spec$cat2_alias)
  keep <- xor(is_cat1, is_cat2)

  if (!any(keep)) {
    skipped_rows[[length(skipped_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      source_col = source_col,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = "No matching tokens for aliases."
    )
    next
  }

  df <- df[keep, , drop = FALSE]
  df$category <- ifelse(is_cat1[keep], spec$cat1, spec$cat2)

  n_cat1 <- sum(df$category == spec$cat1)
  n_cat2 <- sum(df$category == spec$cat2)
  n_total <- nrow(df)

  if (n_cat1 < min_per_category || n_cat2 < min_per_category) {
    skipped_rows[[length(skipped_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      source_col = source_col,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = paste0(
        "Too few tokens per category after filtering: ",
        spec$cat1, "=", n_cat1, ", ",
        spec$cat2, "=", n_cat2, "."
      )
    )
    next
  }

  global_err <- NULL
  global_fit <- tryCatch(
    estimate_jsd_global_compat(
      data = df,
      features = feature_cols,
      category_col = "category",
      n_boot = n_boot_global,
      min_tokens = min_tokens_global,
      bw = bw,
      eval_on = eval_on
    ),
    error = function(e) {
      global_err <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(global_fit)) {
    global_rows[[length(global_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      source_col = source_col,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      n_tokens_total = n_total,
      n_cat1 = n_cat1,
      n_cat2 = n_cat2,
      n_boot = 0L,
      jsd_point = NA_real_,
      jsd_mean = NA_real_,
      jsd_sd = NA_real_,
      jsd_low = NA_real_,
      jsd_high = NA_real_,
      status = "error",
      note = global_err
    )
  } else {
    global_rows[[length(global_rows) + 1L]] <- global_fit %>%
      mutate(
        contrast_id = spec$contrast_id,
        domain = spec$domain,
        source_col = source_col,
        cat1 = spec$cat1,
        cat2 = spec$cat2,
        n_tokens_total = n_total,
        n_cat1 = n_cat1,
        n_cat2 = n_cat2,
        status = "ok",
        note = NA_character_
      ) %>%
      select(
        contrast_id, domain, source_col, cat1, cat2,
        n_tokens_total, n_cat1, n_cat2, n_boot,
        jsd_point, jsd_mean, jsd_sd, jsd_low, jsd_high,
        status, note
      )
  }

  if (!group_col %in% names(df)) {
    next
  }

  group_err <- NULL
  group_fit <- tryCatch(
    estimate_jsd_group_compat(
      data = df,
      features = feature_cols,
      category_col = "category",
      group_col = group_col,
      n_boot = n_boot_group,
      min_tokens = min_tokens_group,
      bw = bw,
      eval_on = eval_on
    ),
    error = function(e) {
      group_err <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(group_fit)) {
    group_rows[[length(group_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      source_col = source_col,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      group = NA_character_,
      n_tokens = NA_integer_,
      n_boot = 0L,
      jsd_point = NA_real_,
      jsd_mean = NA_real_,
      jsd_sd = NA_real_,
      jsd_low = NA_real_,
      jsd_high = NA_real_,
      status = "error",
      note = group_err
    )
  } else if (!nrow(group_fit)) {
    group_rows[[length(group_rows) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      source_col = source_col,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      group = NA_character_,
      n_tokens = NA_integer_,
      n_boot = 0L,
      jsd_point = NA_real_,
      jsd_mean = NA_real_,
      jsd_sd = NA_real_,
      jsd_low = NA_real_,
      jsd_high = NA_real_,
      status = "empty",
      note = paste0("No groups met min_tokens=", min_tokens_group, ".")
    )
  } else {
    group_rows[[length(group_rows) + 1L]] <- group_fit %>%
      mutate(
        contrast_id = spec$contrast_id,
        domain = spec$domain,
        source_col = source_col,
        cat1 = spec$cat1,
        cat2 = spec$cat2,
        status = "ok",
        note = NA_character_,
        group = as.character(group)
      ) %>%
      select(
        contrast_id, domain, source_col, cat1, cat2,
        group, n_tokens, n_boot,
        jsd_point, jsd_mean, jsd_sd, jsd_low, jsd_high,
        status, note
      )
  }
}

global_out <- if (length(global_rows)) bind_rows(global_rows) else empty_global
group_out <- if (length(group_rows)) bind_rows(group_rows) else empty_group
skipped_out <- if (length(skipped_rows)) bind_rows(skipped_rows) else empty_skipped

global_path <- file.path(out_dir, "boundary_jsd_global.csv")
group_path <- file.path(out_dir, "boundary_jsd_by_group.csv")
skipped_path <- file.path(out_dir, "boundary_jsd_skipped.csv")

write.csv(global_out, global_path, row.names = FALSE)
write.csv(group_out, group_path, row.names = FALSE)
write.csv(skipped_out, skipped_path, row.names = FALSE)

message("Finished boundary JSD tests.")
message("Global results: ", global_path)
message("By-group results: ", group_path)
message("Skipped contrasts: ", skipped_path)
