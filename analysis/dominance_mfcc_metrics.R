## ============================================================
## Dominance analysis for MFCC vowel-pair metrics
## ============================================================

set.seed(20260303)

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

safe_read <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

scale_numeric <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

all_subsets <- function(x) {
  x <- as.character(x)
  out <- list(character(0))
  if (!length(x)) return(out)
  for (k in seq_along(x)) {
    cmb <- combn(x, k, simplify = FALSE)
    out <- c(out, cmb)
  }
  out
}

fit_lm_stats <- function(df, outcome, predictors) {
  if (length(predictors) == 0) {
    form <- stats::as.formula(paste(outcome, "~ 1"))
  } else {
    form <- stats::as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  }
  fit <- stats::lm(form, data = df)
  smry <- summary(fit)
  list(
    fit = fit,
    r2 = unname(smry$r.squared),
    adj_r2 = unname(smry$adj.r.squared),
    rmse = sqrt(mean(stats::residuals(fit)^2))
  )
}

pairwise_complete_dominance <- function(detail, predictors) {
  rows <- list()
  idx <- 1L
  for (a in predictors) {
    for (b in predictors) {
      if (identical(a, b)) next
      sub_a <- detail[detail$predictor == a, c("subset_key", "subset_size", "delta_r2")]
      sub_b <- detail[detail$predictor == b, c("subset_key", "subset_size", "delta_r2")]
      names(sub_a)[names(sub_a) == "delta_r2"] <- "delta_a"
      names(sub_b)[names(sub_b) == "delta_r2"] <- "delta_b"
      merged <- merge(sub_a, sub_b, by = c("subset_key", "subset_size"), all = FALSE)
      if (!nrow(merged)) next
      rows[[idx]] <- tibble(
        predictor = a,
        comparator = b,
        complete_dominates = all(merged$delta_a >= merged$delta_b),
        min_delta_margin = min(merged$delta_a - merged$delta_b),
        mean_delta_margin = mean(merged$delta_a - merged$delta_b)
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(rows)
}

dominance_analysis <- function(df, outcome, predictors, analysis_name, base_predictors = character(0)) {
  full_terms <- c(base_predictors, predictors)
  base_fit <- fit_lm_stats(df, outcome, base_predictors)
  full_stats <- fit_lm_stats(df, outcome, full_terms)
  added_r2 <- full_stats$r2 - base_fit$r2

  detail_rows <- list()
  conditional_rows <- list()
  summary_rows <- list()
  detail_idx <- 1L
  cond_idx <- 1L
  sum_idx <- 1L

  for (pred in predictors) {
    others <- setdiff(predictors, pred)
    subsets <- all_subsets(others)
    deltas <- numeric(length(subsets))
    subset_sizes <- integer(length(subsets))
    subset_keys <- character(length(subsets))

    for (i in seq_along(subsets)) {
      subset_terms <- subsets[[i]]
      without_terms <- c(base_predictors, subset_terms)
      with_terms <- c(base_predictors, subset_terms, pred)
      r2_without <- fit_lm_stats(df, outcome, without_terms)$r2
      r2_with <- fit_lm_stats(df, outcome, with_terms)$r2
      delta <- r2_with - r2_without
      subset_sizes[i] <- length(subset_terms)
      subset_keys[i] <- if (length(subset_terms)) paste(sort(subset_terms), collapse = "+") else "(none)"
      deltas[i] <- delta

      detail_rows[[detail_idx]] <- tibble(
        analysis = analysis_name,
        predictor = pred,
        subset_key = subset_keys[i],
        subset_size = subset_sizes[i],
        delta_r2 = delta
      )
      detail_idx <- detail_idx + 1L
    }

    conditional_means <- numeric(0)
    for (k in sort(unique(subset_sizes))) {
      cond_val <- mean(deltas[subset_sizes == k])
      conditional_means <- c(conditional_means, cond_val)
      conditional_rows[[cond_idx]] <- tibble(
        analysis = analysis_name,
        predictor = pred,
        model_size_without_predictor = k,
        conditional_dominance = cond_val
      )
      cond_idx <- cond_idx + 1L
    }

    general <- mean(conditional_means)
    summary_rows[[sum_idx]] <- tibble(
      analysis = analysis_name,
      predictor = pred,
      general_dominance = general,
      share_of_added_r2 = if (isTRUE(all.equal(added_r2, 0))) NA_real_ else general / added_r2,
      base_model_r2 = base_fit$r2,
      full_model_r2 = full_stats$r2,
      added_r2_over_base = added_r2,
      full_model_adj_r2 = full_stats$adj_r2,
      full_model_rmse = full_stats$rmse
    )
    sum_idx <- sum_idx + 1L
  }

  detail_df <- dplyr::bind_rows(detail_rows)
  conditional_df <- dplyr::bind_rows(conditional_rows)
  summary_df <- dplyr::bind_rows(summary_rows) %>% arrange(desc(general_dominance))

  complete_df <- pairwise_complete_dominance(detail_df, predictors) %>%
    mutate(analysis = analysis_name, .before = 1L)

  coef_df <- summary(full_stats$fit)$coefficients
  coef_tbl <- tibble(
    analysis = analysis_name,
    term = rownames(coef_df),
    estimate = coef_df[, "Estimate"],
    std_error = coef_df[, "Std. Error"],
    t_value = coef_df[, "t value"],
    p_value = coef_df[, "Pr(>|t|)"]
  )

  model_tbl <- tibble(
    analysis = analysis_name,
    outcome = outcome,
    n_rows = nrow(df),
    base_model_r2 = base_fit$r2,
    full_model_r2 = full_stats$r2,
    full_model_adj_r2 = full_stats$adj_r2,
    full_model_rmse = full_stats$rmse,
    added_r2_over_base = full_stats$r2 - base_fit$r2
  )

  list(
    summary = summary_df,
    conditional = conditional_df,
    detail = detail_df,
    complete = complete_df,
    coefficients = coef_tbl,
    model = model_tbl
  )
}

metric_path <- Sys.getenv(
  "MFCC_METRIC_PATH",
  unset = file.path(script_dir, "..", "comparison_outputs", "pairwise_metrics_mfcc13.csv")
)
out_dir <- Sys.getenv(
  "MFCC_DOMINANCE_OUT_DIR",
  unset = file.path(script_dir, "..", "comparison_outputs")
)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

raw_df <- safe_read(metric_path)
required_cols <- c("n_tokens", "jsd", "pillai", "bhatt", "overlap")
missing_cols <- setdiff(required_cols, names(raw_df))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- raw_df %>%
  mutate(
    n_tokens = as.numeric(n_tokens),
    jsd = as.numeric(jsd),
    pillai = as.numeric(pillai),
    bhatt = as.numeric(bhatt),
    overlap = as.numeric(overlap),
    sep_from_overlap = 1 - overlap
  ) %>%
  filter(stats::complete.cases(n_tokens, jsd, pillai, bhatt, overlap, sep_from_overlap)) %>%
  mutate(
    jsd_z = scale_numeric(jsd),
    pillai_z = scale_numeric(pillai),
    bhatt_z = scale_numeric(bhatt),
    n_tokens_z = scale_numeric(n_tokens)
  )

predictors <- c("jsd_z", "pillai_z", "bhatt_z")

unadjusted <- dominance_analysis(
  df = df,
  outcome = "sep_from_overlap",
  predictors = predictors,
  analysis_name = "mfcc_unadjusted"
)

adjusted <- dominance_analysis(
  df = df,
  outcome = "sep_from_overlap",
  predictors = predictors,
  analysis_name = "mfcc_adjusted_n_tokens",
  base_predictors = c("n_tokens_z")
)

summary_df <- bind_rows(unadjusted$summary, adjusted$summary) %>%
  mutate(
    predictor = dplyr::recode(
      predictor,
      jsd_z = "JSD",
      pillai_z = "Pillai",
      bhatt_z = "Bhattacharyya"
    )
  )

conditional_df <- bind_rows(unadjusted$conditional, adjusted$conditional) %>%
  mutate(
    predictor = dplyr::recode(
      predictor,
      jsd_z = "JSD",
      pillai_z = "Pillai",
      bhatt_z = "Bhattacharyya"
    )
  )

complete_df <- bind_rows(unadjusted$complete, adjusted$complete) %>%
  mutate(
    predictor = dplyr::recode(predictor, jsd_z = "JSD", pillai_z = "Pillai", bhatt_z = "Bhattacharyya"),
    comparator = dplyr::recode(comparator, jsd_z = "JSD", pillai_z = "Pillai", bhatt_z = "Bhattacharyya")
  )

coef_df <- bind_rows(unadjusted$coefficients, adjusted$coefficients) %>%
  mutate(
    term = dplyr::recode(
      term,
      `(Intercept)` = "Intercept",
      jsd_z = "JSD",
      pillai_z = "Pillai",
      bhatt_z = "Bhattacharyya",
      n_tokens_z = "n_tokens"
    )
  )

model_df <- bind_rows(unadjusted$model, adjusted$model)

detail_df <- bind_rows(unadjusted$detail, adjusted$detail) %>%
  mutate(
    predictor = dplyr::recode(predictor, jsd_z = "JSD", pillai_z = "Pillai", bhatt_z = "Bhattacharyya")
  )

readable_lines <- c(
  sprintf("Rows analyzed: %d", nrow(df)),
  "",
  "General dominance (share of explained R^2):"
)

for (analysis_name in unique(summary_df$analysis)) {
  readable_lines <- c(readable_lines, paste0("- ", analysis_name))
  block <- summary_df %>% filter(analysis == analysis_name) %>% arrange(desc(general_dominance))
  for (i in seq_len(nrow(block))) {
    readable_lines <- c(
      readable_lines,
      sprintf(
        "  %s: dominance=%.6f, share=%.4f, added_R2=%.4f, full_R2=%.4f, adj_R2=%.4f",
        block$predictor[i],
        block$general_dominance[i],
        block$share_of_added_r2[i],
        block$added_r2_over_base[i],
        block$full_model_r2[i],
        block$full_model_adj_r2[i]
      )
    )
  }
  cmp <- complete_df %>% filter(analysis == analysis_name, complete_dominates)
  if (nrow(cmp)) {
    readable_lines <- c(readable_lines, "  Complete dominance:")
    for (i in seq_len(nrow(cmp))) {
      readable_lines <- c(
        readable_lines,
        sprintf(
          "    %s > %s (min margin %.6f, mean margin %.6f)",
          cmp$predictor[i], cmp$comparator[i], cmp$min_delta_margin[i], cmp$mean_delta_margin[i]
        )
      )
    }
  }
}

write.csv(summary_df, file.path(out_dir, "mfcc_metric_dominance_summary.csv"), row.names = FALSE)
write.csv(conditional_df, file.path(out_dir, "mfcc_metric_dominance_conditional.csv"), row.names = FALSE)
write.csv(complete_df, file.path(out_dir, "mfcc_metric_dominance_complete.csv"), row.names = FALSE)
write.csv(coef_df, file.path(out_dir, "mfcc_metric_dominance_coefficients.csv"), row.names = FALSE)
write.csv(model_df, file.path(out_dir, "mfcc_metric_dominance_models.csv"), row.names = FALSE)
write.csv(detail_df, file.path(out_dir, "mfcc_metric_dominance_detail.csv"), row.names = FALSE)
writeLines(readable_lines, file.path(out_dir, "mfcc_metric_dominance_summary.txt"))

message("Wrote MFCC dominance outputs to: ", out_dir)
