#!/usr/bin/env Rscript

# PB52 MMO-framework analysis for /E/ ~ /I/ in normalized F1/F2.
#
# This is the acoustically faithful version for the poster: Peterson & Barney
# has F1/F2 and repeated productions by speaker. It is not the full Smith et al.
# covariate design because PB52 does not provide lexical frequency or
# phonological context; vowel identity is effectively tied to the elicited word.

suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(phonTools)
  library(posterior)
})

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE)
}

analysis_dir <- dirname(script_path())
output_dir <- file.path(analysis_dir, "data")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    fit = FALSE,
    chains = 4L,
    iter = 2000L,
    warmup = 1000L,
    posterior_draws = 1000L,
    seed = 20260625L,
    backend = "cmdstanr",
    run_label = "",
    adapt_delta = 0.95,
    max_treedepth = 12L,
    population_only = FALSE
  )

  for (arg in args) {
    if (arg == "--fit") {
      out$fit <- TRUE
    } else if (arg == "--population-only") {
      out$population_only <- TRUE
    } else if (grepl("^--chains=", arg)) {
      out$chains <- as.integer(sub("^--chains=", "", arg))
    } else if (grepl("^--iter=", arg)) {
      out$iter <- as.integer(sub("^--iter=", "", arg))
    } else if (grepl("^--warmup=", arg)) {
      out$warmup <- as.integer(sub("^--warmup=", "", arg))
    } else if (grepl("^--posterior-draws=", arg)) {
      out$posterior_draws <- as.integer(sub("^--posterior-draws=", "", arg))
    } else if (grepl("^--seed=", arg)) {
      out$seed <- as.integer(sub("^--seed=", "", arg))
    } else if (grepl("^--backend=", arg)) {
      out$backend <- sub("^--backend=", "", arg)
    } else if (grepl("^--run-label=", arg)) {
      out$run_label <- sub("^--run-label=", "", arg)
    } else if (grepl("^--adapt-delta=", arg)) {
      out$adapt_delta <- as.numeric(sub("^--adapt-delta=", "", arg))
    } else if (grepl("^--max-treedepth=", arg)) {
      out$max_treedepth <- as.integer(sub("^--max-treedepth=", "", arg))
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  out
}

safe_slug <- function(x) {
  x <- gsub("[^[:alnum:]_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

analysis_slug <- function(run_label) {
  slug <- "pb52_E_I_mmo"
  if (!is.null(run_label) && nzchar(run_label)) {
    slug <- paste(slug, safe_slug(run_label), sep = "_")
  }
  slug
}

lobanov <- function(x) {
  as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

prepare_data <- function() {
  data(pb52, package = "phonTools")
  all_vowels <- pb52
  all_vowels$speaker <- factor(all_vowels$speaker)

  # Lobanov-style speaker normalization over the full vowel space, then filter
  # to the target contrast. This keeps /E/ and /I/ on a speaker-normalized F1/F2
  # scale without estimating normalization from the contrast alone.
  all_vowels$f1_z <- ave(all_vowels$f1, all_vowels$speaker, FUN = lobanov)
  all_vowels$f2_z <- ave(all_vowels$f2, all_vowels$speaker, FUN = lobanov)

  df <- all_vowels[all_vowels$vowel %in% c("E", "I"), , drop = FALSE]
  df <- df[complete.cases(df[, c("type", "speaker", "vowel", "repetition", "f1_z", "f2_z")]), , drop = FALSE]
  df$vowel <- factor(df$vowel, levels = c("E", "I"))
  df$vowel_ipa <- ifelse(df$vowel == "E", "/ɛ/", "/ɪ/")
  df$type <- factor(df$type)
  df$speaker <- factor(df$speaker)
  df$repetition <- factor(df$repetition)
  df
}

build_formula <- function() {
  bf_f1 <- bf(f1_z ~ vowel + type + repetition + (1 + vowel | speaker))
  bf_f2 <- bf(f2_z ~ vowel + type + repetition + (1 + vowel | speaker))
  bf_f1 + bf_f2 + set_rescor(TRUE)
}

build_priors <- function() {
  c(
    prior(normal(0, 2), class = "b", resp = "f1z"),
    prior(normal(0, 2), class = "b", resp = "f2z"),
    prior(student_t(3, 0, 2.5), class = "Intercept", resp = "f1z"),
    prior(student_t(3, 0, 2.5), class = "Intercept", resp = "f2z"),
    prior(exponential(1), class = "sd", resp = "f1z"),
    prior(exponential(1), class = "sd", resp = "f2z"),
    prior(exponential(1), class = "sigma", resp = "f1z"),
    prior(exponential(1), class = "sigma", resp = "f2z"),
    prior(lkj(2), class = "rescor")
  )
}

bhatt_affinity <- function(x, y, ridge = 1e-6) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  mu_x <- colMeans(x)
  mu_y <- colMeans(y)
  cov_x <- cov(x) + diag(ridge, ncol(x))
  cov_y <- cov(y) + diag(ridge, ncol(y))
  cov_mean <- (cov_x + cov_y) / 2
  diff <- matrix(mu_x - mu_y, ncol = 1)
  term1 <- as.numeric(t(diff) %*% solve(cov_mean, diff)) / 8
  det_x <- determinant(cov_x, logarithm = TRUE)$modulus[1]
  det_y <- determinant(cov_y, logarithm = TRUE)$modulus[1]
  det_mean <- determinant(cov_mean, logarithm = TRUE)$modulus[1]
  term2 <- 0.5 * (det_mean - 0.5 * (det_x + det_y))
  exp(-(term1 + term2))
}

posterior_array <- function(pp) {
  if (is.list(pp)) {
    draws <- nrow(pp[[1]])
    obs <- ncol(pp[[1]])
    resp <- length(pp)
    arr <- array(NA_real_, dim = c(draws, obs, resp))
    for (i in seq_along(pp)) arr[, , i] <- pp[[i]]
    return(arr)
  }
  if (length(dim(pp)) == 3L) return(pp)
  stop("Unexpected posterior_predict shape.", call. = FALSE)
}

write_preflight <- function(df, opts, slug) {
  model_data_path <- file.path(output_dir, paste0(slug, "_model_data.csv"))
  summary_path <- file.path(output_dir, paste0(slug, "_preflight_summary.csv"))
  fig_path <- file.path(output_dir, paste0(slug, "_f1f2_tokens.png"))

  write.csv(df, model_data_path, row.names = FALSE)

  summary_df <- data.frame(
    contrast = "/ɛ/ ~ /ɪ/",
    data = "Peterson & Barney 1952 via phonTools::pb52",
    response_space = "speaker-normalized F1/F2",
    fit_requested = opts$fit,
    n_model_rows = nrow(df),
    n_E = sum(df$vowel == "E"),
    n_I = sum(df$vowel == "I"),
    n_speakers = length(unique(df$speaker)),
    n_talker_types = length(unique(df$type)),
    formula = "f1_z/f2_z ~ vowel + type + repetition + (1 + vowel | speaker)",
    limitation = "PB52 supports F1/F2 and speaker controls, but not lexical-frequency or phonological-context controls."
  )
  write.csv(summary_df, summary_path, row.names = FALSE)

  png(fig_path, width = 1600, height = 1200, res = 220)
  print(
    ggplot(df, aes(x = f2_z, y = f1_z, color = vowel_ipa)) +
      geom_point(alpha = 0.68, size = 1.9) +
      stat_ellipse(level = 0.68, linewidth = 0.9) +
      scale_x_reverse() +
      scale_y_reverse() +
      scale_color_manual(values = c("/ɛ/" = "#D7263D", "/ɪ/" = "#15928A")) +
      labs(
        x = "F2, speaker-normalized",
        y = "F1, speaker-normalized",
        color = "Vowel"
      ) +
      theme_minimal(base_size = 16)
  )
  dev.off()

  cat("Wrote", model_data_path, "\n")
  cat("Wrote", summary_path, "\n")
  cat("Wrote", fig_path, "\n")
  print(summary_df)
}

write_diagnostics <- function(fit, slug, opts) {
  diag_path <- file.path(output_dir, paste0(slug, "_diagnostics.csv"))
  draws_summary <- posterior::summarise_draws(as_draws_df(fit))
  np <- tryCatch(nuts_params(fit), error = function(e) NULL)

  if (!is.null(np)) {
    treedepth <- subset(np, Parameter == "treedepth__")$Value
    divergent <- subset(np, Parameter == "divergent__")$Value
    max_td <- max(treedepth)
    td_hits <- sum(treedepth >= opts$max_treedepth)
    divergences <- sum(divergent)
  } else {
    max_td <- NA_real_
    td_hits <- NA_integer_
    divergences <- NA_integer_
  }

  diag_df <- data.frame(
    max_rhat = max(draws_summary$rhat, na.rm = TRUE),
    median_rhat = median(draws_summary$rhat, na.rm = TRUE),
    min_bulk_ess = min(draws_summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draws_summary$ess_tail, na.rm = TRUE),
    max_treedepth_limit = opts$max_treedepth,
    max_treedepth_observed = max_td,
    treedepth_limit_hit_count = td_hits,
    divergences = divergences
  )
  write.csv(diag_df, diag_path, row.names = FALSE)
  cat("Wrote", diag_path, "\n")
  print(diag_df)
  invisible(diag_df)
}

run_fit <- function(df, opts, slug) {
  detected_cores <- parallel::detectCores(logical = TRUE)
  if (length(detected_cores) != 1L || is.na(detected_cores) || detected_cores < 1L) {
    detected_cores <- 1L
  }
  run_cores <- min(opts$chains, detected_cores)

  fit_path <- file.path(output_dir, paste0(slug, "_fit.rds"))
  draws_path <- file.path(output_dir, paste0(slug, "_ba_draws.csv"))
  summary_path <- file.path(output_dir, paste0(slug, "_ba_summary.csv"))
  fig_path <- file.path(output_dir, paste0(slug, "_ba_posterior.png"))

  fit <- brm(
    formula = build_formula(),
    data = df,
    family = gaussian(),
    prior = build_priors(),
    chains = opts$chains,
    iter = opts$iter,
    warmup = opts$warmup,
    seed = opts$seed,
    backend = opts$backend,
    cores = run_cores,
    control = list(adapt_delta = opts$adapt_delta, max_treedepth = opts$max_treedepth)
  )
  saveRDS(fit, fit_path)
  diag_df <- write_diagnostics(fit, slug, opts)

  template <- df
  new_e <- template
  new_i <- template
  new_e$vowel <- factor("E", levels = c("E", "I"))
  new_i$vowel <- factor("I", levels = c("E", "I"))

  re_formula <- if (opts$population_only) NA else NULL
  pp_e <- posterior_array(posterior_predict(
    fit,
    newdata = new_e,
    re_formula = re_formula,
    ndraws = opts$posterior_draws,
    seed = opts$seed + 1L
  ))
  pp_i <- posterior_array(posterior_predict(
    fit,
    newdata = new_i,
    re_formula = re_formula,
    ndraws = opts$posterior_draws,
    seed = opts$seed + 2L
  ))

  n_draws <- min(dim(pp_e)[1], dim(pp_i)[1])
  ba <- numeric(n_draws)
  for (draw in seq_len(n_draws)) {
    x <- cbind(pp_e[draw, , 1], pp_e[draw, , 2])
    y <- cbind(pp_i[draw, , 1], pp_i[draw, , 2])
    ba[draw] <- bhatt_affinity(x, y)
  }

  ba_df <- data.frame(draw = seq_along(ba), bhattacharyya_affinity = ba)
  write.csv(ba_df, draws_path, row.names = FALSE)

  ba_summary <- data.frame(
    contrast = "/ɛ/ ~ /ɪ/",
    data = "PB52",
    response_space = "speaker-normalized F1/F2",
    conditional_on_speakers = !opts$population_only,
    n_model_rows = nrow(df),
    n_posterior_draws = length(ba),
    ba_mean = mean(ba),
    ba_median = median(ba),
    ba_q025 = unname(quantile(ba, 0.025)),
    ba_q975 = unname(quantile(ba, 0.975)),
    max_rhat = diag_df$max_rhat,
    min_bulk_ess = diag_df$min_bulk_ess,
    max_treedepth_observed = diag_df$max_treedepth_observed,
    treedepth_limit_hit_count = diag_df$treedepth_limit_hit_count,
    divergences = diag_df$divergences
  )
  write.csv(ba_summary, summary_path, row.names = FALSE)

  png(fig_path, width = 1800, height = 1200, res = 220)
  print(
    ggplot(ba_df, aes(x = bhattacharyya_affinity)) +
      geom_histogram(bins = 40, fill = "#2D6A8E", color = "white", linewidth = 0.2) +
      geom_vline(xintercept = ba_summary$ba_median, color = "#D7263D", linewidth = 1) +
      coord_cartesian(xlim = c(0, 1)) +
      labs(
        x = "Bhattacharyya affinity",
        y = "Posterior draws"
      ) +
      theme_minimal(base_size = 16)
  )
  dev.off()

  cat("Wrote", fit_path, "\n")
  cat("Wrote", draws_path, "\n")
  cat("Wrote", summary_path, "\n")
  cat("Wrote", fig_path, "\n")
  print(ba_summary)
}

opts <- parse_args()
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
slug <- analysis_slug(opts$run_label)
df <- prepare_data()
write_preflight(df, opts, slug)

if (!opts$fit) {
  cat("Preflight only. Re-run with --fit to launch the PB52 brms model.\n")
  quit(save = "no", status = 0)
}

run_fit(df, opts, slug)
