#!/usr/bin/env Rscript

# Compare empirical JSD with PB52 MMO-framework posterior overlap.

suppressPackageStartupMessages({
  library(phonJSD)
  library(phonTools)
  library(dplyr)
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
out_dir <- file.path(analysis_dir, "data")
mmo_summary_path <- file.path(out_dir, "pb52_E_I_mmo_ba_summary.csv")
comparison_out <- file.path(out_dir, "pb52_E_I_jsd_mmo_comparison.csv")
metrics_out <- file.path(out_dir, "pb52_E_I_jsd_mmo_supporting_metrics.csv")

lobanov <- function(x) {
  as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

set.seed(20260625)
data(pb52, package = "phonTools")

df <- pb52 |>
  mutate(speaker = factor(speaker)) |>
  group_by(speaker) |>
  mutate(
    f1_z = lobanov(f1),
    f2_z = lobanov(f2)
  ) |>
  ungroup() |>
  filter(vowel %in% c("E", "I")) |>
  mutate(vowel = droplevels(vowel))

jsd <- estimate_jsd(
  data = df,
  features = c("f1_z", "f2_z"),
  category_col = "vowel",
  do_boot = TRUE,
  n_boot = 1000,
  min_tokens = 20
)

kde_overlap <- estimate_overlap(
  data = df,
  features = c("f1_z", "f2_z"),
  category_col = "vowel",
  min_tokens = 20
)

bhatt <- estimate_bhatt(
  data = df,
  features = c("f1_z", "f2_z"),
  category_col = "vowel",
  min_tokens = 20
)

mmo <- read.csv(mmo_summary_path, stringsAsFactors = FALSE)

comparison <- bind_rows(
  tibble(
    measure = "JSD",
    model = "Empirical KDE",
    estimand = "separation",
    point = jsd$jsd_point,
    low = jsd$jsd_low,
    high = jsd$jsd_high,
    interval = "95% bootstrap CI",
    n_tokens = jsd$n_tokens,
    note = "JSD is an information-theoretic separation index: 0 = identical, 1 = maximally distinct."
  ),
  tibble(
    measure = "1 - MMO BA",
    model = "PB52 MMO-framework",
    estimand = "separation",
    point = 1 - mmo$ba_median,
    low = 1 - mmo$ba_q975,
    high = 1 - mmo$ba_q025,
    interval = "95% posterior CrI",
    n_tokens = mmo$n_model_rows,
    note = "MMO returns Bhattacharyya affinity overlap, so 1 - BA is shown to put it on a separation scale."
  ),
  tibble(
    measure = "KDE overlap",
    model = "Empirical KDE",
    estimand = "overlap",
    point = kde_overlap$overlap,
    low = NA_real_,
    high = NA_real_,
    interval = NA_character_,
    n_tokens = kde_overlap$n_tokens,
    note = "Reference empirical overlap used for the main PB52 metric comparison."
  ),
  tibble(
    measure = "MMO BA",
    model = "PB52 MMO-framework",
    estimand = "overlap",
    point = mmo$ba_median,
    low = mmo$ba_q025,
    high = mmo$ba_q975,
    interval = "95% posterior CrI",
    n_tokens = mmo$n_model_rows,
    note = "Posterior predictive Bhattacharyya affinity from the brms mixed-effects model."
  )
)

supporting <- tibble(
  n_tokens = jsd$n_tokens,
  jsd_point = jsd$jsd_point,
  jsd_low = jsd$jsd_low,
  jsd_high = jsd$jsd_high,
  kde_overlap = kde_overlap$overlap,
  bhatt_affinity_empirical = bhatt$bhatt_affinity,
  mmo_ba_median = mmo$ba_median,
  mmo_ba_low = mmo$ba_q025,
  mmo_ba_high = mmo$ba_q975,
  mmo_separation = 1 - mmo$ba_median,
  mmo_separation_low = 1 - mmo$ba_q975,
  mmo_separation_high = 1 - mmo$ba_q025
)

write.csv(comparison, comparison_out, row.names = FALSE)
write.csv(supporting, metrics_out, row.names = FALSE)

cat("Wrote", comparison_out, "\n")
cat("Wrote", metrics_out, "\n")
print(comparison)
print(supporting)
