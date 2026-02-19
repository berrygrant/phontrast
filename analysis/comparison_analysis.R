## ============================================================
## LabPhon 2026: Pairwise metric comparison in 13D MFCC space
## ============================================================

set.seed(2026)

# Core packages
library(phonJSD)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(tidyquant)

## ============================================================
## 0. Input configuration
## ============================================================

# Provide a path to a CSV or RDS that contains vowel labels and either:
# - MFCC columns named mfcc1..mfcc13 (mode = "mfcc_ready"), or
# - audio paths + segment boundaries so we can extract MFCCs (mode = "extract_mfcc")
#
# Example:
# data_path <- "/path/to/your/mfcc_data.rds"

data_path <- Sys.getenv("MFCC_DATA_PATH", unset = "")
category_col <- "vowel"
mode <- "mfcc_ready"  # "mfcc_ready" | "extract_mfcc"

# If mode == "extract_mfcc", set these columns
file_col  <- "file"
start_col <- "t_start"   # seconds
end_col   <- "t_end"     # seconds

# MFCC configuration
numcep <- 13
mfcc_cols <- paste0("mfcc", 1:numcep)

# Sampling controls for tractable pairwise KDE on large corpora
min_tokens_per_vowel <- 200
max_per_vowel <- 500

# KDE settings for speed (diagonal Scott's rule + subsampled eval points)
bw_method <- "scott"      # "scott" | "Hpi" | "Hscv" | "Hpi.diag"
eval_n <- 200             # number of pooled points to evaluate KDE on

# Bhattacharyya ridge constant for numerical stability
bhatt_eps <- 1e-2

if (data_path == "") {
  stop("Set `data_path` to a CSV or RDS with your vowel data.")
}

# Load data
if (grepl("\\.rds$", data_path, ignore.case = TRUE)) {
  df <- readRDS(data_path)
} else {
  df <- read.csv(data_path, stringsAsFactors = FALSE)
}

# Extract MFCCs if needed
if (mode == "extract_mfcc") {
  df <- extract_mfcc(
    data      = df,
    file_col  = file_col,
    start_col = start_col,
    end_col   = end_col,
    numcep    = numcep,
    prefix    = "mfcc"
  )
}

# Validate required columns
missing_cols <- setdiff(c(category_col, mfcc_cols), names(df))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Drop incomplete rows
keep_cols <- c(category_col, mfcc_cols)
df <- df[stats::complete.cases(df[, keep_cols, drop = FALSE]), , drop = FALSE]

# Filter to sufficiently frequent vowels and sample per vowel
df <- df %>%
  group_by(.data[[category_col]]) %>%
  filter(dplyr::n() >= min_tokens_per_vowel) %>%
  group_modify(~ {
    if (nrow(.x) > max_per_vowel) {
      dplyr::slice_sample(.x, n = max_per_vowel)
    } else {
      .x
    }
  }) %>%
  ungroup()

## ============================================================
## 1. Pairwise metrics (all vowel contrasts)
## ============================================================

scott_H <- function(X) {
  d <- ncol(X)
  n <- nrow(X)
  sds <- apply(X, 2, stats::sd)
  h <- n ^ (-1 / (d + 4))
  diag((h * sds) ^ 2)
}

jsd_kde_fast <- function(df_pair,
                         features,
                         category_col,
                         bw_method = "scott",
                         eval_n = 200) {
  levs <- unique(df_pair[[category_col]])
  if (length(levs) != 2) return(NA_real_)
  d1 <- df_pair[df_pair[[category_col]] == levs[1], , drop = FALSE]
  d2 <- df_pair[df_pair[[category_col]] == levs[2], , drop = FALSE]
  X1 <- as.matrix(d1[, features, drop = FALSE])
  X2 <- as.matrix(d2[, features, drop = FALSE])
  X_all <- as.matrix(df_pair[, features, drop = FALSE])
  if (nrow(X_all) > eval_n) {
    X_all <- X_all[sample.int(nrow(X_all), eval_n), , drop = FALSE]
  }
  H1 <- switch(bw_method,
               scott = scott_H(X1),
               Hpi = ks::Hpi(X1),
               Hscv = ks::Hscv(X1),
               Hpi.diag = ks::Hpi.diag(X1))
  H2 <- switch(bw_method,
               scott = scott_H(X2),
               Hpi = ks::Hpi(X2),
               Hscv = ks::Hscv(X2),
               Hpi.diag = ks::Hpi.diag(X2))
  kde1 <- ks::kde(x = X1, H = H1, eval.points = X_all)
  kde2 <- ks::kde(x = X2, H = H2, eval.points = X_all)
  p <- as.numeric(kde1$estimate); q <- as.numeric(kde2$estimate)
  p <- p / sum(p); q <- q / sum(q)
  jsd(p, q)
}

overlap_kde_fast <- function(df_pair,
                             features,
                             category_col,
                             bw_method = "scott",
                             eval_n = 200) {
  levs <- unique(df_pair[[category_col]])
  if (length(levs) != 2) return(NA_real_)
  d1 <- df_pair[df_pair[[category_col]] == levs[1], , drop = FALSE]
  d2 <- df_pair[df_pair[[category_col]] == levs[2], , drop = FALSE]
  X1 <- as.matrix(d1[, features, drop = FALSE])
  X2 <- as.matrix(d2[, features, drop = FALSE])
  X_all <- as.matrix(df_pair[, features, drop = FALSE])
  if (nrow(X_all) > eval_n) {
    X_all <- X_all[sample.int(nrow(X_all), eval_n), , drop = FALSE]
  }
  H1 <- switch(bw_method,
               scott = scott_H(X1),
               Hpi = ks::Hpi(X1),
               Hscv = ks::Hscv(X1),
               Hpi.diag = ks::Hpi.diag(X1))
  H2 <- switch(bw_method,
               scott = scott_H(X2),
               Hpi = ks::Hpi(X2),
               Hscv = ks::Hscv(X2),
               Hpi.diag = ks::Hpi.diag(X2))
  kde1 <- ks::kde(x = X1, H = H1, eval.points = X_all)
  kde2 <- ks::kde(x = X2, H = H2, eval.points = X_all)
  p <- as.numeric(kde1$estimate); q <- as.numeric(kde2$estimate)
  p <- p / sum(p); q <- q / sum(q)
  sum(pmin(p, q))
}

bhatt_mvnorm_fast <- function(df_pair,
                              features,
                              category_col,
                              eps = 1e-2) {
  levs <- unique(df_pair[[category_col]])
  if (length(levs) != 2) return(NA_real_)
  X1 <- as.matrix(df_pair[df_pair[[category_col]] == levs[1], features, drop = FALSE])
  X2 <- as.matrix(df_pair[df_pair[[category_col]] == levs[2], features, drop = FALSE])
  mu1 <- colMeans(X1); mu2 <- colMeans(X2)
  S1 <- stats::cov(X1); S2 <- stats::cov(X2)
  if (any(!is.finite(S1)) || any(!is.finite(S2))) return(NA_real_)
  d <- ncol(X1)
  S1 <- S1 + diag(eps, d)
  S2 <- S2 + diag(eps, d)
  S <- (S1 + S2) / 2

  invS <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(invS)) return(NA_real_)

  diff <- matrix(mu2 - mu1, ncol = 1)
  term1 <- 0.125 * t(diff) %*% invS %*% diff

  detS  <- determinant(S, logarithm = TRUE)
  detS1 <- determinant(S1, logarithm = TRUE)
  detS2 <- determinant(S2, logarithm = TRUE)
  if (detS$sign <= 0 || detS1$sign <= 0 || detS2$sign <= 0) return(NA_real_)

  term2 <- 0.5 * (detS$modulus - 0.5 * (detS1$modulus + detS2$modulus))
  as.numeric(term1 + term2)
}

compute_metrics_pair <- function(df_pair,
                                 features,
                                 category_col = "vowel",
                                 bhatt_eps = 1e-6,
                                 bw_method = "scott",
                                 eval_n = 200) {
  tibble(
    jsd = tryCatch(
      jsd_kde_fast(df_pair, features, category_col, bw_method, eval_n),
      error = function(e) NA_real_
    ),
    pillai = tryCatch(
      estimate_pillai(df_pair, features, category_col = category_col)$pillai,
      error = function(e) NA_real_
    ),
    bhatt = tryCatch(
      bhatt_mvnorm_fast(df_pair, features, category_col, bhatt_eps),
      error = function(e) NA_real_
    ),
    overlap = tryCatch(
      overlap_kde_fast(df_pair, features, category_col, bw_method, eval_n),
      error = function(e) NA_real_
    )
  )
}

vowels <- sort(unique(df[[category_col]]))
pairs <- combn(vowels, 2, simplify = FALSE)

## ============================================================
## 2. Outputs (set up early for incremental writes)
## ============================================================

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

out_dir <- file.path(script_dir, "comparison_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pair_results_path <- file.path(out_dir, "pairwise_metrics_mfcc13.csv")

# Load any existing results to allow resume
done_pairs <- NULL
if (file.exists(pair_results_path)) {
  existing <- read.csv(pair_results_path, stringsAsFactors = FALSE)
  done_pairs <- paste(existing$v1, existing$v2, sep = "_")
} else {
  # write header
  write.csv(
    data.frame(
      v1 = character(0),
      v2 = character(0),
      n_tokens = integer(0),
      jsd = numeric(0),
      pillai = numeric(0),
      bhatt = numeric(0),
      overlap = numeric(0)
    ),
    pair_results_path,
    row.names = FALSE
  )
}

pair_results <- list()
counter <- 0

for (p in pairs) {
  v1 <- p[1]
  v2 <- p[2]
  key <- paste(v1, v2, sep = "_")
  if (!is.null(done_pairs) && key %in% done_pairs) {
    next
  }

  df_pair <- df %>%
    filter(.data[[category_col]] %in% c(v1, v2)) %>%
    droplevels()

  metrics <- compute_metrics_pair(
    df_pair,
    features = mfcc_cols,
    category_col = category_col,
    bhatt_eps = bhatt_eps,
    bw_method = bw_method,
    eval_n = eval_n
  )

  row <- tibble(v1 = v1, v2 = v2, n_tokens = nrow(df_pair)) %>%
    bind_cols(metrics)

  pair_results[[length(pair_results) + 1]] <- row
  counter <- counter + 1

  # flush every 25 pairs
  if (counter %% 25 == 0) {
    out <- dplyr::bind_rows(pair_results)
    suppressWarnings(
      write.table(out, pair_results_path, sep = ",", row.names = FALSE,
                  col.names = FALSE, append = TRUE)
    )
    pair_results <- list()
  }
}

# flush remainder
if (length(pair_results)) {
  out <- dplyr::bind_rows(pair_results)
  suppressWarnings(
    write.table(out, pair_results_path, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE)
  )
}

# Read full results
pair_results <- read.csv(pair_results_path, stringsAsFactors = FALSE)

## ============================================================
## 3. Correlations with overlap
## ============================================================

cor_tbl <- tibble(
  metric = c("JSD", "Pillai", "Bhatt"),
  pearson = c(
    {
      ok <- is.finite(pair_results$jsd) & is.finite(pair_results$overlap)
      if (sum(ok) < 3) NA_real_ else cor(pair_results$jsd[ok], pair_results$overlap[ok])
    },
    {
      ok <- is.finite(pair_results$pillai) & is.finite(pair_results$overlap)
      if (sum(ok) < 3) NA_real_ else cor(pair_results$pillai[ok], pair_results$overlap[ok])
    },
    {
      ok <- is.finite(pair_results$bhatt) & is.finite(pair_results$overlap)
      if (sum(ok) < 3) NA_real_ else cor(pair_results$bhatt[ok], pair_results$overlap[ok])
    }
  ),
  spearman = c(
    {
      ok <- is.finite(pair_results$jsd) & is.finite(pair_results$overlap)
      if (sum(ok) < 3) NA_real_ else cor(pair_results$jsd[ok], pair_results$overlap[ok], method = "spearman")
    },
    {
      ok <- is.finite(pair_results$pillai) & is.finite(pair_results$overlap)
      if (sum(ok) < 3) NA_real_ else cor(pair_results$pillai[ok], pair_results$overlap[ok], method = "spearman")
    },
    {
      ok <- is.finite(pair_results$bhatt) & is.finite(pair_results$overlap)
      if (sum(ok) < 3) NA_real_ else cor(pair_results$bhatt[ok], pair_results$overlap[ok], method = "spearman")
    }
  )
)

## ============================================================
## 4. Outputs
## ============================================================

write.csv(cor_tbl, file.path(out_dir, "metric_overlap_correlations.csv"), row.names = FALSE)

## ============================================================
## 5. Visualizations (tidyquant palettes)
## ============================================================

# 4a. Metric vs overlap (scatter + linear fit)
plot_df <- pair_results %>%
  pivot_longer(cols = c(jsd, pillai, bhatt), names_to = "metric", values_to = "value")

p_scatter <- ggplot(plot_df, aes(x = value, y = overlap, color = metric)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(
    x = "Metric value",
    y = "KDE overlap",
    title = "Metric vs overlap across all vowel pairs (13D MFCC)"
  ) +
  theme_tq() +
  scale_color_tq()

# 4b. Correlation summary (Pearson + Spearman)
cor_long <- cor_tbl %>%
  pivot_longer(cols = c(pearson, spearman), names_to = "type", values_to = "r")

p_corr <- ggplot(cor_long, aes(x = metric, y = r, fill = type)) +
  geom_col(position = "dodge") +
  labs(
    x = "Metric",
    y = "Correlation with overlap",
    title = "Correlation of metrics with overlap (13D MFCC)"
  ) +
  theme_tq() +
  scale_fill_tq()

# 4c. Overlap heatmap (pairwise)
pair_sym <- pair_results %>%
  select(v1, v2, overlap) %>%
  bind_rows(pair_results %>% transmute(v1 = v2, v2 = v1, overlap = overlap)) %>%
  bind_rows(tibble(v1 = vowels, v2 = vowels, overlap = 1))

p_heat <- ggplot(pair_sym, aes(x = v1, y = v2, fill = overlap)) +
  geom_tile(color = "white", linewidth = 0.2) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL,
    title = "Pairwise overlap (13D MFCC)",
    fill = "Overlap"
  ) +
  theme_tq() +
  scale_fill_gradientn(colors = tidyquant::palette_light())

# Save plots
ggsave(file.path(out_dir, "scatter_metric_vs_overlap.png"), p_scatter, width = 10, height = 4, dpi = 300)
ggsave(file.path(out_dir, "correlation_summary.png"), p_corr, width = 7, height = 4, dpi = 300)
ggsave(file.path(out_dir, "overlap_heatmap.png"), p_heat, width = 6, height = 6, dpi = 300)

## ============================================================
## 6. Quick diagnostics (top / bottom overlap pairs)
## ============================================================

pair_results %>%
  mutate(pair = paste(v1, v2, sep = "-")) %>%
  arrange(overlap) %>%
  slice_head(n = 5)

pair_results %>%
  mutate(pair = paste(v1, v2, sep = "-")) %>%
  arrange(desc(overlap)) %>%
  slice_head(n = 5)
