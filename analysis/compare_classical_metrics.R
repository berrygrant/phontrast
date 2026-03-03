## ============================================================
## Compare JSD to classical metrics on consonant contrasts
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

parse_int_env <- function(name, default, min_value = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  val <- suppressWarnings(as.integer(raw))
  if (is.na(val) || val < min_value) return(as.integer(default))
  val
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

jsd_1d <- function(x1, x2, n_grid = 512L) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 2L || length(x2) < 2L) {
    return(list(jsd = NA_real_, overlap = NA_real_))
  }
  rng <- range(c(x1, x2))
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(list(jsd = 0, overlap = 1))
  }
  d1 <- stats::density(x1, n = n_grid, from = rng[1], to = rng[2])
  d2 <- stats::density(x2, n = n_grid, from = rng[1], to = rng[2])
  p <- pmax(d1$y, 0)
  q <- pmax(d2$y, 0)
  p <- p / sum(p)
  q <- q / sum(q)
  list(
    jsd = jsd(p, q),
    overlap = sum(pmin(p, q))
  )
}

density_overlap_nd <- function(X1, X2) {
  X_all <- rbind(X1, X2)
  H1 <- ks::Hpi.diag(X1)
  H2 <- ks::Hpi.diag(X2)
  kde1 <- ks::kde(x = X1, H = H1, eval.points = X_all)
  kde2 <- ks::kde(x = X2, H = H2, eval.points = X_all)
  p <- as.numeric(kde1$estimate)
  q <- as.numeric(kde2$estimate)
  p <- p / sum(p)
  q <- q / sum(q)
  sum(pmin(p, q))
}

jsd_nd <- function(df, features, category_col = "category", bw = "Hpi.diag", eval_on = "pooled") {
  args <- list(
    data = df,
    features = features,
    group = category_col
  )
  if (has_arg(jsd_kde_nd, "bw")) args$bw <- bw
  if (has_arg(jsd_kde_nd, "eval_on")) args$eval_on <- eval_on
  do.call(jsd_kde_nd, args)
}

pillai_value <- function(df, features, category_col = "category") {
  Y <- as.matrix(df[, features, drop = FALSE])
  cat <- factor(df[[category_col]])
  if (nlevels(cat) != 2L) return(NA_real_)
  m <- stats::manova(Y ~ cat)
  s <- summary(m, test = "Pillai")
  as.numeric(s$stats[1, "Pillai"])
}

bhatt_mvnorm <- function(X1, X2, eps = 1e-6) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if (is.null(dim(X1))) X1 <- matrix(X1, ncol = 1)
  if (is.null(dim(X2))) X2 <- matrix(X2, ncol = 1)

  mu1 <- colMeans(X1)
  mu2 <- colMeans(X2)
  S1 <- stats::cov(X1)
  S2 <- stats::cov(X2)
  if (!is.matrix(S1)) S1 <- matrix(S1, nrow = ncol(X1), ncol = ncol(X1))
  if (!is.matrix(S2)) S2 <- matrix(S2, nrow = ncol(X2), ncol = ncol(X2))

  d <- ncol(X1)
  S1 <- S1 + diag(eps, d)
  S2 <- S2 + diag(eps, d)
  S <- (S1 + S2) / 2

  invS <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(invS)) return(NA_real_)

  diff <- matrix(mu2 - mu1, ncol = 1)
  term1 <- as.numeric(0.125 * t(diff) %*% invS %*% diff)

  detS <- tryCatch(determinant(S, logarithm = TRUE), error = function(e) NULL)
  detS1 <- tryCatch(determinant(S1, logarithm = TRUE), error = function(e) NULL)
  detS2 <- tryCatch(determinant(S2, logarithm = TRUE), error = function(e) NULL)
  if (is.null(detS) || is.null(detS1) || is.null(detS2)) return(NA_real_)
  if (detS$sign <= 0 || detS1$sign <= 0 || detS2$sign <= 0) return(NA_real_)

  term2 <- as.numeric(0.5 * (detS$modulus - 0.5 * (detS1$modulus + detS2$modulus)))
  as.numeric(term1 + term2)
}

safe_mahal <- function(mu_diff, S, eps = 1e-6) {
  if (!is.matrix(S)) S <- matrix(S, nrow = length(mu_diff), ncol = length(mu_diff))
  invS <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(invS)) {
    invS <- tryCatch(solve(S + diag(eps, nrow(S))), error = function(e) NULL)
  }
  if (is.null(invS)) return(NA_real_)
  as.numeric(sqrt(t(mu_diff) %*% invS %*% mu_diff))
}

lda_accuracy <- function(df, features, category_col = "category", test_prop = 0.3) {
  df <- df %>% mutate(.id = row_number())
  n <- nrow(df)
  n_test <- max(4L, floor(test_prop * n))
  if (n_test >= n) return(NA_real_)
  test_ids <- sample(df$.id, n_test)
  train <- df %>% filter(!.id %in% test_ids)
  test <- df %>% filter(.id %in% test_ids)
  form <- as.formula(paste(category_col, "~", paste(features, collapse = " + ")))
  fit <- tryCatch(MASS::lda(form, data = train), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  pred <- tryCatch(predict(fit, newdata = test)$class, error = function(e) NULL)
  if (is.null(pred)) return(NA_real_)
  mean(pred == test[[category_col]])
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
min_per_category <- parse_int_env("CLASSICAL_MIN_PER_CATEGORY", 40L, 2L)
max_per_category <- parse_int_env("CLASSICAL_MAX_PER_CATEGORY", 700L, 10L)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- safe_read(data_path)
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

rows <- list()
skipped <- list()

for (spec in contrast_specs) {
  feat <- spec$features
  missing_feat <- setdiff(feat, names(df))
  if (length(missing_feat)) {
    skipped[[length(skipped) + 1L]] <- tibble(
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
    skipped[[length(skipped) + 1L]] <- tibble(
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      cat1 = spec$cat1,
      cat2 = spec$cat2,
      status = "skipped",
      note = paste0("Too few tokens after filtering: ", spec$cat1, "=", n1, ", ", spec$cat2, "=", n2)
    )
    next
  }

  X1 <- as_numeric_matrix(d %>% filter(.data$category == spec$cat1), feat)
  X2 <- as_numeric_matrix(d %>% filter(.data$category == spec$cat2), feat)
  X <- as_numeric_matrix(d, feat)

  if (length(feat) == 1L) {
    jsd_ov <- jsd_1d(X1[, 1], X2[, 1])
    jsd_val <- jsd_ov$jsd
    overlap <- jsd_ov$overlap
  } else {
    jsd_val <- tryCatch(jsd_nd(d, feat, "category", bw, eval_on), error = function(e) NA_real_)
    overlap <- tryCatch(density_overlap_nd(X1, X2), error = function(e) NA_real_)
  }

  pillai <- tryCatch(pillai_value(d, feat, "category"), error = function(e) NA_real_)
  bhatt <- tryCatch(bhatt_mvnorm(X1, X2), error = function(e) NA_real_)
  euclid <- tryCatch(sqrt(sum((colMeans(X2) - colMeans(X1))^2)), error = function(e) NA_real_)
  mahal <- tryCatch({
    S <- stats::cov(X)
    safe_mahal(colMeans(X2) - colMeans(X1), S)
  }, error = function(e) NA_real_)
  lda_acc <- tryCatch(lda_accuracy(d, feat, "category", 0.3), error = function(e) NA_real_)

  rows[[length(rows) + 1L]] <- tibble(
    contrast_id = spec$contrast_id,
    domain = spec$domain,
    cat1 = spec$cat1,
    cat2 = spec$cat2,
    n_cat1 = n1,
    n_cat2 = n2,
    feature_set = paste(feat, collapse = ","),
    jsd = jsd_val,
    pillai = pillai,
    bhatt_dist = bhatt,
    overlap = overlap,
    sep_from_overlap = ifelse(is.finite(overlap), 1 - overlap, NA_real_),
    euclid = euclid,
    mahal = mahal,
    lda_acc = lda_acc
  )
}

metrics_df <- if (length(rows)) bind_rows(rows) else tibble()
skipped_df <- if (length(skipped)) bind_rows(skipped) else tibble()

comparison_path <- file.path(out_dir, "classical_metric_comparison.csv")
skipped_path <- file.path(out_dir, "classical_metric_comparison_skipped.csv")
write.csv(metrics_df, comparison_path, row.names = FALSE)
write.csv(skipped_df, skipped_path, row.names = FALSE)

if (nrow(metrics_df)) {
  scoring_df <- metrics_df %>%
    mutate(expected_sep = ifelse(.data$domain == "stop_voicing", 0, 1))

  metric_cols <- c("jsd", "pillai", "bhatt_dist", "sep_from_overlap", "euclid", "mahal", "lda_acc")
  rank_rows <- lapply(metric_cols, function(m) {
    v <- scoring_df[[m]]
    ok <- is.finite(v) & is.finite(scoring_df$expected_sep)
    rho <- if (sum(ok) >= 3L) suppressWarnings(cor(v[ok], scoring_df$expected_sep[ok], method = "spearman")) else NA_real_
    tibble(metric = m, spearman_with_expected_sep = rho)
  })
  rank_df <- bind_rows(rank_rows) %>%
    arrange(desc(.data$spearman_with_expected_sep))
  ranking_path <- file.path(out_dir, "classical_metric_rankings.csv")
  write.csv(rank_df, ranking_path, row.names = FALSE)
  message("Metric ranking: ", ranking_path)
}

message("Finished metric comparison.")
message("Using local source: ", using_local_source)
message("Comparison table: ", comparison_path)
message("Skipped: ", skipped_path)
