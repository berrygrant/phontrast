#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || identical(x, "")) y else x
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript analysis/validation/run_validation_metrics.R \\",
      "    --input token_features.csv --out-dir outputs/validation \\",
      "    [--feature-sets feature_sets.csv | --features f1,f2,f3] \\",
      "    [--domain tone] [--domain-col domain] [--category-col category] \\",
      "    [--control-col control_group] [--min-per-category 20] \\",
      "    [--bw scott.diag] [--eval-on pooled_sample] [--eval-n 300] \\",
      "    [--engine fast_diag] [--cv-folds 5] [--seed 2026]",
      "",
      "Writes validation_metrics.csv, validation_skipped_contrasts.csv, and",
      "validation_run_summary.csv to --out-dir.",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    key <- gsub("-", "_", key)
    if (identical(key, "help")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

split_feature_string <- function(x) {
  x <- paste(x, collapse = " ")
  out <- unlist(strsplit(x, "[,;[:space:]]+"), use.names = FALSE)
  out[nzchar(out)]
}

read_input_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "rds")) {
    obj <- readRDS(path)
    if (!is.data.frame(obj)) {
      stop("RDS input must contain a data frame.", call. = FALSE)
    }
    return(as.data.frame(obj, stringsAsFactors = FALSE))
  }
  if (identical(ext, "tsv")) {
    return(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (identical(ext, "csv")) {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop("Unsupported input extension: ", ext, ". Use CSV, TSV, or RDS.", call. = FALSE)
}

read_feature_sets <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "tsv")) {
    fs <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    fs <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  required <- c("feature_set", "features")
  missing <- setdiff(required, names(fs))
  if (length(missing)) {
    stop("Feature-set file is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  fs <- fs[nzchar(fs$feature_set) & nzchar(fs$features), , drop = FALSE]
  if (!nrow(fs)) {
    stop("Feature-set file contains no usable rows.", call. = FALSE)
  }
  lapply(seq_len(nrow(fs)), function(i) {
    list(
      name = as.character(fs$feature_set[[i]]),
      features = split_feature_string(fs$features[[i]])
    )
  })
}

infer_numeric_features <- function(df, args) {
  metadata_cols <- unique(c(
    "token_id", "domain", "language", "source_corpus", "speaker",
    "category", "control_group", "file", "start", "end", "word", "phone",
    "syllable", "rime", "preceding_phone", "following_phone",
    "measurement_method", "quality_flag",
    args$domain_col, args$category_col, args$control_col
  ))
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  setdiff(numeric_cols, metadata_cols)
}

load_feature_sets <- function(df, args) {
  if (!is.null(args$feature_sets)) {
    feature_sets <- read_feature_sets(args$feature_sets)
  } else if (!is.null(args$features)) {
    feature_sets <- list(list(name = "all_features", features = split_feature_string(args$features)))
  } else {
    inferred <- infer_numeric_features(df, args)
    if (!length(inferred)) {
      stop(
        "No feature sets supplied and no numeric acoustic feature columns could be inferred.",
        call. = FALSE
      )
    }
    feature_sets <- list(list(name = "all_numeric_features", features = inferred))
  }

  for (fs in feature_sets) {
    if (!length(fs$features)) {
      stop("Feature set `", fs$name, "` has no features.", call. = FALSE)
    }
    missing <- setdiff(fs$features, names(df))
    if (length(missing)) {
      stop(
        "Feature set `", fs$name, "` references missing columns: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    non_numeric <- fs$features[!vapply(df[fs$features], is.numeric, logical(1))]
    if (length(non_numeric)) {
      stop(
        "Feature set `", fs$name, "` contains non-numeric columns: ",
        paste(non_numeric, collapse = ", "),
        call. = FALSE
      )
    }
  }
  feature_sets
}

complete_finite_cases <- function(df, cols) {
  keep <- stats::complete.cases(df[, cols, drop = FALSE])
  numeric_cols <- cols[vapply(df[, cols, drop = FALSE], is.numeric, logical(1))]
  if (length(numeric_cols)) {
    finite <- Reduce(
      `&`,
      lapply(numeric_cols, function(col) is.finite(df[[col]])),
      init = rep(TRUE, nrow(df))
    )
    keep <- keep & finite
  }
  df[keep, , drop = FALSE]
}

rank_auc <- function(labels, scores) {
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]
  scores <- scores[ok]
  if (length(unique(labels)) != 2L) {
    return(NA_real_)
  }
  positive <- sort(unique(labels))[[2L]]
  y <- labels == positive
  n_pos <- sum(y)
  n_neg <- sum(!y)
  if (!n_pos || !n_neg) {
    return(NA_real_)
  }
  ranks <- rank(scores, ties.method = "average")
  as.numeric((sum(ranks[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
}

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  folds <- integer(length(y))
  for (level in unique(y)) {
    idx <- which(y == level)
    idx <- sample(idx)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

scale_by_train <- function(train_x, test_x) {
  center <- vapply(train_x, mean, numeric(1), na.rm = TRUE)
  scale <- vapply(train_x, stats::sd, numeric(1), na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = as.data.frame(scale(train_x, center = center, scale = scale)),
    test = as.data.frame(scale(test_x, center = center, scale = scale))
  )
}

classifier_reference <- function(df, features, category_col, cv_folds = 5L, seed = 2026L) {
  df <- complete_finite_cases(df, c(category_col, features))
  cats <- sort(unique(as.character(df[[category_col]])))
  if (length(cats) != 2L) {
    return(list(
      values = data.frame(
        classifier_auc = NA_real_,
        classifier_balanced_accuracy = NA_real_,
        classifier_accuracy = NA_real_,
        classifier_n_folds = 0L,
        classifier_model = "binomial_glm"
      ),
      note = "classifier skipped: contrast does not have exactly two categories"
    ))
  }
  y <- factor(as.character(df[[category_col]]), levels = cats)
  counts <- table(y)
  k <- min(as.integer(cv_folds), as.integer(min(counts)))
  if (k < 2L) {
    return(list(
      values = data.frame(
        classifier_auc = NA_real_,
        classifier_balanced_accuracy = NA_real_,
        classifier_accuracy = NA_real_,
        classifier_n_folds = 0L,
        classifier_model = "binomial_glm"
      ),
      note = "classifier skipped: fewer than two tokens in at least one category"
    ))
  }

  folds <- make_stratified_folds(y, k, seed)
  scores <- rep(NA_real_, length(y))
  predictions <- rep(NA_character_, length(y))
  positive <- cats[[2L]]

  for (fold in seq_len(k)) {
    train_idx <- folds != fold
    test_idx <- folds == fold
    scaled <- scale_by_train(df[train_idx, features, drop = FALSE], df[test_idx, features, drop = FALSE])
    train_df <- scaled$train
    train_df$.y <- as.integer(y[train_idx] == positive)
    test_df <- scaled$test
    fit <- tryCatch(
      suppressWarnings(stats::glm(.y ~ ., data = train_df, family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    fold_scores <- tryCatch(
      suppressWarnings(stats::predict(fit, newdata = test_df, type = "response")),
      error = function(e) rep(NA_real_, sum(test_idx))
    )
    scores[test_idx] <- as.numeric(fold_scores)
    predictions[test_idx] <- ifelse(scores[test_idx] >= 0.5, positive, cats[[1L]])
  }

  ok <- !is.na(predictions) & is.finite(scores)
  if (!any(ok)) {
    return(list(
      values = data.frame(
        classifier_auc = NA_real_,
        classifier_balanced_accuracy = NA_real_,
        classifier_accuracy = NA_real_,
        classifier_n_folds = k,
        classifier_model = "binomial_glm"
      ),
      note = "classifier failed: all folds failed to fit or predict"
    ))
  }

  y_ok <- y[ok]
  pred_ok <- factor(predictions[ok], levels = cats)
  sensitivity <- mean(pred_ok[y_ok == positive] == positive)
  specificity <- mean(pred_ok[y_ok != positive] != positive)
  values <- data.frame(
    classifier_auc = rank_auc(as.character(y_ok), scores[ok]),
    classifier_balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    classifier_accuracy = mean(pred_ok == y_ok),
    classifier_n_folds = k,
    classifier_model = "binomial_glm"
  )
  list(values = values, note = NA_character_)
}

empty_skipped <- function() {
  data.frame(
    domain = character(),
    language = character(),
    source_corpus = character(),
    control_group = character(),
    cat1 = character(),
    cat2 = character(),
    feature_set = character(),
    n_tokens = integer(),
    n_cat1 = integer(),
    n_cat2 = integer(),
    status = character(),
    note = character(),
    stringsAsFactors = FALSE
  )
}

empty_metrics <- function() {
  data.frame(
    domain = character(),
    language = character(),
    source_corpus = character(),
    control_group = character(),
    cat1 = character(),
    cat2 = character(),
    feature_set = character(),
    feature_count = integer(),
    features = character(),
    n_cat1 = integer(),
    n_cat2 = integer(),
    scope = character(),
    n_tokens = integer(),
    pillai = numeric(),
    pillai_p_value = numeric(),
    bhatt_dist = numeric(),
    bhatt_affinity = numeric(),
    jsd = numeric(),
    js_distance = numeric(),
    mahalanobis_dist = numeric(),
    percent_overlap = numeric(),
    metric_mode = character(),
    classifier_auc = numeric(),
    classifier_balanced_accuracy = numeric(),
    classifier_accuracy = numeric(),
    classifier_n_folds = integer(),
    classifier_model = character(),
    warning_note = character(),
    stringsAsFactors = FALSE
  )
}

make_skip <- function(meta, feature_set, n_tokens, n_cat1, n_cat2, note) {
  data.frame(
    domain = meta$domain,
    language = meta$language,
    source_corpus = meta$source_corpus,
    control_group = meta$control_group,
    cat1 = meta$cat1,
    cat2 = meta$cat2,
    feature_set = feature_set,
    n_tokens = as.integer(n_tokens),
    n_cat1 = as.integer(n_cat1),
    n_cat2 = as.integer(n_cat2),
    status = "skipped",
    note = note,
    stringsAsFactors = FALSE
  )
}

safe_value <- function(df, col, default = "unspecified") {
  if (!col %in% names(df)) {
    return(default)
  }
  values <- unique(as.character(df[[col]]))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) {
    return(default)
  }
  if (length(values) == 1L) {
    return(values)
  }
  "multiple"
}

jsd_overlap_fallback <- function(df,
                                 features,
                                 category_col,
                                 min_tokens,
                                 bw,
                                 eval_on,
                                 eval_n,
                                 eval_seed,
                                 engine,
                                 chunk_size) {
  jsd_out <- phonJSD::estimate_jsd(
    data = df,
    features = features,
    category_col = category_col,
    min_tokens = min_tokens,
    do_boot = FALSE,
    bw = bw,
    eval_on = eval_on,
    eval_n = eval_n,
    eval_seed = eval_seed,
    engine = engine,
    chunk_size = chunk_size
  )

  overlap_out <- tryCatch(
    phonJSD::estimate_overlap(
      data = df,
      features = features,
      category_col = category_col,
      min_tokens = min_tokens,
      bw = bw,
      eval_on = eval_on,
      eval_n = eval_n,
      eval_seed = eval_seed,
      engine = engine,
      chunk_size = chunk_size
    ),
    error = function(e) NULL
  )
  overlap <- if (is.null(overlap_out) || !nrow(overlap_out)) NA_real_ else overlap_out$overlap[[1L]]

  data.frame(
    scope = "global",
    n_tokens = jsd_out$n_tokens[[1L]],
    pillai = NA_real_,
    pillai_p_value = NA_real_,
    bhatt_dist = NA_real_,
    bhatt_affinity = NA_real_,
    jsd = jsd_out$jsd_point[[1L]],
    js_distance = sqrt(jsd_out$jsd_point[[1L]]),
    mahalanobis_dist = NA_real_,
    percent_overlap = overlap,
    metric_mode = "jsd_overlap_fallback",
    stringsAsFactors = FALSE
  )
}

run_one_contrast <- function(df,
                             feature_set,
                             category_col,
                             meta,
                             min_per_category,
                             min_tokens,
                             bw,
                             eval_on,
                             eval_n,
                             eval_seed,
                             engine,
                             chunk_size,
                             eps,
                             cv_folds,
                             seed) {
  features <- feature_set$features
  df <- complete_finite_cases(df, c(category_col, features))
  counts <- table(factor(as.character(df[[category_col]]), levels = c(meta$cat1, meta$cat2)))
  n_cat1 <- as.integer(counts[[meta$cat1]])
  n_cat2 <- as.integer(counts[[meta$cat2]])
  n_tokens <- nrow(df)

  if (n_tokens < min_tokens) {
    return(list(
      metrics = NULL,
      skipped = make_skip(meta, feature_set$name, n_tokens, n_cat1, n_cat2, "not enough complete tokens")
    ))
  }
  if (n_cat1 < min_per_category || n_cat2 < min_per_category) {
    return(list(
      metrics = NULL,
      skipped = make_skip(meta, feature_set$name, n_tokens, n_cat1, n_cat2, "not enough complete tokens per category")
    ))
  }

  metric_warnings <- character()
  metrics <- tryCatch(
    withCallingHandlers(
      phonJSD::compare_overlap_metrics(
        data = df,
        features = features,
        category_col = category_col,
        min_tokens = min_tokens,
        bw = bw,
        eval_on = eval_on,
        eval_n = eval_n,
        eval_seed = eval_seed,
        engine = engine,
        chunk_size = chunk_size,
        eps = eps,
        output = "wide",
        do_boot = FALSE,
        progress = FALSE
      ),
      warning = function(w) {
        metric_warnings <<- c(metric_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  metric_mode_note <- NA_character_
  if (inherits(metrics, "error")) {
    compare_error <- conditionMessage(metrics)
    metrics <- tryCatch(
      jsd_overlap_fallback(
        df = df,
        features = features,
        category_col = category_col,
        min_tokens = min_tokens,
        bw = bw,
        eval_on = eval_on,
        eval_n = eval_n,
        eval_seed = eval_seed,
        engine = engine,
        chunk_size = chunk_size
      ),
      error = function(e) e
    )
    if (inherits(metrics, "error")) {
      return(list(
        metrics = NULL,
        skipped = make_skip(
          meta,
          feature_set$name,
          n_tokens,
          n_cat1,
          n_cat2,
          paste(
            "metric estimation failed:",
            compare_error,
            "| fallback failed:",
            conditionMessage(metrics)
          )
        )
      ))
    }
    metric_mode_note <- paste(
      "classical metric comparison failed; reported JSD/overlap fallback:",
      compare_error
    )
  } else {
    metrics$metric_mode <- "all_metrics"
  }
  if (!nrow(metrics)) {
    return(list(
      metrics = NULL,
      skipped = make_skip(
        meta,
        feature_set$name,
        n_tokens,
        n_cat1,
        n_cat2,
        paste(c("metric estimation returned no rows", metric_warnings), collapse = " | ")
      )
    ))
  }

  reference <- classifier_reference(
    df = df,
    features = features,
    category_col = category_col,
    cv_folds = cv_folds,
    seed = seed
  )

  out <- cbind(
    data.frame(
      domain = meta$domain,
      language = meta$language,
      source_corpus = meta$source_corpus,
      control_group = meta$control_group,
      cat1 = meta$cat1,
      cat2 = meta$cat2,
      feature_set = feature_set$name,
      feature_count = length(features),
      features = paste(features, collapse = ";"),
      n_cat1 = n_cat1,
      n_cat2 = n_cat2,
      stringsAsFactors = FALSE
    ),
    metrics,
    reference$values
  )
  notes <- c(metric_warnings, metric_mode_note, reference$note)
  notes <- notes[!is.na(notes) & nzchar(notes)]
  out$warning_note <- if (length(notes)) paste(notes, collapse = " | ") else NA_character_

  list(metrics = out, skipped = NULL)
}

enumerate_contrasts <- function(df, args) {
  category_col <- args$category_col
  control_col <- args$control_col
  if (!category_col %in% names(df)) {
    stop("Missing category column: ", category_col, call. = FALSE)
  }

  if (!is.null(control_col) && nzchar(control_col) && control_col %in% names(df)) {
    control_values <- sort(unique(as.character(df[[control_col]])))
    control_values <- control_values[!is.na(control_values) & nzchar(control_values)]
    groups <- setNames(
      lapply(control_values, function(value) df[as.character(df[[control_col]]) == value, , drop = FALSE]),
      control_values
    )
  } else {
    groups <- list(pooled = df)
  }

  out <- list()
  for (control_value in names(groups)) {
    df_g <- groups[[control_value]]
    cats <- sort(unique(as.character(df_g[[category_col]])))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) < 2L) {
      next
    }
    pairs <- utils::combn(cats, 2, simplify = FALSE)
    for (pair in pairs) {
      df_pair <- df_g[as.character(df_g[[category_col]]) %in% pair, , drop = FALSE]
      meta <- list(
        domain = safe_value(df_pair, args$domain_col, args$domain %||% "all"),
        language = safe_value(df_pair, "language"),
        source_corpus = safe_value(df_pair, "source_corpus"),
        control_group = control_value,
        cat1 = pair[[1L]],
        cat2 = pair[[2L]]
      )
      out[[length(out) + 1L]] <- list(data = df_pair, meta = meta)
    }
  }
  out
}

write_outputs <- function(metrics, skipped, out_dir, summary) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  metrics_path <- file.path(out_dir, "validation_metrics.csv")
  skipped_path <- file.path(out_dir, "validation_skipped_contrasts.csv")
  summary_path <- file.path(out_dir, "validation_run_summary.csv")

  if (is.null(metrics) || !nrow(metrics)) {
    metrics <- empty_metrics()
  }
  if (is.null(skipped) || !nrow(skipped)) {
    skipped <- empty_skipped()
  }

  write.csv(metrics, metrics_path, row.names = FALSE)
  write.csv(skipped, skipped_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE)
  invisible(list(metrics_path = metrics_path, skipped_path = skipped_path, summary_path = summary_path))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(0L))
  }

  if (is.null(args$input)) {
    usage()
    stop("Missing --input.", call. = FALSE)
  }

  args$out_dir <- args$out_dir %||% file.path("analysis", "validation", "outputs")
  args$domain_col <- args$domain_col %||% "domain"
  args$category_col <- args$category_col %||% "category"
  args$control_col <- args$control_col %||% "control_group"
  args$min_per_category <- as.integer(args$min_per_category %||% 20L)
  args$min_tokens <- as.integer(args$min_tokens %||% (2L * args$min_per_category))
  args$bw <- args$bw %||% "scott.diag"
  args$eval_on <- args$eval_on %||% "pooled_sample"
  args$eval_n <- as.integer(args$eval_n %||% 300L)
  args$eval_seed <- as.integer(args$eval_seed %||% 2026L)
  args$engine <- args$engine %||% "fast_diag"
  args$chunk_size <- as.integer(args$chunk_size %||% 1000L)
  args$eps <- as.numeric(args$eps %||% 1e-6)
  args$cv_folds <- as.integer(args$cv_folds %||% 5L)
  args$seed <- as.integer(args$seed %||% 2026L)

  if (!requireNamespace("phonJSD", quietly = TRUE)) {
    stop("The phonJSD package must be installed or available on .libPaths().", call. = FALSE)
  }

  df <- read_input_table(args$input)
  if (!nrow(df)) {
    stop("Input table has no rows.", call. = FALSE)
  }

  if (!is.null(args$domain) && args$domain_col %in% names(df)) {
    df <- df[as.character(df[[args$domain_col]]) == args$domain, , drop = FALSE]
    if (!nrow(df)) {
      stop("No rows remain after filtering ", args$domain_col, " == ", args$domain, ".", call. = FALSE)
    }
  }

  feature_sets <- load_feature_sets(df, args)
  contrasts <- enumerate_contrasts(df, args)
  if (!length(contrasts)) {
    stop("No binary contrasts could be enumerated.", call. = FALSE)
  }

  set.seed(args$seed)
  metric_rows <- list()
  skipped_rows <- list()
  counter <- 0L
  total <- length(contrasts) * length(feature_sets)

  for (contrast in contrasts) {
    for (feature_set in feature_sets) {
      counter <- counter + 1L
      message(
        "[", counter, "/", total, "] ",
        contrast$meta$control_group, ": ",
        contrast$meta$cat1, " ~ ", contrast$meta$cat2,
        " / ", feature_set$name
      )
      result <- run_one_contrast(
        df = contrast$data,
        feature_set = feature_set,
        category_col = args$category_col,
        meta = contrast$meta,
        min_per_category = args$min_per_category,
        min_tokens = args$min_tokens,
        bw = args$bw,
        eval_on = args$eval_on,
        eval_n = args$eval_n,
        eval_seed = args$eval_seed,
        engine = args$engine,
        chunk_size = args$chunk_size,
        eps = args$eps,
        cv_folds = args$cv_folds,
        seed = args$seed + counter
      )
      if (!is.null(result$metrics)) {
        metric_rows[[length(metric_rows) + 1L]] <- result$metrics
      }
      if (!is.null(result$skipped)) {
        skipped_rows[[length(skipped_rows) + 1L]] <- result$skipped
      }
    }
  }

  metrics <- if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame()
  skipped <- if (length(skipped_rows)) do.call(rbind, skipped_rows) else empty_skipped()
  summary <- data.frame(
    input = args$input,
    domain = args$domain %||% "all",
    n_input_rows = nrow(df),
    n_feature_sets = length(feature_sets),
    n_contrasts = length(contrasts),
    n_attempted = total,
    n_successful = nrow(metrics),
    n_skipped = nrow(skipped),
    bw = args$bw,
    eval_on = args$eval_on,
    eval_n = args$eval_n,
    engine = args$engine,
    min_per_category = args$min_per_category,
    min_tokens = args$min_tokens,
    cv_folds = args$cv_folds,
    seed = args$seed,
    stringsAsFactors = FALSE
  )

  paths <- write_outputs(metrics, skipped, args$out_dir, summary)
  message("Wrote: ", paths$metrics_path)
  message("Wrote: ", paths$skipped_path)
  message("Wrote: ", paths$summary_path)
  invisible(0L)
}

if (identical(environment(), globalenv())) {
  main()
}
