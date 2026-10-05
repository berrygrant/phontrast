# ---- The kernel family, scored on one shared density estimate ----------------
#
# phontrast()'s distributional metrics -- Jensen-Shannon divergence and
# distance, proportional overlap and its complement total variation, and the
# matched-kernel Bhattacharyya affinity / distance and Hellinger distance --
# are all read off a single estimate of the two category densities per
# comparison. Scoring every member on the same densities is what makes them
# comparable as measures rather than as estimators; it also halves the kernel
# work of a bootstrap, which recomputes this block on every resample.

.kernel_family_columns <- function() {
  c("jsd", "js_distance", "percent_overlap", "total_variation",
    "bhatt_kde_dist", "bhatt_kde_affinity", "hellinger")
}

.kernel_family_vector <- function(jsd, overlap, bc) {
  c(
    jsd = jsd,
    js_distance = sqrt(jsd),
    percent_overlap = overlap,
    total_variation = 1 - overlap,
    bhatt_kde_dist = -log(bc),
    bhatt_kde_affinity = bc,
    hellinger = sqrt(max(1 - bc, 0))
  )
}

.kernel_family_from_mc <- function(mc, loo) {
  .kernel_family_vector(
    jsd = .jsd_mc(mc, loo = loo),
    overlap = .overlap_mc(mc),
    bc = .bhatt_mc(mc, loo = loo)
  )
}

.kernel_family_point <- function(df, features, category_col, bw, eval_on,
                                 eval_n, eval_seed, engine, chunk_size, method,
                                 density, mc_n, bw_scale, loo = TRUE,
                                 metric = "phontrast()") {
  if (identical(density, "mvnorm")) {
    mc <- .mvnorm_mc_pair(
      data = df, features = features, category_col = category_col,
      mc_n = mc_n, eval_seed = eval_seed, metric = metric
    )
    return(.kernel_family_from_mc(mc, loo = FALSE))
  }
  if (identical(method, "mc")) {
    mc <- .kde_mc_pair(
      data = df, features = features, category_col = category_col, bw = bw,
      eval_n = eval_n, eval_seed = eval_seed, engine = engine,
      chunk_size = chunk_size, metric = metric, bw_scale = bw_scale
    )
    return(.kernel_family_from_mc(mc, loo = loo))
  }
  # Legacy self-normalized sample-point estimates on a shared evaluation set.
  dens <- .kde_density_pair(
    data = df, features = features, category_col = category_col, bw = bw,
    eval_on = eval_on, eval_n = eval_n, eval_seed = eval_seed, engine = engine,
    chunk_size = chunk_size, metric = metric, bw_scale = bw_scale
  )
  p <- dens$p / sum(dens$p)
  q <- dens$q / sum(dens$q)
  .kernel_family_vector(
    jsd = jsd(dens$p, dens$q),
    overlap = min(max(sum(pmin(p, q)), 0), 1),
    bc = min(max(sum(sqrt(p * q)), 0), 1)
  )
}

.estimate_kernel_family <- function(data, features, category_col,
                                    group_col = NULL, min_tokens = 20, bw,
                                    eval_on, eval_n, eval_seed, engine,
                                    chunk_size, method, density, mc_n,
                                    bw_scale, loo = TRUE) {
  cols <- .kernel_family_columns()
  point <- function(df) {
    .kernel_family_point(
      df = df, features = features, category_col = category_col, bw = bw,
      eval_on = eval_on, eval_n = eval_n, eval_seed = eval_seed,
      engine = engine, chunk_size = chunk_size, method = method,
      density = density, mc_n = mc_n, bw_scale = bw_scale, loo = loo
    )
  }

  if (is.null(group_col)) {
    df <- .metric_data(data, c(category_col, features))
    n <- nrow(df)
    if (n < min_tokens) {
      stop(
        "phontrast(): Not enough tokens after removing missing values. Got ",
        n, ", need at least ", min_tokens, ".",
        call. = FALSE
      )
    }
    .two_levels(df[[category_col]], "category_col")
    out <- data.frame(scope = "global", n_tokens = n, stringsAsFactors = FALSE)
    out[cols] <- as.list(point(df))
    return(out)
  }

  group_col <- .check_group_cols(group_col)
  .check_columns(data, c(group_col, category_col, features))
  df <- .metric_data(data, c(group_col, category_col, features))
  out <- lapply(.split_groups(df, group_col), function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens || .observed_n_categories(df_g[[category_col]]) != 2L) {
      return(NULL)
    }
    vals <- tryCatch(
      point(df_g),
      error = function(e) stats::setNames(rep(NA_real_, length(cols)), cols)
    )
    row <- data.frame(
      scope = "group",
      group = .group_label(df_g, group_col),
      n_tokens = n_tok,
      stringsAsFactors = FALSE
    )
    row[cols] <- as.list(vals)
    row
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(
      scope = character(), group = character(), n_tokens = integer(),
      stringsAsFactors = FALSE
    )
    for (col in cols) out[[col]] <- numeric()
  }
  rownames(out) <- NULL
  .warn_failed_groups(out, "jsd", "phontrast()")
}

# ---- Euclidean distance between standardized category means -----------------

.euclidean_means <- function(data, features, category_col) {
  # Distance between the two category centroids after dividing each feature
  # by the standard deviation of the pooled two-category sample. Standardizing
  # on the pooled sample (not the within-category spread) bounds each feature's
  # standardized mean difference by 2 for equal category sizes, so the
  # distance lies in [0, 2 * sqrt(d)).
  .check_columns(data, c(category_col, features))
  data <- .metric_data(data, c(category_col, features))
  .check_numeric_features(data, features)
  levs <- .two_levels(data[[category_col]], "category_col")
  .check_two_category_sample_size(
    data, category_col, 2L, "Euclidean distance of standardized means"
  )

  X <- as.matrix(data[, features, drop = FALSE])
  s <- apply(X, 2, stats::sd)
  if (any(!is.finite(s)) || any(s <= 0)) {
    stop(
      "Euclidean distance of standardized means failed: a feature has zero or ",
      "non-finite spread across the pooled sample.",
      call. = FALSE
    )
  }
  mu1 <- colMeans(X[data[[category_col]] == levs[1], , drop = FALSE])
  mu2 <- colMeans(X[data[[category_col]] == levs[2], , drop = FALSE])
  sqrt(sum(((mu2 - mu1) / s)^2))
}

.estimate_euclidean <- function(data, features, category_col, group_col = NULL,
                                min_tokens = 20) {
  .check_positive_count(min_tokens, "min_tokens")

  if (is.null(group_col)) {
    df <- .metric_data(data, c(category_col, features))
    n <- nrow(df)
    if (n < min_tokens) {
      stop("Not enough tokens after removing missing values. Got ",
           n, ", need at least ", min_tokens, ".")
    }
    return(data.frame(
      scope = "global",
      n_tokens = n,
      euclidean_dist = .euclidean_means(df, features, category_col),
      stringsAsFactors = FALSE
    ))
  }

  group_col <- .check_group_cols(group_col)
  .check_columns(data, c(group_col, category_col, features))
  df <- .metric_data(data, c(group_col, category_col, features))
  out <- lapply(.split_groups(df, group_col), function(df_g) {
    n_tok <- nrow(df_g)
    if (n_tok < min_tokens || .observed_n_categories(df_g[[category_col]]) != 2L) {
      return(NULL)
    }
    dist <- tryCatch(
      .euclidean_means(df_g, features, category_col),
      error = function(e) NA_real_
    )
    data.frame(
      scope = "group",
      group = .group_label(df_g, group_col),
      n_tokens = n_tok,
      euclidean_dist = dist,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(
      scope = character(),
      group = character(),
      n_tokens = integer(),
      euclidean_dist = numeric(),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  .warn_failed_groups(out, "euclidean_dist", "estimate_euclidean()")
}
