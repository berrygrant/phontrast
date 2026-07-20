.check_columns <- function(data, cols, arg = "columns") {
  missing <- setdiff(cols, names(data))
  if (length(missing)) {
    stop(
      "`", arg, "` must exist in `data`: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

.metric_data <- function(data, cols) {
  .check_columns(data, cols)
  out <- data[stats::complete.cases(data[, cols, drop = FALSE]), cols, drop = FALSE]
  numeric_cols <- cols[vapply(out[, cols, drop = FALSE], is.numeric, logical(1))]
  if (length(numeric_cols)) {
    finite_rows <- Reduce(
      `&`,
      lapply(out[, numeric_cols, drop = FALSE], is.finite)
    )
    out <- out[finite_rows, , drop = FALSE]
  }
  out
}

.observed_category_counts <- function(x) {
  x <- x[!is.na(x)]
  if (is.factor(x)) {
    x <- droplevels(x)
  }
  table(x)
}

.observed_n_categories <- function(x) {
  length(.observed_category_counts(x))
}

.two_levels <- function(x, arg = "category") {
  x <- x[!is.na(x)]
  if (is.factor(x)) {
    x <- as.character(droplevels(x))
  }
  levs <- unique(x)
  if (length(levs) != 2L) {
    stop("`", arg, "` must have exactly two non-missing values.", call. = FALSE)
  }
  levs
}

.check_numeric_features <- function(data, features) {
  if (!length(features)) {
    stop("`features` must contain at least one column.", call. = FALSE)
  }
  is_num <- vapply(data[, features, drop = FALSE], is.numeric, logical(1))
  if (!all(is_num)) {
    stop(
      "All `features` must be numeric. Non-numeric feature(s): ",
      paste(features[!is_num], collapse = ", "),
      call. = FALSE
    )
  }
}

.check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single finite number between 0 and 1.", call. = FALSE)
  }
}

.check_positive_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != as.integer(x)) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }
}

.check_group_cols <- function(group_col, arg = "group_col") {
  if (is.null(group_col)) {
    return(NULL)
  }
  if (!is.character(group_col) || !length(group_col) ||
      anyNA(group_col) || any(!nzchar(group_col))) {
    stop("`", arg, "` must be a character vector of one or more column names.", call. = FALSE)
  }
  if (anyDuplicated(group_col)) {
    stop("`", arg, "` must not contain duplicate column names.", call. = FALSE)
  }
  group_col
}

.validate_metric_inputs <- function(data, features, category_col, group_col = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(category_col) || length(category_col) != 1L ||
      is.na(category_col) || !nzchar(category_col)) {
    stop("`category_col` must be a single column name.", call. = FALSE)
  }
  if (!is.character(features) || !length(features) || anyNA(features) ||
      any(!nzchar(features))) {
    stop("`features` must be a non-empty character vector of column names.", call. = FALSE)
  }
  group_col <- .check_group_cols(group_col)
  clash <- intersect(features, c(category_col, group_col))
  if (length(clash)) {
    stop(
      "`features` must not overlap with `category_col`/`group_col`: ",
      paste(clash, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(group_col) && category_col %in% group_col) {
    stop("`category_col` must not also appear in `group_col`.", call. = FALSE)
  }
  invisible(TRUE)
}

.warn_failed_groups <- function(out, value_col, fn) {
  if (nrow(out) > 0L && value_col %in% names(out)) {
    failed <- sum(is.na(out[[value_col]]))
    if (failed) {
      warning(
        fn, ": ", failed, " of ", nrow(out),
        " group(s) could not be estimated and were returned as NA.",
        call. = FALSE
      )
    }
  }
  out
}

.group_label <- function(data, group_col) {
  group_col <- .check_group_cols(group_col)
  if (length(group_col) == 1L) {
    return(data[[group_col]][1])
  }

  vals <- vapply(group_col, function(col) {
    as.character(data[[col]][1])
  }, character(1))
  paste(paste0(group_col, "=", vals), collapse = " | ")
}

.match_kde_engine <- function(engine) {
  engine <- match.arg(engine, c("ks", "fast_diag", "fast_diagonal"))
  if (identical(engine, "fast_diagonal")) {
    return("fast_diag")
  }
  engine
}

.check_ridge_eps <- function(eps, arg = "eps") {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`", arg, "` must be a single positive finite number.", call. = FALSE)
  }
}

.empty_group_pillai <- function() {
  tibble::tibble(
    group = character(),
    n_tokens = integer(),
    pillai = numeric(),
    p_value = numeric()
  )
}

.empty_group_jsd <- function() {
  tibble::tibble(group = character(), n_tokens = integer(), jsd = numeric())
}

.empty_group_bhatt <- function() {
  data.frame(
    group = character(),
    n_tokens = integer(),
    bhatt_dist = numeric(),
    bhatt_affinity = numeric(),
    stringsAsFactors = FALSE
  )
}

.empty_estimate_bhatt_group <- function() {
  data.frame(
    scope = character(),
    group = character(),
    n_tokens = integer(),
    bhatt_dist = numeric(),
    bhatt_affinity = numeric(),
    stringsAsFactors = FALSE
  )
}

.kde_min_category_tokens <- function(n_features) {
  max(2L, n_features + 1L)
}

.check_two_category_sample_size <- function(data, category_col, min_per_category, metric) {
  counts <- .observed_category_counts(data[[category_col]])
  if (length(counts) != 2L || any(counts < min_per_category)) {
    stop(
      metric, " requires at least ", min_per_category,
      " finite observations in each category after removing missing values.",
      call. = FALSE
    )
  }
}

.split_groups <- function(data, group_col) {
  group_col <- .check_group_cols(group_col)
  if (length(group_col) == 1L) {
    return(split(data, data[[group_col]], drop = TRUE))
  }

  group_key <- interaction(data[, group_col, drop = FALSE], drop = TRUE, sep = "\r")
  groups <- split(data, group_key, drop = TRUE)
  names(groups) <- vapply(groups, .group_label, character(1), group_col = group_col)
  groups
}

.select_univariate_bandwidth <- function(x, bw) {
  if (identical(bw, "scott.diag")) {
    h <- stats::sd(x) * length(x) ^ (-1 / 5)
    if (!is.finite(h) || h <= 0) {
      spread <- max(stats::sd(x), diff(range(x)), 1, na.rm = TRUE)
      h <- spread / 10
    }
    return(h)
  }

  selector <- switch(
    bw,
    Hpi = stats::bw.SJ,
    Hscv = stats::bw.ucv,
    Hpi.diag = stats::bw.nrd0
  )

  h <- tryCatch(selector(x), error = function(e) NA_real_)
  if (!is.finite(h) || h <= 0) {
    h <- tryCatch(stats::density(x)$bw, error = function(e) NA_real_)
  }
  if (!is.finite(h) || h <= 0) {
    spread <- max(stats::sd(x), diff(range(x)), 1, na.rm = TRUE)
    h <- spread / 10
  }
  h
}

.scott_diag_bandwidth <- function(x) {
  d <- ncol(x)
  n <- nrow(x)
  sds <- apply(x, 2, stats::sd)
  bad <- !is.finite(sds) | sds <= 0
  if (any(bad)) {
    spreads <- apply(x[, bad, drop = FALSE], 2, function(col) {
      max(stats::sd(col), diff(range(col)), 1, na.rm = TRUE)
    })
    sds[bad] <- spreads
  }
  h <- n ^ (-1 / (d + 4))
  diag((h * sds) ^ 2, nrow = d, ncol = d)
}

.kde_1d_values <- function(x, eval_points, h, chunk_size = 5000L) {
  out <- numeric(length(eval_points))
  starts <- seq.int(1L, length(eval_points), by = chunk_size)
  for (start in starts) {
    stop <- min(start + chunk_size - 1L, length(eval_points))
    z <- outer(eval_points[start:stop], x, `-`) / h
    out[start:stop] <- rowMeans(stats::dnorm(z)) / h
  }
  out
}

.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

.kde_diag_gaussian_values <- function(x, eval_points, H, chunk_size = 1000L) {
  if (!is.matrix(H) || nrow(H) != ncol(H) || nrow(H) != ncol(x)) {
    stop("Diagonal KDE engine received an invalid bandwidth matrix.", call. = FALSE)
  }
  off_diag <- H
  diag(off_diag) <- 0
  if (any(abs(off_diag) > sqrt(.Machine$double.eps))) {
    stop(
      "`engine = \"fast_diag\"` requires a diagonal bandwidth matrix. ",
      "Use `bw = \"scott.diag\"` or `bw = \"Hpi.diag\"`.",
      call. = FALSE
    )
  }

  variances <- diag(H)
  if (any(!is.finite(variances)) || any(variances <= 0)) {
    stop("KDE bandwidth matrix must have positive finite diagonal entries.", call. = FALSE)
  }

  x <- as.matrix(x)
  eval_points <- as.matrix(eval_points)
  log_density <- numeric(nrow(eval_points))
  inv_variances <- 1 / variances
  starts <- seq.int(1L, nrow(eval_points), by = chunk_size)

  for (start in starts) {
    stop <- min(start + chunk_size - 1L, nrow(eval_points))
    eval_chunk <- eval_points[start:stop, , drop = FALSE]
    log_kernel <- matrix(0, nrow = nrow(eval_chunk), ncol = nrow(x))
    for (j in seq_len(ncol(x))) {
      diff <- outer(eval_chunk[, j], x[, j], `-`)
      log_kernel <- log_kernel - 0.5 * diff * diff * inv_variances[j]
    }
    log_density[start:stop] <- apply(log_kernel, 1, .logsumexp) - log(nrow(x))
  }

  scale <- max(log_density)
  if (!is.finite(scale)) {
    return(rep(0, length(log_density)))
  }
  exp(log_density - scale)
}

.sample_kde_eval_points <- function(eval_pts, eval_n = NULL, eval_seed = NULL) {
  if (is.null(eval_n)) {
    return(eval_pts)
  }
  .check_positive_count(eval_n, "eval_n")
  eval_n <- as.integer(eval_n)
  if (nrow(eval_pts) <= eval_n) {
    return(eval_pts)
  }

  if (!is.null(eval_seed)) {
    if (!is.numeric(eval_seed) || length(eval_seed) != 1L ||
        !is.finite(eval_seed) || eval_seed != as.integer(eval_seed)) {
      stop("`eval_seed` must be NULL or a single finite integer.", call. = FALSE)
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(eval_seed))
  }

  idx <- sample.int(nrow(eval_pts), eval_n)
  eval_pts[idx, , drop = FALSE]
}

.select_multivariate_bandwidth <- function(x, bw, label) {
  tryCatch(
    switch(
      bw,
      Hpi = ks::Hpi(x),
      Hscv = ks::Hscv(x),
      Hpi.diag = ks::Hpi.diag(x),
      scott.diag = .scott_diag_bandwidth(x)
    ),
    error = function(e) {
      stop(
        "KDE bandwidth selection failed for category `", label, "`. ",
        "Check that the category has more observations than feature dimensions ",
        "and that feature columns are not constant or collinear. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

.kde_density_pair <- function(data,
                              features,
                              category_col,
                              bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                              eval_on = c("pooled", "group1", "group2", "pooled_sample"),
                              eval_n = NULL,
                              eval_seed = NULL,
                              engine = c("ks", "fast_diag", "fast_diagonal"),
                              chunk_size = 1000L,
                              metric = "KDE") {
  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)
  engine <- .match_kde_engine(engine)
  if (!is.null(eval_n)) {
    .check_positive_count(eval_n, "eval_n")
  }
  if (identical(eval_on, "pooled_sample") && is.null(eval_n)) {
    stop("`eval_n` must be supplied when `eval_on = \"pooled_sample\"`.", call. = FALSE)
  }
  .check_positive_count(chunk_size, "chunk_size")

  .check_columns(data, c(category_col, features))
  data <- .metric_data(data, c(category_col, features))
  .check_numeric_features(data, features)

  levs <- .two_levels(data[[category_col]], "category_col")
  n_features <- length(features)
  .check_two_category_sample_size(
    data,
    category_col,
    .kde_min_category_tokens(n_features),
    metric
  )

  d1 <- data[data[[category_col]] == levs[1], , drop = FALSE]
  d2 <- data[data[[category_col]] == levs[2], , drop = FALSE]

  X1 <- as.matrix(d1[, features, drop = FALSE])
  X2 <- as.matrix(d2[, features, drop = FALSE])
  X_all <- as.matrix(data[, features, drop = FALSE])

  eval_source <- if (identical(eval_on, "pooled_sample")) "pooled" else eval_on
  eval_pts <- switch(
    eval_source,
    pooled = X_all,
    group1 = X1,
    group2 = X2
  )
  eval_pts <- .sample_kde_eval_points(eval_pts, eval_n = eval_n, eval_seed = eval_seed)

  if (n_features == 1L) {
    x1 <- as.numeric(X1[, 1])
    x2 <- as.numeric(X2[, 1])
    eval_vec <- as.numeric(eval_pts[, 1])
    h1 <- .select_univariate_bandwidth(x1, bw)
    h2 <- .select_univariate_bandwidth(x2, bw)
    p <- .kde_1d_values(x1, eval_vec, h1)
    q <- .kde_1d_values(x2, eval_vec, h2)
  } else {
    if (identical(engine, "fast_diag") &&
        !bw %in% c("Hpi.diag", "scott.diag")) {
      stop(
        "`engine = \"fast_diag\"` requires `bw = \"scott.diag\"` or ",
        "`bw = \"Hpi.diag\"` for multivariate KDE.",
        call. = FALSE
      )
    }

    H1 <- .select_multivariate_bandwidth(X1, bw, levs[1])
    H2 <- .select_multivariate_bandwidth(X2, bw, levs[2])

    if (identical(engine, "fast_diag")) {
      p <- .kde_diag_gaussian_values(X1, eval_pts, H1, chunk_size = chunk_size)
      q <- .kde_diag_gaussian_values(X2, eval_pts, H2, chunk_size = chunk_size)
    } else {
      kde1 <- tryCatch(
        ks::kde(x = X1, H = H1, eval.points = eval_pts),
        error = function(e) {
          stop(
            "KDE estimation failed for category `", levs[1], "`. Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      kde2 <- tryCatch(
        ks::kde(x = X2, H = H2, eval.points = eval_pts),
        error = function(e) {
          stop(
            "KDE estimation failed for category `", levs[2], "`. Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      p <- as.numeric(kde1$estimate)
      q <- as.numeric(kde2$estimate)
    }
  }

  if (any(!is.finite(p)) || any(!is.finite(q)) || sum(p) <= 0 || sum(q) <= 0) {
    stop("KDE returned invalid density estimates.", call. = FALSE)
  }

  list(p = p, q = q, levels = levs, data = data)
}

# ---- Monte-Carlo plug-in KDE estimator ----------------------------------
# Consistent estimator of the continuous Jensen-Shannon divergence / overlap:
# evaluate each category's KDE at that category's own observations and average
# the true log density ratio against the mixture. Dimension-agnostic (unlike a
# grid) and unbiased in the limit (unlike the self-normalized sample-point plug
# -in used by `.kde_density_pair()`, kept as the `method = "legacy"` path).

.log_add_exp <- function(a, b) {
  n <- max(length(a), length(b))
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)
  m <- pmax(a, b)
  out <- m
  # Only evaluate log1p where the max is finite; both -Inf stays -Inf. Indexing
  # (not ifelse) avoids computing NaN intermediates that emit spurious warnings.
  fin <- is.finite(m)
  out[fin] <- m[fin] + log1p(exp(-abs(a[fin] - b[fin])))
  out
}

.log_sub_exp <- function(a, b) {
  # log(exp(a) - exp(b)); -Inf where a <= b (e.g., an isolated leave-one-out
  # point). Evaluate the log only on strictly-greater elements so `log1p` is
  # never handed a value <= -1 (which would emit a "NaNs produced" warning).
  n <- max(length(a), length(b))
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)
  out <- rep(-Inf, n)
  ok <- a > b
  ok[is.na(ok)] <- FALSE
  out[ok] <- a[ok] + log1p(-exp(b[ok] - a[ok]))
  out
}

.kde_diag_log_density <- function(x, eval_points, H, chunk_size = 1000L) {
  variances <- diag(H)
  d <- ncol(x)
  log_norm <- -0.5 * (d * log(2 * pi) + sum(log(variances)))
  x <- as.matrix(x)
  eval_points <- as.matrix(eval_points)
  log_density <- numeric(nrow(eval_points))
  inv_variances <- 1 / variances
  starts <- seq.int(1L, nrow(eval_points), by = chunk_size)
  for (start in starts) {
    stop <- min(start + chunk_size - 1L, nrow(eval_points))
    eval_chunk <- eval_points[start:stop, , drop = FALSE]
    log_kernel <- matrix(0, nrow = nrow(eval_chunk), ncol = nrow(x))
    for (j in seq_len(ncol(x))) {
      diff <- outer(eval_chunk[, j], x[, j], `-`)
      log_kernel <- log_kernel - 0.5 * diff * diff * inv_variances[j]
    }
    log_density[start:stop] <- apply(log_kernel, 1, .logsumexp) - log(nrow(x))
  }
  log_density + log_norm
}

.kde_kh0 <- function(bwspec, d) {
  # K_H(0): the kernel's value at the origin (self-contribution), for LOO.
  if (is.matrix(bwspec)) {
    (2 * pi) ^ (-d / 2) * det(bwspec) ^ (-0.5)
  } else {
    stats::dnorm(0) / bwspec
  }
}

.kde_eval_logdens <- function(train, eval, bwspec, engine, chunk_size, label) {
  if (!is.matrix(bwspec)) {
    dens <- .kde_1d_values(as.numeric(train[, 1]), as.numeric(eval[, 1]), bwspec)
    return(log(dens))
  }
  if (identical(engine, "fast_diag")) {
    return(.kde_diag_log_density(train, eval, bwspec, chunk_size = chunk_size))
  }
  kde <- tryCatch(
    ks::kde(x = train, H = bwspec, eval.points = eval),
    error = function(e) {
      stop(
        "KDE estimation failed for category `", label, "`. Original error: ",
        conditionMessage(e), call. = FALSE
      )
    }
  )
  dens <- as.numeric(kde$estimate)
  dens[dens < 0] <- 0
  log(dens)
}

.select_kde_bandwidth <- function(train, bw, engine, n_features, label) {
  if (n_features == 1L) {
    return(.select_univariate_bandwidth(as.numeric(train[, 1]), bw))
  }
  if (identical(engine, "fast_diag") && !bw %in% c("Hpi.diag", "scott.diag")) {
    stop(
      "`engine = \"fast_diag\"` requires `bw = \"scott.diag\"` or ",
      "`bw = \"Hpi.diag\"` for multivariate KDE.",
      call. = FALSE
    )
  }
  .select_multivariate_bandwidth(train, bw, label)
}

.kde_mc_pair <- function(data,
                         features,
                         category_col,
                         bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                         eval_n = NULL,
                         eval_seed = NULL,
                         engine = c("ks", "fast_diag", "fast_diagonal"),
                         chunk_size = 1000L,
                         metric = "KDE") {
  bw <- match.arg(bw)
  engine <- .match_kde_engine(engine)
  if (!is.null(eval_n)) {
    .check_positive_count(eval_n, "eval_n")
  }
  .check_positive_count(chunk_size, "chunk_size")

  .check_columns(data, c(category_col, features))
  data <- .metric_data(data, c(category_col, features))
  .check_numeric_features(data, features)

  levs <- .two_levels(data[[category_col]], "category_col")
  n_features <- length(features)
  .check_two_category_sample_size(
    data, category_col, .kde_min_category_tokens(n_features), metric
  )

  X1 <- as.matrix(data[data[[category_col]] == levs[1], features, drop = FALSE])
  X2 <- as.matrix(data[data[[category_col]] == levs[2], features, drop = FALSE])
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  # KDEs are trained on the full samples; evaluation points may be subsampled
  # for speed (leave-one-out below still uses the full training size n1/n2).
  X1e <- .sample_kde_eval_points(X1, eval_n = eval_n, eval_seed = eval_seed)
  X2e <- .sample_kde_eval_points(X2, eval_n = eval_n, eval_seed = eval_seed)

  bw1 <- .select_kde_bandwidth(X1, bw, engine, n_features, levs[1])
  bw2 <- .select_kde_bandwidth(X2, bw, engine, n_features, levs[2])

  out <- list(
    logp1 = .kde_eval_logdens(X1, X1e, bw1, engine, chunk_size, levs[1]),
    logq1 = .kde_eval_logdens(X2, X1e, bw2, engine, chunk_size, levs[2]),
    logp2 = .kde_eval_logdens(X1, X2e, bw1, engine, chunk_size, levs[1]),
    logq2 = .kde_eval_logdens(X2, X2e, bw2, engine, chunk_size, levs[2]),
    n1 = n1, n2 = n2,
    kh0_1 = .kde_kh0(bw1, n_features),
    kh0_2 = .kde_kh0(bw2, n_features),
    levels = levs, data = data
  )
  if (any(!is.finite(out$logp1)) && any(!is.finite(out$logp2))) {
    stop("KDE returned invalid density estimates.", call. = FALSE)
  }
  out
}

.loo_logdens <- function(log_dens, n, kh0) {
  # leave-one-out log density at a KDE's own training points
  .log_sub_exp(log(n) + log_dens, log(kh0)) - log(n - 1)
}

.jsd_mc <- function(mc, loo = TRUE) {
  ln2 <- log(2)
  logp1 <- if (isTRUE(loo)) .loo_logdens(mc$logp1, mc$n1, mc$kh0_1) else mc$logp1
  logm1 <- log(0.5) + .log_add_exp(logp1, mc$logq1)
  t1 <- (logp1 - logm1) / ln2

  logq2 <- if (isTRUE(loo)) .loo_logdens(mc$logq2, mc$n2, mc$kh0_2) else mc$logq2
  logm2 <- log(0.5) + .log_add_exp(mc$logp2, logq2)
  t2 <- (logq2 - logm2) / ln2

  t1 <- t1[is.finite(t1)]
  t2 <- t2[is.finite(t2)]
  if (!length(t1) || !length(t2)) {
    stop("Monte-Carlo JSD: no usable evaluation points.", call. = FALSE)
  }
  min(max(0.5 * mean(t1) + 0.5 * mean(t2), 0), 1)
}

.overlap_mc <- function(mc) {
  # OVL = integral of min(p, q); estimate each half with that group's own
  # samples via min(1, cross-density / self-density).
  o1 <- pmin(1, exp(mc$logq1 - mc$logp1))
  o2 <- pmin(1, exp(mc$logp2 - mc$logq2))
  o1 <- o1[is.finite(o1)]
  o2 <- o2[is.finite(o2)]
  if (!length(o1) || !length(o2)) {
    stop("Monte-Carlo overlap: no usable evaluation points.", call. = FALSE)
  }
  min(max(0.5 * mean(o1) + 0.5 * mean(o2), 0), 1)
}
