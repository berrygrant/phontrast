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

.two_levels <- function(x, arg = "category") {
  levs <- unique(x)
  levs <- levs[!is.na(levs)]
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
  counts <- table(data[[category_col]])
  if (length(counts) != 2L || any(counts < min_per_category)) {
    stop(
      metric, " requires at least ", min_per_category,
      " finite observations in each category after removing missing values.",
      call. = FALSE
    )
  }
}

.select_univariate_bandwidth <- function(x, bw) {
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

.select_multivariate_bandwidth <- function(x, bw, label) {
  tryCatch(
    switch(
      bw,
      Hpi = ks::Hpi(x),
      Hscv = ks::Hscv(x),
      Hpi.diag = ks::Hpi.diag(x)
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
                              bw = c("Hpi", "Hscv", "Hpi.diag"),
                              eval_on = c("pooled", "group1", "group2"),
                              metric = "KDE") {
  bw <- match.arg(bw)
  eval_on <- match.arg(eval_on)

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

  eval_pts <- switch(
    eval_on,
    pooled = X_all,
    group1 = X1,
    group2 = X2
  )

  if (n_features == 1L) {
    x1 <- as.numeric(X1[, 1])
    x2 <- as.numeric(X2[, 1])
    eval_vec <- as.numeric(eval_pts[, 1])
    h1 <- .select_univariate_bandwidth(x1, bw)
    h2 <- .select_univariate_bandwidth(x2, bw)
    p <- .kde_1d_values(x1, eval_vec, h1)
    q <- .kde_1d_values(x2, eval_vec, h2)
  } else {
    H1 <- .select_multivariate_bandwidth(X1, bw, levs[1])
    H2 <- .select_multivariate_bandwidth(X2, bw, levs[2])

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

  if (any(!is.finite(p)) || any(!is.finite(q)) || sum(p) <= 0 || sum(q) <= 0) {
    stop("KDE returned invalid density estimates.", call. = FALSE)
  }

  list(p = p, q = q, levels = levs, data = data)
}
