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
