#' Distribution-aware contrast plot for two categories
#'
#' The flagship visualization for a two-category contrast. Unlike a decorative
#' scatterplot, \code{plot_contrast()} draws the \emph{same density model the
#' distributional metrics are computed from}: under \code{density = "kde"} it
#' shows highest-density regions of each category's kernel density estimate
#' (same bandwidth selection as \code{jsd_kde_nd()}); under
#' \code{density = "mvnorm"} it shows coverage ellipses of the fitted
#' multivariate normals used by the parametric backend. The pointwise minimum
#' of the two densities -- the mass that the proportional-overlap metric
#' integrates -- is shaded, so the overlap itself is visible rather than
#' implied.
#'
#' With \code{annotate = TRUE} (default) the panel is labelled with the
#' Jensen-Shannon divergence and proportional overlap computed by
#' \code{phontrast()} under the same \code{density}, \code{bw}, \code{mc_n},
#' and \code{eval_seed} settings, and the caption records the estimator
#' configuration. The full annotation table is attached to the returned plot
#' as \code{attr(p, "contrast_metrics")}.
#'
#' @param data Data frame with category labels and one or two numeric features.
#' @param features One or two numeric feature columns. One feature gives
#'   density curves; two give a feature-space plot with density regions. For
#'   higher-dimensional spaces, plot a projection with
#'   \code{plot_category_pca()} (metrics should still be computed on the full
#'   space).
#' @param category_col String; category column with exactly two observed
#'   categories.
#' @param group_col Optional character vector of grouping columns; one panel
#'   per group, with per-group densities and annotations.
#' @param density Density model to draw and to use for annotations:
#'   \code{"kde"} (default) or \code{"mvnorm"}. Matches the \code{density}
#'   argument of the metric functions.
#' @param bw Bandwidth selection method for \code{density = "kde"}; same
#'   options as \code{jsd_kde_nd()}.
#' @param bw_scale Positive multiplier on the selected kernel bandwidth for
#'   \code{density = "kde"} (default \code{1}; see \code{jsd_kde_nd()}). The
#'   drawn regions, the shaded overlap, and the annotations all use the same
#'   scaled bandwidth, so \code{bw_scale = 0.5} and \code{2} show the halved
#'   and doubled smoothing of the bandwidth check in \code{rank_contrasts()}.
#' @param levels Numeric vector of probability levels in (0, 1) for the drawn
#'   regions: highest-density regions under \code{"kde"}, coverage ellipses
#'   under \code{"mvnorm"}.
#' @param points Logical; show observed tokens (2D points, 1D rug).
#' @param overlap Logical; shade the pointwise minimum of the two category
#'   densities (a ribbon in 1D, a soft raster in 2D). Shading strength is
#'   normalized across panels, so lighter panels genuinely overlap less.
#' @param annotate Logical; label each panel with Jensen-Shannon divergence
#'   and proportional overlap computed under the plotted density model.
#' @param n_boot Number of bootstrap resamples for annotation confidence
#'   intervals; \code{0} (default) annotates point estimates only.
#' @param conf_level Confidence level for bootstrap intervals.
#' @param min_tokens Minimum tokens per group; smaller groups are dropped with
#'   a warning (same convention as the metric functions).
#' @param mc_n Monte-Carlo sample size for \code{density = "mvnorm"}
#'   annotations; passed to the metric functions.
#' @param eval_seed Optional integer seed passed to the metric functions so
#'   annotated values are reproducible.
#' @param grid_n Grid resolution for density evaluation: points per axis.
#'   Default 512 for one feature, 151 for two.
#' @param point_alpha Point (or rug) transparency.
#' @param point_size Point size for two-feature plots.
#' @param reverse_x,reverse_y Logical; reverse an axis (e.g. F2 by F1 vowel
#'   space convention).
#' @param facet_scales Scales passed to \code{ggplot2::facet_wrap()} when
#'   \code{group_col} is supplied.
#'
#' @return A \pkg{ggplot2} plot object. When \code{annotate = TRUE}, the
#'   \code{phontrast()} table behind the labels is attached as
#'   \code{attr(p, "contrast_metrics")}.
#'
#' @examples
#' set.seed(2026)
#' vowels <- data.frame(
#'   vowel = rep(c("ih", "eh"), each = 60),
#'   f1 = c(rnorm(60, 500, 55), rnorm(60, 565, 60)),
#'   f2 = c(rnorm(60, 1980, 150), rnorm(60, 1870, 155))
#' )
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Two-feature contrast in vowel-space orientation, KDE regions.
#'   plot_contrast(vowels, c("f2", "f1"), "vowel",
#'                 reverse_x = TRUE, reverse_y = TRUE)
#'
#'   # One-feature contrast with the overlap ribbon.
#'   plot_contrast(vowels, "f1", "vowel")
#'
#'   # The same contrast under the multivariate-normal backend.
#'   plot_contrast(vowels, c("f2", "f1"), "vowel", density = "mvnorm",
#'                 eval_seed = 2026, reverse_x = TRUE, reverse_y = TRUE)
#' }
#' @export
plot_contrast <- function(data,
                          features,
                          category_col,
                          group_col = NULL,
                          density = c("kde", "mvnorm"),
                          bw = c("Hpi", "Hscv", "Hpi.diag", "scott.diag"),
                          bw_scale = 1,
                          levels = c(0.5, 0.8, 0.95),
                          points = TRUE,
                          overlap = TRUE,
                          annotate = TRUE,
                          n_boot = 0,
                          conf_level = 0.95,
                          min_tokens = 20,
                          mc_n = 10000L,
                          eval_seed = NULL,
                          grid_n = NULL,
                          point_alpha = 0.55,
                          point_size = 1.6,
                          reverse_x = FALSE,
                          reverse_y = FALSE,
                          facet_scales = c("fixed", "free", "free_x", "free_y")) {
  .require_ggplot2()
  density <- match.arg(density)
  bw <- match.arg(bw)
  .check_bw_scale(bw_scale)
  facet_scales <- match.arg(facet_scales)
  .check_bool(points, "points")
  .check_bool(overlap, "overlap")
  .check_bool(annotate, "annotate")
  .check_bool(reverse_x, "reverse_x")
  .check_bool(reverse_y, "reverse_y")
  .check_plot_number(point_alpha, "point_alpha", lower = 0, upper = 1)
  .check_plot_number(point_size, "point_size", lower = 0)
  .check_positive_count(min_tokens, "min_tokens")
  .check_conf_level(conf_level)
  .check_eval_seed(eval_seed)
  .check_contrast_levels(levels)
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot < 0 || n_boot != as.integer(n_boot)) {
    stop("`n_boot` must be a single non-negative integer.", call. = FALSE)
  }
  n_boot <- as.integer(n_boot)

  if (!is.character(features) || !length(features) || length(features) > 2L) {
    stop(
      "`features` must be one or two column names. For higher-dimensional ",
      "spaces use plot_category_pca() for display and compute metrics on the ",
      "full feature set.",
      call. = FALSE
    )
  }
  d <- length(features)
  grid_n <- .check_contrast_grid_n(grid_n, d)

  keep_cols <- c(group_col, category_col, features)
  .check_columns(data, keep_cols)
  data <- .metric_data(data, keep_cols)
  .check_numeric_features(data, features)
  levs <- .two_levels(data[[category_col]], "category_col")

  # ---- split into panels and compute layer data per panel ----
  grouped <- !is.null(group_col)
  if (grouped) {
    group_col <- .check_group_cols(group_col)
    groups <- .split_groups(data, group_col)
  } else {
    groups <- list(data)
  }

  min_cat <- .kde_min_category_tokens(d)
  layer <- list(curves = list(), regions = list(), ov = list(), tokens = list())
  kept_labels <- character(0)
  skipped <- character(0)
  kept_data <- list()

  for (df_g in groups) {
    label <- if (grouped) .group_label(df_g, group_col) else "global"
    counts <- table(factor(df_g[[category_col]], levels = levs))
    if (nrow(df_g) < min_tokens || any(counts < min_cat)) {
      skipped <- c(skipped, label)
      next
    }
    comp <- .contrast_panel_data(
      df_g = df_g, features = features, category_col = category_col,
      levs = levs, density = density, bw = bw, levels = levels,
      grid_n = grid_n, label = label, bw_scale = bw_scale
    )
    layer$curves[[label]] <- comp$curves
    layer$regions[[label]] <- comp$regions
    layer$ov[[label]] <- comp$ov
    layer$tokens[[label]] <- comp$tokens
    kept_labels <- c(kept_labels, label)
    kept_data[[label]] <- df_g
  }

  if (!length(kept_labels)) {
    stop(
      "plot_contrast(): no group met `min_tokens` = ", min_tokens,
      " with at least ", min_cat, " tokens per category.",
      call. = FALSE
    )
  }
  if (length(skipped)) {
    shown <- utils::head(skipped, 5L)
    warning(
      "plot_contrast(): skipped ", length(skipped),
      " group(s) below `min_tokens` or without both categories: ",
      paste(shown, collapse = ", "),
      if (length(skipped) > 5L) ", ..." else "",
      call. = FALSE
    )
  }

  curves <- do.call(rbind, layer$curves)
  regions <- do.call(rbind, layer$regions)
  ov <- do.call(rbind, layer$ov)
  tokens <- do.call(rbind, layer$tokens)
  rownames(curves) <- rownames(ov) <- rownames(tokens) <- NULL
  if (!is.null(regions)) rownames(regions) <- NULL

  # In 2D the overlap drives an alpha channel: normalize across panels so
  # shading strength is comparable (a fainter panel genuinely overlaps less),
  # and drop near-zero cells. In 1D the overlap is a ribbon on the density
  # axis and must keep its raw scale.
  if (d == 2L && !is.null(ov) && nrow(ov)) {
    ov_max <- max(ov$ov, na.rm = TRUE)
    if (is.finite(ov_max) && ov_max > 0) {
      ov$ov <- ov$ov / ov_max
    }
    ov <- ov[ov$ov > 0.004, , drop = FALSE]
  }

  plot_df <- do.call(rbind, kept_data)

  # ---- assemble ----
  p <- if (d == 1L) {
    .build_contrast_1d(
      curves = curves, ov = ov, tokens = tokens, feature = features[[1]],
      category_col = category_col, points = points, overlap = overlap,
      point_alpha = point_alpha
    )
  } else {
    .build_contrast_2d(
      regions = regions, ov = ov, tokens = tokens, features = features,
      category_col = category_col, points = points, overlap = overlap,
      point_alpha = point_alpha, point_size = point_size
    )
  }

  # ---- accountability: annotate with the metrics this model produces ----
  ann <- NULL
  if (isTRUE(annotate)) {
    ann <- tryCatch(
      .contrast_annotation_table(
        data = plot_df, features = features, category_col = category_col,
        group_col = group_col, density = density, bw = bw,
        bw_scale = bw_scale,
        min_tokens = min_tokens, mc_n = mc_n, eval_seed = eval_seed,
        n_boot = n_boot, conf_level = conf_level
      ),
      error = function(e) {
        warning(
          "plot_contrast(): metric annotation failed and was skipped: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(ann)) {
      lab_df <- .contrast_annotation_labels(ann, grouped = grouped)
      if (!is.null(lab_df) && nrow(lab_df)) {
        # Anchor at the visual top-left; +/-Inf swap sides under reversed
        # scales, so the anchor follows the reverse_* flags.
        lab_df$.x <- if (isTRUE(reverse_x)) Inf else -Inf
        lab_df$.y <- if (isTRUE(reverse_y) && d == 2L) -Inf else Inf
        p <- p + ggplot2::geom_label(
          data = lab_df,
          ggplot2::aes(
            x = .data[[".x"]], y = .data[[".y"]], label = .data[[".lab"]]
          ),
          inherit.aes = FALSE,
          hjust = 0, vjust = 1,
          size = 2.9, lineheight = 1.1,
          label.size = 0, fill = "white", alpha = 0.75,
          color = .phontrast_ink
        )
      }
    }
  }

  if (isTRUE(reverse_x)) p <- p + ggplot2::scale_x_reverse()
  if (isTRUE(reverse_y) && d == 2L) p <- p + ggplot2::scale_y_reverse()
  if (grouped) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(.data[[".group"]]), scales = facet_scales
    )
  }

  p <- p + ggplot2::labs(
    caption = .contrast_caption(
      density = density, bw = bw, bw_scale = bw_scale,
      levels = if (d == 2L) levels else NULL,
      mc_n = mc_n, eval_seed = eval_seed, n_boot = n_boot,
      conf_level = conf_level, data = plot_df,
      category_col = category_col, levs = levs,
      grouped = grouped, n_groups = length(kept_labels)
    )
  )

  p <- .phontrast_style(p, fill = (d == 1L))
  if (!is.null(ann)) {
    attr(p, "contrast_metrics") <- ann
  }
  p
}

# ---- validation helpers ----------------------------------------------------

.check_contrast_levels <- function(levels) {
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) ||
      any(levels <= 0) || any(levels >= 1) || anyDuplicated(levels)) {
    stop(
      "`levels` must be distinct probabilities strictly between 0 and 1.",
      call. = FALSE
    )
  }
  invisible(sort(levels))
}

.check_contrast_grid_n <- function(grid_n, d) {
  if (is.null(grid_n)) {
    return(if (d == 1L) 512L else 151L)
  }
  .check_positive_count(grid_n, "grid_n")
  as.integer(grid_n)
}

# ---- per-panel density computation -----------------------------------------

# Matches the mvnorm backend's default ridge (`.mvnorm_mc_pair()`).
.contrast_ridge <- 1e-6

.contrast_panel_data <- function(df_g, features, category_col, levs,
                                 density, bw, levels, grid_n, label,
                                 bw_scale = 1) {
  if (length(features) == 1L) {
    .contrast_panel_1d(
      df_g = df_g, feature = features[[1]], category_col = category_col,
      levs = levs, density = density, bw = bw, grid_n = grid_n, label = label,
      bw_scale = bw_scale
    )
  } else {
    .contrast_panel_2d(
      df_g = df_g, features = features, category_col = category_col,
      levs = levs, density = density, bw = bw, levels = levels,
      grid_n = grid_n, label = label, bw_scale = bw_scale
    )
  }
}

.contrast_panel_1d <- function(df_g, feature, category_col, levs,
                               density, bw, grid_n, label, bw_scale = 1) {
  x1 <- df_g[df_g[[category_col]] == levs[1], feature]
  x2 <- df_g[df_g[[category_col]] == levs[2], feature]

  if (identical(density, "kde")) {
    # Same univariate bandwidth selection (and scaling) as the 1-D metric path.
    h1 <- .scale_bandwidth(.select_univariate_bandwidth(x1, bw), bw_scale)
    h2 <- .scale_bandwidth(.select_univariate_bandwidth(x2, bw), bw_scale)
    pad <- 3 * max(h1, h2)
    grid <- seq(min(x1, x2) - pad, max(x1, x2) + pad, length.out = grid_n)
    d1 <- .kde_1d_values(x1, grid, h1)
    d2 <- .kde_1d_values(x2, grid, h2)
  } else {
    fit1 <- .mvn_fit_category(matrix(x1, ncol = 1), .contrast_ridge)
    fit2 <- .mvn_fit_category(matrix(x2, ncol = 1), .contrast_ridge)
    s1 <- sqrt(fit1$S[1, 1]); s2 <- sqrt(fit2$S[1, 1])
    pad <- 3.5 * max(s1, s2)
    grid <- seq(min(x1, x2) - pad, max(x1, x2) + pad, length.out = grid_n)
    d1 <- stats::dnorm(grid, fit1$mu, s1)
    d2 <- stats::dnorm(grid, fit2$mu, s2)
  }

  curves <- data.frame(
    .group = label,
    category = factor(rep(levs, each = length(grid)), levels = levs),
    x = c(grid, grid),
    dens = c(d1, d2),
    stringsAsFactors = FALSE
  )
  ov <- data.frame(.group = label, x = grid, ov = pmin(d1, d2))
  tokens <- data.frame(
    .group = label,
    category = factor(df_g[[category_col]], levels = levs),
    x = df_g[[feature]],
    stringsAsFactors = FALSE
  )
  list(curves = curves, regions = NULL, ov = ov, tokens = tokens)
}

.contrast_panel_2d <- function(df_g, features, category_col, levs,
                               density, bw, levels, grid_n, label,
                               bw_scale = 1) {
  X1 <- as.matrix(df_g[df_g[[category_col]] == levs[1], features, drop = FALSE])
  X2 <- as.matrix(df_g[df_g[[category_col]] == levs[2], features, drop = FALSE])

  if (identical(density, "kde")) {
    # Same bandwidth selectors (and scaling) and ks evaluation as the KDE
    # metric path.
    H1 <- .scale_bandwidth(.select_multivariate_bandwidth(X1, bw, levs[1]), bw_scale)
    H2 <- .scale_bandwidth(.select_multivariate_bandwidth(X2, bw, levs[2]), bw_scale)
    pad <- 3 * sqrt(pmax(diag(H1), diag(H2)))
    grid <- .contrast_grid_2d(rbind(X1, X2), pad, grid_n)
    z1 <- .kde_grid_values(X1, H1, grid)
    z2 <- .kde_grid_values(X2, H2, grid)
    cell <- diff(grid$gx[1:2]) * diff(grid$gy[1:2])
    regions <- rbind(
      .hdr_contour_paths(grid, z1, .hdr_levels(z1, cell, levels),
                         levs[1], label),
      .hdr_contour_paths(grid, z2, .hdr_levels(z2, cell, levels),
                         levs[2], label)
    )
  } else {
    fit1 <- .mvn_fit_category(X1, .contrast_ridge)
    fit2 <- .mvn_fit_category(X2, .contrast_ridge)
    if (!isTRUE(tryCatch({ chol(fit1$S); chol(fit2$S); TRUE },
                         error = function(e) FALSE))) {
      stop(
        "plot_contrast(): a category covariance is not positive definite; ",
        "features may be collinear.",
        call. = FALSE
      )
    }
    pad <- 2.9 * sqrt(pmax(diag(fit1$S), diag(fit2$S)))
    grid <- .contrast_grid_2d(rbind(X1, X2), pad, grid_n)
    z1 <- exp(.mvn_logdens(grid$mat, fit1$mu, fit1$S))
    z2 <- exp(.mvn_logdens(grid$mat, fit2$mu, fit2$S))
    regions <- rbind(
      .gaussian_ellipse_paths(fit1$mu, fit1$S, levels, levs[1], label),
      .gaussian_ellipse_paths(fit2$mu, fit2$S, levels, levs[2], label)
    )
  }

  if (!is.null(regions)) {
    regions$category <- factor(regions$category, levels = levs)
  }
  ov <- data.frame(
    .group = label,
    x = grid$mat[, 1], y = grid$mat[, 2],
    ov = pmin(z1, z2),
    .w = diff(grid$gx[1:2]),
    .h = diff(grid$gy[1:2])
  )
  tokens <- data.frame(
    .group = label,
    category = factor(df_g[[category_col]], levels = levs),
    x = df_g[[features[[1]]]],
    y = df_g[[features[[2]]]],
    stringsAsFactors = FALSE
  )
  list(curves = NULL, regions = regions, ov = ov, tokens = tokens)
}

.contrast_grid_2d <- function(X, pad, grid_n) {
  gx <- seq(min(X[, 1]) - pad[1], max(X[, 1]) + pad[1], length.out = grid_n)
  gy <- seq(min(X[, 2]) - pad[2], max(X[, 2]) + pad[2], length.out = grid_n)
  list(gx = gx, gy = gy, mat = as.matrix(expand.grid(gx, gy)))
}

.kde_grid_values <- function(X, H, grid) {
  # Exact KDE evaluation at the grid points via ks::kde(), the same engine the
  # KDE metrics use (`engine = "ks"`).
  fhat <- ks::kde(x = X, H = H, eval.points = grid$mat)
  pmax(as.numeric(fhat$estimate), 0)
}

.hdr_levels <- function(z, cell_area, probs) {
  # Density thresholds whose superlevel sets enclose `probs` of the (grid-
  # discretized) probability mass: highest-density regions.
  o <- order(z, decreasing = TRUE)
  cum_mass <- cumsum(z[o]) * cell_area
  total <- cum_mass[length(cum_mass)]
  vapply(probs, function(p) {
    z[o][which(cum_mass >= p * total)[1]]
  }, numeric(1))
}

.hdr_contour_paths <- function(grid, z, breaks, category, label) {
  zmat <- matrix(z, nrow = length(grid$gx))
  lines <- grDevices::contourLines(grid$gx, grid$gy, zmat, levels = breaks)
  if (!length(lines)) {
    return(NULL)
  }
  out <- lapply(seq_along(lines), function(i) {
    seg <- lines[[i]]
    data.frame(
      .group = label,
      category = category,
      .piece = paste(label, category, i, sep = "\r"),
      x = seg$x, y = seg$y,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

.gaussian_ellipse_paths <- function(mu, S, levels, category, label, n = 181L) {
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(theta), sin(theta))
  R <- chol(S)
  out <- lapply(seq_along(levels), function(i) {
    r <- sqrt(stats::qchisq(levels[i], df = 2))
    pts <- sweep(r * (circle %*% R), 2, mu, "+")
    data.frame(
      .group = label,
      category = category,
      .piece = paste(label, category, i, sep = "\r"),
      x = pts[, 1], y = pts[, 2],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# ---- plot assembly ---------------------------------------------------------

.build_contrast_1d <- function(curves, ov, tokens, feature, category_col,
                               points, overlap, point_alpha) {
  p <- ggplot2::ggplot()

  if (isTRUE(overlap) && !is.null(ov) && nrow(ov)) {
    p <- p + ggplot2::geom_area(
      data = ov,
      ggplot2::aes(x = .data[["x"]], y = .data[["ov"]]),
      fill = "grey30", alpha = 0.35, color = NA
    )
  }

  p <- p +
    ggplot2::geom_ribbon(
      data = curves,
      ggplot2::aes(
        x = .data[["x"]], ymin = 0, ymax = .data[["dens"]],
        fill = .data[["category"]]
      ),
      alpha = 0.16, color = NA, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = curves,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["dens"]], color = .data[["category"]]
      ),
      linewidth = 0.7
    )

  if (isTRUE(points) && !is.null(tokens) && nrow(tokens)) {
    p <- p + ggplot2::geom_rug(
      data = tokens,
      ggplot2::aes(x = .data[["x"]], color = .data[["category"]]),
      alpha = point_alpha, sides = "b", show.legend = FALSE
    )
  }

  p + ggplot2::labs(
    x = feature, y = "density",
    color = category_col, fill = category_col
  )
}

.build_contrast_2d <- function(regions, ov, tokens, features, category_col,
                               points, overlap, point_alpha, point_size) {
  p <- ggplot2::ggplot()

  if (isTRUE(overlap) && !is.null(ov) && nrow(ov)) {
    # geom_tile with explicit cell sizes: robust under reversed axes and
    # per-panel grids, where geom_raster's even-spacing check misfires.
    p <- p +
      ggplot2::geom_tile(
        data = ov,
        ggplot2::aes(
          x = .data[["x"]], y = .data[["y"]], alpha = .data[["ov"]],
          width = .data[[".w"]], height = .data[[".h"]]
        ),
        fill = "grey25"
      ) +
      ggplot2::scale_alpha_continuous(range = c(0, 0.5), guide = "none")
  }

  if (isTRUE(points) && !is.null(tokens) && nrow(tokens)) {
    p <- p + ggplot2::geom_point(
      data = tokens,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], color = .data[["category"]]
      ),
      alpha = point_alpha, size = point_size
    )
  }

  if (!is.null(regions) && nrow(regions)) {
    p <- p + ggplot2::geom_path(
      data = regions,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        color = .data[["category"]], group = .data[[".piece"]]
      ),
      linewidth = 0.55, show.legend = FALSE
    )
  }

  p + ggplot2::labs(
    x = features[[1]], y = features[[2]], color = category_col
  )
}

# ---- accountability layer --------------------------------------------------

.contrast_annotation_table <- function(data, features, category_col, group_col,
                                       density, bw, min_tokens, mc_n,
                                       eval_seed, n_boot, conf_level,
                                       bw_scale = 1) {
  # The annotation values are not recomputed ad hoc: they are the phontrast()
  # estimates under the same density model that the plot draws.
  phontrast(
    data = data,
    features = features,
    category_col = category_col,
    group_col = group_col,
    metrics = c("jsd", "overlap"),
    min_tokens = min_tokens,
    bw = bw,
    bw_scale = bw_scale,
    density = density,
    mc_n = mc_n,
    eval_seed = eval_seed,
    do_boot = n_boot > 0L,
    n_boot = if (n_boot > 0L) n_boot else 1000,
    conf_level = conf_level,
    progress = FALSE,
    output = "wide"
  )
}

.contrast_annotation_labels <- function(ann, grouped) {
  if (is.null(ann) || !nrow(ann)) {
    return(NULL)
  }
  fmt <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "NA")
  has_ci <- all(c("jsd_ci_lower", "jsd_ci_upper") %in% names(ann))
  jsd_part <- paste0(
    "JSD ", fmt(ann$jsd),
    if (has_ci) {
      ifelse(
        is.finite(ann$jsd_ci_lower) & is.finite(ann$jsd_ci_upper),
        paste0(" [", fmt(ann$jsd_ci_lower), ", ", fmt(ann$jsd_ci_upper), "]"),
        ""
      )
    } else ""
  )
  has_ov_ci <- all(
    c("percent_overlap_ci_lower", "percent_overlap_ci_upper") %in% names(ann)
  )
  ov_part <- paste0(
    "overlap ", fmt(ann$percent_overlap),
    if (has_ov_ci) {
      ifelse(
        is.finite(ann$percent_overlap_ci_lower) &
          is.finite(ann$percent_overlap_ci_upper),
        paste0(
          " [", fmt(ann$percent_overlap_ci_lower), ", ",
          fmt(ann$percent_overlap_ci_upper), "]"
        ),
        ""
      )
    } else ""
  )
  data.frame(
    .group = if (grouped) as.character(ann$group) else "global",
    .lab = paste0(jsd_part, "\n", ov_part),
    stringsAsFactors = FALSE
  )
}

.contrast_caption <- function(density, bw, levels, mc_n, eval_seed, n_boot,
                              conf_level, data, category_col, levs, grouped,
                              n_groups, bw_scale = 1) {
  # `levels = NULL` means no regions are drawn (the 1-D plot), so the caption
  # omits the region clause rather than describing layers that are not there.
  region_lab <- if (is.null(levels)) {
    ""
  } else {
    paste0(
      "; regions: ", paste0(round(100 * levels), "%", collapse = "/"),
      if (identical(density, "kde")) " highest-density" else " coverage"
    )
  }
  model <- if (identical(density, "kde")) {
    paste0(
      "density: kde (bw = ", bw,
      if (bw_scale != 1) paste0(" x", format(bw_scale)) else "",
      ")", region_lab
    )
  } else {
    paste0(
      "density: mvnorm (mc_n = ", format(mc_n, big.mark = ","),
      if (!is.null(eval_seed)) paste0(", seed = ", eval_seed) else "",
      ")", region_lab
    )
  }
  n_lab <- if (grouped) {
    paste0("N = ", nrow(data), " tokens in ", n_groups, " group(s)")
  } else {
    counts <- table(factor(data[[category_col]], levels = levs))
    paste0("n: ", paste(levs, as.integer(counts), collapse = ", "))
  }
  boot_lab <- if (n_boot > 0L) {
    paste0("; CI: ", round(100 * conf_level), "% bootstrap (",
           n_boot, " resamples)")
  } else {
    ""
  }
  paste0(model, "; ", n_lab, boot_lab)
}
