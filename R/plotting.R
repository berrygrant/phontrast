#' Plot overlap metric comparisons
#'
#' Visualizes the output of \code{phontrast()} with \pkg{ggplot2}.
#' The input may be either wide or long output. By default, values are plotted on
#' a separation-oriented scale so overlap-oriented metrics are transformed as
#' \code{1 - estimate}.
#'
#' @param metrics Data frame returned by \code{phontrast()}.
#' @param value Scale to plot: \code{"separation"} plots
#'   \code{separation_value}; \code{"estimate"} plots raw metric estimates.
#' @param metric Optional character vector of metric display names to include.
#' @param group_col Optional column to use on the x-axis. Defaults to
#'   \code{"group"} when present, then \code{"scope"}.
#' @param show_ci Logical; if \code{TRUE}, draw confidence intervals when
#'   \code{ci_lower}/\code{ci_upper} columns are present.
#' @param facet Logical; if \code{TRUE}, facet by metric.
#' @param sort Logical; if \code{TRUE}, order comparisons by their mean plotted
#'   value.
#'
#' @return A \pkg{ggplot2} plot object.
#' @examples
#' set.seed(2026)
#' vowels <- data.frame(
#'   speaker = rep(c("s01", "s02"), each = 60),
#'   vowel = rep(rep(c("ih", "eh"), each = 30), 2),
#'   f1 = c(rnorm(30, 500, 55), rnorm(30, 560, 60),
#'          rnorm(30, 510, 60), rnorm(30, 575, 65)),
#'   f2 = c(rnorm(30, 1980, 150), rnorm(30, 1880, 155),
#'          rnorm(30, 1960, 160), rnorm(30, 1840, 165))
#' )
#'
#' metrics <- phontrast(
#'   vowels,
#'   features = c("f1", "f2"),
#'   category_col = "vowel",
#'   group_col = "speaker",
#'   output = "long"
#' )
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_overlap_metrics(metrics)
#' }
#' @export
plot_overlap_metrics <- function(metrics,
                                 value = c("separation", "estimate"),
                                 metric = NULL,
                                 group_col = NULL,
                                 show_ci = TRUE,
                                 facet = TRUE,
                                 sort = TRUE) {
  .require_ggplot2()
  value <- match.arg(value)
  .check_bool(show_ci, "show_ci")
  .check_bool(facet, "facet")
  .check_bool(sort, "sort")

  plot_data <- .as_overlap_metrics_long(metrics)
  if (!nrow(plot_data)) {
    stop("`metrics` has no rows to plot.", call. = FALSE)
  }

  if (!is.null(metric)) {
    plot_data <- plot_data[plot_data$metric %in% metric, , drop = FALSE]
    if (!nrow(plot_data)) {
      stop("No rows matched `metric`.", call. = FALSE)
    }
  }

  if (identical(value, "separation")) {
    if (!"separation_value" %in% names(plot_data)) {
      if (!all(c("estimate", "orientation") %in% names(plot_data))) {
        stop(
          "`metrics` must contain `separation_value` or both `estimate` and `orientation`.",
          call. = FALSE
        )
      }
      plot_data$separation_value <- ifelse(
        plot_data$orientation == "overlap",
        1 - plot_data$estimate,
        plot_data$estimate
      )
    }
    plot_data$.plot_value <- plot_data$separation_value
    y_lab <- "Separation-oriented value"
  } else {
    if (!"estimate" %in% names(plot_data)) {
      stop("`metrics` must contain an `estimate` column.", call. = FALSE)
    }
    plot_data$.plot_value <- plot_data$estimate
    y_lab <- "Metric estimate"
  }

  x_col <- .choose_metric_x_col(plot_data, group_col)
  if (identical(x_col, ".comparison") && !".comparison" %in% names(plot_data)) {
    plot_data$.comparison <- seq_len(nrow(plot_data))
  }
  x_lab <- if (identical(x_col, ".comparison")) "comparison" else x_col
  plot_data$.plot_x <- as.character(plot_data[[x_col]])
  if (isTRUE(sort)) {
    plot_data$.plot_x <- .order_discrete_axis(plot_data$.plot_x, plot_data$.plot_value)
  }

  plot_data$metric <- factor(plot_data$metric, levels = unique(plot_data$metric))

  has_ci <- isTRUE(show_ci) && all(c("ci_lower", "ci_upper") %in% names(plot_data))
  if (has_ci) {
    ci <- .plot_metric_ci(plot_data, value)
    plot_data$.ci_lower <- ci$lower
    plot_data$.ci_upper <- ci$upper
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[".plot_x"]],
      y = .data[[".plot_value"]],
      color = .data[["orientation"]]
    )
  )

  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[".ci_lower"]], ymax = .data[[".ci_upper"]]),
      width = 0.18,
      alpha = 0.75,
      na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    ggplot2::labs(x = x_lab, y = y_lab, color = "Direction") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(facet)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[["metric"]]), scales = "free_y")
  } else {
    p <- p + ggplot2::aes(shape = .data[["metric"]])
  }

  p
}

#' Plot phonological categories in acoustic space
#'
#' Creates a \pkg{ggplot2} visualization of one or two acoustic dimensions from
#' a token-level table. One-dimensional inputs are plotted as density curves;
#' two-dimensional inputs are plotted as category-colored scatterplots with
#' optional normal ellipses.
#'
#' @param data Data frame containing category labels and acoustic features.
#' @param features One or two numeric feature columns to plot.
#' @param category_col String; category column.
#' @param group_col Optional grouping column used for facets.
#' @param points Logical; if \code{TRUE}, show observed tokens.
#' @param ellipses Logical; if \code{TRUE}, add normal ellipses for two-feature
#'   plots when enough observations are available.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#' @param reverse_x Logical; if \code{TRUE}, reverse the x-axis.
#' @param reverse_y Logical; if \code{TRUE}, reverse the y-axis.
#' @param equal_axes Logical; if \code{TRUE}, use a fixed coordinate ratio for
#'   two-feature plots.
#' @param facet_scales Scales passed to \code{ggplot2::facet_wrap()} when
#'   \code{group_col} is supplied.
#'
#' @return A \pkg{ggplot2} plot object.
#' @examples
#' set.seed(2026)
#' vowels <- data.frame(
#'   vowel = rep(c("ih", "eh"), each = 40),
#'   f1 = c(rnorm(40, 500, 55), rnorm(40, 565, 60)),
#'   f2 = c(rnorm(40, 1980, 150), rnorm(40, 1870, 155))
#' )
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_category_space(vowels, features = "f1", category_col = "vowel")
#'   plot_category_space(
#'     vowels,
#'     features = c("f2", "f1"),
#'     category_col = "vowel",
#'     reverse_x = TRUE,
#'     reverse_y = TRUE
#'   )
#' }
#' @export
plot_category_space <- function(data,
                                features,
                                category_col,
                                group_col = NULL,
                                points = TRUE,
                                ellipses = TRUE,
                                point_alpha = 0.65,
                                point_size = 1.8,
                                reverse_x = FALSE,
                                reverse_y = FALSE,
                                equal_axes = FALSE,
                                facet_scales = c("fixed", "free", "free_x", "free_y")) {
  .require_ggplot2()
  facet_scales <- match.arg(facet_scales)
  .check_bool(points, "points")
  .check_bool(ellipses, "ellipses")
  .check_bool(reverse_x, "reverse_x")
  .check_bool(reverse_y, "reverse_y")
  .check_bool(equal_axes, "equal_axes")
  .check_plot_number(point_alpha, "point_alpha", lower = 0, upper = 1)
  .check_plot_number(point_size, "point_size", lower = 0)

  if (!is.character(features) || !length(features) || length(features) > 2L) {
    stop("`features` must be one or two column names.", call. = FALSE)
  }
  keep_cols <- c(category_col, features, group_col)
  .check_columns(data, keep_cols)
  plot_data <- .metric_data(data, keep_cols)
  .check_numeric_features(plot_data, features)
  if (.observed_n_categories(plot_data[[category_col]]) < 2L) {
    stop("`category_col` must contain at least two observed categories.", call. = FALSE)
  }

  plot_data[[category_col]] <- .drop_unused_levels(plot_data[[category_col]])
  if (!is.null(group_col)) {
    plot_data[[group_col]] <- .drop_unused_levels(plot_data[[group_col]])
  }

  if (length(features) == 1L) {
    p <- .plot_category_space_1d(
      plot_data = plot_data,
      feature = features[[1]],
      category_col = category_col,
      group_col = group_col,
      points = points,
      point_alpha = point_alpha,
      reverse_x = reverse_x,
      facet_scales = facet_scales
    )
  } else {
    p <- .plot_category_space_2d(
      plot_data = plot_data,
      features = features,
      category_col = category_col,
      group_col = group_col,
      points = points,
      ellipses = ellipses,
      point_alpha = point_alpha,
      point_size = point_size,
      reverse_x = reverse_x,
      reverse_y = reverse_y,
      equal_axes = equal_axes,
      facet_scales = facet_scales
    )
  }

  p
}

#' Plot a PCA projection of multidimensional category space
#'
#' Projects an arbitrary multidimensional acoustic feature set onto two
#' principal components and visualizes the result with \pkg{ggplot2}. This is a
#' diagnostic plot for high-dimensional workflows: metric estimates should still
#' be computed on the full feature set when that is the intended analysis.
#'
#' @param data Data frame containing category labels and acoustic features.
#' @param features Character vector of two or more numeric feature columns used
#'   for PCA.
#' @param category_col String; category column.
#' @param group_col Optional grouping column used for facets.
#' @param components Two positive integers giving principal components to plot.
#' @param center,scale. Passed to \code{stats::prcomp()}.
#' @param points Logical; if \code{TRUE}, show projected observations.
#' @param ellipses Logical; if \code{TRUE}, add normal ellipses when enough
#'   observations are available.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#' @param equal_axes Logical; if \code{TRUE}, use a fixed coordinate ratio.
#' @param facet_scales Scales passed to \code{ggplot2::facet_wrap()} when
#'   \code{group_col} is supplied.
#'
#' @return A \pkg{ggplot2} plot object. The fitted \code{prcomp} object and
#'   variance-explained table are stored as \code{"pca"} and
#'   \code{"variance_explained"} attributes.
#' @examples
#' set.seed(2026)
#' features <- paste0("mfcc", 1:5)
#' vowels <- data.frame(
#'   vowel = rep(c("ih", "eh"), each = 50),
#'   matrix(rnorm(100 * length(features)), ncol = length(features))
#' )
#' names(vowels)[-1] <- features
#' vowels[vowels$vowel == "eh", features[1:2]] <- vowels[vowels$vowel == "eh", features[1:2]] + 0.8
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_category_pca(vowels, features = features, category_col = "vowel")
#' }
#' @export
plot_category_pca <- function(data,
                              features,
                              category_col,
                              group_col = NULL,
                              components = c(1L, 2L),
                              center = TRUE,
                              scale. = TRUE,
                              points = TRUE,
                              ellipses = TRUE,
                              point_alpha = 0.65,
                              point_size = 1.8,
                              equal_axes = TRUE,
                              facet_scales = c("fixed", "free", "free_x", "free_y")) {
  .require_ggplot2()
  facet_scales <- match.arg(facet_scales)
  .check_bool(center, "center")
  .check_bool(scale., "scale.")
  .check_bool(points, "points")
  .check_bool(ellipses, "ellipses")
  .check_bool(equal_axes, "equal_axes")
  .check_plot_number(point_alpha, "point_alpha", lower = 0, upper = 1)
  .check_plot_number(point_size, "point_size", lower = 0)

  if (!is.character(features) || length(features) < 2L) {
    stop("`features` must contain at least two numeric feature columns.", call. = FALSE)
  }
  if (!is.numeric(components) || length(components) != 2L ||
      any(!is.finite(components)) || any(components < 1) ||
      any(components != as.integer(components))) {
    stop("`components` must contain two positive integers.", call. = FALSE)
  }
  components <- as.integer(components)

  keep_cols <- c(category_col, features, group_col)
  .check_columns(data, keep_cols)
  plot_data <- .metric_data(data, keep_cols)
  .check_numeric_features(plot_data, features)
  if (.observed_n_categories(plot_data[[category_col]]) < 2L) {
    stop("`category_col` must contain at least two observed categories.", call. = FALSE)
  }

  max_component <- min(nrow(plot_data), length(features))
  if (any(components > max_component)) {
    stop(
      "`components` must be no larger than the number of available principal components (",
      max_component, ").",
      call. = FALSE
    )
  }

  fit <- stats::prcomp(plot_data[, features, drop = FALSE], center = center, scale. = scale.)
  scores <- as.data.frame(fit$x[, components, drop = FALSE])
  score_cols <- paste0("PC", components)
  names(scores) <- score_cols
  scores[[category_col]] <- .drop_unused_levels(plot_data[[category_col]])
  if (!is.null(group_col)) {
    scores[[group_col]] <- .drop_unused_levels(plot_data[[group_col]])
  }

  variance_explained <- (fit$sdev^2) / sum(fit$sdev^2)
  axis_labels <- paste0(
    score_cols,
    " (",
    format(round(100 * variance_explained[components], 1), nsmall = 1),
    "%)"
  )

  p <- .plot_category_space_2d(
    plot_data = scores,
    features = score_cols,
    category_col = category_col,
    group_col = group_col,
    points = points,
    ellipses = ellipses,
    point_alpha = point_alpha,
    point_size = point_size,
    reverse_x = FALSE,
    reverse_y = FALSE,
    equal_axes = equal_axes,
    facet_scales = facet_scales
  ) +
    ggplot2::labs(
      x = axis_labels[[1]],
      y = axis_labels[[2]]
    )

  attr(p, "pca") <- fit
  attr(p, "variance_explained") <- data.frame(
    component = paste0("PC", seq_along(variance_explained)),
    variance_explained = variance_explained,
    stringsAsFactors = FALSE
  )
  p
}

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "`ggplot2` is required for plotting. Install it with install.packages(\"ggplot2\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.check_bool <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

.check_plot_number <- function(x, arg, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop("`", arg, "` must be a single finite number.", call. = FALSE)
  }
}

.as_overlap_metrics_long <- function(metrics) {
  if (!is.data.frame(metrics)) {
    stop("`metrics` must be a data frame returned by phontrast().", call. = FALSE)
  }
  if (all(c("metric", "estimate") %in% names(metrics))) {
    return(as.data.frame(metrics))
  }
  if (any(.compare_metric_columns() %in% names(metrics))) {
    return(as.data.frame(.comparison_long(metrics)))
  }
  stop(
    "`metrics` must be wide or long output from phontrast().",
    call. = FALSE
  )
}

.choose_metric_x_col <- function(plot_data, group_col = NULL) {
  if (!is.null(group_col)) {
    if (!group_col %in% names(plot_data)) {
      stop("`group_col` must name a column in `metrics`.", call. = FALSE)
    }
    return(group_col)
  }
  if ("group" %in% names(plot_data)) {
    return("group")
  }
  if ("scope" %in% names(plot_data)) {
    return("scope")
  }
  plot_data$.comparison <- seq_len(nrow(plot_data))
  ".comparison"
}

.order_discrete_axis <- function(x, y) {
  means <- stats::aggregate(y, list(x = x), mean, na.rm = TRUE)
  means <- means[order(means$x, decreasing = FALSE), , drop = FALSE]
  means <- means[order(means[, 2], decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  factor(x, levels = means$x)
}

.plot_metric_ci <- function(plot_data, value) {
  lower <- plot_data$ci_lower
  upper <- plot_data$ci_upper

  if (identical(value, "separation") && "orientation" %in% names(plot_data)) {
    overlap <- plot_data$orientation == "overlap"
    lower_out <- lower
    upper_out <- upper
    lower_out[overlap] <- 1 - upper[overlap]
    upper_out[overlap] <- 1 - lower[overlap]
    lower <- lower_out
    upper <- upper_out
  }

  data.frame(lower = lower, upper = upper)
}

.drop_unused_levels <- function(x) {
  if (is.factor(x)) {
    return(droplevels(x))
  }
  x
}

.plot_category_space_1d <- function(plot_data,
                                    feature,
                                    category_col,
                                    group_col = NULL,
                                    points = TRUE,
                                    point_alpha = 0.65,
                                    reverse_x = FALSE,
                                    facet_scales = "fixed") {
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[feature]],
      color = .data[[category_col]],
      fill = .data[[category_col]]
    )
  ) +
    ggplot2::geom_density(alpha = 0.25, na.rm = TRUE) +
    ggplot2::labs(x = feature, y = "Density", color = category_col, fill = category_col) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (isTRUE(points)) {
    p <- p + ggplot2::geom_rug(alpha = point_alpha, sides = "b", na.rm = TRUE)
  }
  if (isTRUE(reverse_x)) {
    p <- p + ggplot2::scale_x_reverse()
  }
  if (!is.null(group_col)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[group_col]]), scales = facet_scales)
  }

  p
}

.plot_category_space_2d <- function(plot_data,
                                    features,
                                    category_col,
                                    group_col = NULL,
                                    points = TRUE,
                                    ellipses = TRUE,
                                    point_alpha = 0.65,
                                    point_size = 1.8,
                                    reverse_x = FALSE,
                                    reverse_y = FALSE,
                                    equal_axes = FALSE,
                                    facet_scales = "fixed") {
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[features[[1]]]],
      y = .data[[features[[2]]]],
      color = .data[[category_col]]
    )
  )

  if (isTRUE(points)) {
    p <- p + ggplot2::geom_point(alpha = point_alpha, size = point_size, na.rm = TRUE)
  }
  if (isTRUE(ellipses)) {
    ellipse_data <- .filter_ellipse_data(plot_data, category_col, group_col)
    if (nrow(ellipse_data)) {
      p <- p + ggplot2::stat_ellipse(
        data = ellipse_data,
        type = "norm",
        linewidth = 0.7,
        show.legend = FALSE,
        na.rm = TRUE
      )
    }
  }

  p <- p +
    ggplot2::labs(x = features[[1]], y = features[[2]], color = category_col) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (isTRUE(reverse_x)) {
    p <- p + ggplot2::scale_x_reverse()
  }
  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (isTRUE(equal_axes)) {
    p <- p + ggplot2::coord_equal()
  }
  if (!is.null(group_col)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[group_col]]), scales = facet_scales)
  }

  p
}

.filter_ellipse_data <- function(plot_data, category_col, group_col = NULL) {
  if (is.null(group_col)) {
    counts <- table(plot_data[[category_col]])
    keep <- names(counts)[counts >= 3L]
    return(plot_data[plot_data[[category_col]] %in% keep, , drop = FALSE])
  }

  key <- paste(plot_data[[group_col]], plot_data[[category_col]], sep = "\r")
  counts <- table(key)
  keep <- names(counts)[counts >= 3L]
  plot_data[key %in% keep, , drop = FALSE]
}
