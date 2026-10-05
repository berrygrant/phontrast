#' Rank-agreement plot for a speaker ranking
#'
#' Draws step 3 of the ranking protocol behind \code{rank_contrasts()}: each
#' speaker's Pillai percentile rank against its Jensen-Shannon distance
#' percentile rank, the identity line on which the two measures agree, and the
#' inspection band of \code{margin} either side of it. Speakers outside the
#' band are flagged; speakers the bandwidth check set aside are crossed;
#' speakers whose flag may not be read at this sample size and dimensionality
#' are hollow. Speakers at the \eqn{\sqrt{JSD}} ceiling have no rank and are
#' listed in the caption instead of drawn.
#'
#' @param ranking A \code{phontrast_ranking} object from
#'   \code{rank_contrasts()}.
#' @param label Which speakers to name on the plot: \code{"flagged"}
#'   (default: flagged or set-aside speakers), \code{"all"}, or
#'   \code{"none"}.
#' @param base_size Base font size passed to \code{theme_phontrast()}.
#'
#' @return A \pkg{ggplot2} plot object.
#'
#' @seealso \code{rank_contrasts()}, \code{inspect_contrast()}.
#'
#' @examples
#' set.seed(2026)
#' gaps <- c(30, 60, 90, 120, 150, 180, 210, 240)
#' cohort <- do.call(rbind, lapply(seq_along(gaps), function(i) {
#'   data.frame(
#'     speaker = sprintf("s%02d", i),
#'     vowel = rep(c("ih", "eh"), each = 50),
#'     f1 = c(rnorm(50, 500, 55), rnorm(50, 500 + gaps[i], 55)),
#'     f2 = c(rnorm(50, 1900, 120), rnorm(50, 1900 - gaps[i], 120))
#'   )
#' }))
#' ranking <- rank_contrasts(cohort, c("f1", "f2"), "vowel", "speaker",
#'                           bw_check = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_rank_agreement(ranking, label = "all")
#' }
#' @export
plot_rank_agreement <- function(ranking,
                                label = c("flagged", "all", "none"),
                                base_size = 12) {
  .require_ggplot2()
  label <- match.arg(label)
  .check_ranking(ranking)
  pa <- attr(ranking, "protocol")
  margin <- pa$margin

  df <- as.data.frame(ranking)
  df$status <- .ranking_status(df)
  drawn <- df[!is.na(df$pr_jsd) & !is.na(df$pr_pillai), , drop = FALSE]

  band <- data.frame(x = c(0, 1), lo = c(0, 1) - margin, hi = c(0, 1) + margin)
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = .data[["x"]], ymin = .data[["lo"]], ymax = .data[["hi"]]),
      fill = "grey88", alpha = 0.7
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = .phontrast_ink_muted, linewidth = 0.45
    )

  if (nrow(drawn)) {
    p <- p + ggplot2::geom_point(
      data = drawn,
      ggplot2::aes(
        x = .data[["pr_pillai"]], y = .data[["pr_jsd"]],
        colour = .data[["status"]], shape = .data[["status"]]
      ),
      size = 2.7, stroke = 0.9
    )
    to_label <- switch(
      label,
      none = drawn[0, , drop = FALSE],
      all = drawn,
      flagged = drawn[drawn$status %in% c("flagged", "set aside"), , drop = FALSE]
    )
    if (nrow(to_label)) {
      p <- p + ggplot2::geom_text(
        data = to_label,
        ggplot2::aes(
          x = .data[["pr_pillai"]], y = .data[["pr_jsd"]],
          label = .data[["group"]]
        ),
        nudge_y = 0.035, vjust = 0, size = 3, colour = .phontrast_ink,
        show.legend = FALSE
      )
    }
  }

  p <- p +
    ggplot2::scale_colour_manual(values = .ranking_status_colours(), drop = FALSE) +
    ggplot2::scale_shape_manual(values = .ranking_status_shapes(), drop = FALSE) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Pillai percentile rank",
      y = expression(sqrt(JSD) ~ "percentile rank"),
      colour = "status", shape = "status",
      caption = .ranking_caption(df, pa)
    ) +
    theme_phontrast(base_size = base_size)
  p
}

#' Plot a rank_contrasts() result directly
#'
#' \code{rank_contrasts()} results carry the class \code{"phontrast_ranking"},
#' so \code{plot(ranking)} draws the rank-agreement plot of
#' \code{plot_rank_agreement()} and \code{ggplot2::autoplot(ranking)} returns
#' it unprinted.
#'
#' @param x,object A \code{phontrast_ranking} object.
#' @param ... Passed on to \code{plot_rank_agreement()}.
#'
#' @return \code{plot()} draws the plot and returns it invisibly;
#'   \code{autoplot()} returns the \pkg{ggplot2} object unprinted.
#' @export
plot.phontrast_ranking <- function(x, ...) {
  p <- plot_rank_agreement(x, ...)
  print(p)
  invisible(p)
}

#' @rdname plot.phontrast_ranking
#' @exportS3Method ggplot2::autoplot
autoplot.phontrast_ranking <- function(object, ...) {
  plot_rank_agreement(object, ...)
}

#' Inspect one speaker's contrast across bandwidths
#'
#' Redraws a single speaker from a \code{rank_contrasts()} ranking with
#' \code{plot_contrast()}'s distribution-aware layers, one panel per kernel
#' bandwidth, so a flagged or set-aside speaker can be inspected the way the
#' protocol asks: at half, at the selected, and at twice the diagonal Scott
#' bandwidth (the smoothing used by the bandwidth check). Each panel is
#' labelled with the Jensen-Shannon distance and shared mass at that
#' bandwidth and with the speaker's Pillai trace, which takes no bandwidth;
#' the subtitle restates the reported ranks, the rank difference, and the
#' outcome of the bandwidth check. The tokens come from the ranking itself
#' (\code{attr(ranking, "protocol")$data}), so no further data is needed.
#'
#' @param ranking A \code{phontrast_ranking} object from
#'   \code{rank_contrasts()}.
#' @param group String; the speaker (a value of \code{ranking$group}) to
#'   inspect.
#' @param features Features to display: by default the ranking's own features,
#'   which must number one or two. For a ranking in more than two dimensions,
#'   name two of its features here; the annotated metrics are still computed
#'   on all of the ranking's features.
#' @param bw_scales Positive bandwidth multipliers, one panel each (default
#'   \code{c(0.5, 1, 2)}).
#' @param levels,points,overlap,grid_n,point_alpha,point_size,reverse_x,reverse_y
#'   As in \code{plot_contrast()}.
#'
#' @return A \pkg{ggplot2} plot object. The per-panel metrics are attached as
#'   \code{attr(p, "inspect_metrics")}.
#'
#' @seealso \code{rank_contrasts()}, \code{plot_rank_agreement()},
#'   \code{plot_contrast()}.
#'
#' @examples
#' set.seed(2026)
#' gaps <- c(30, 60, 90, 120, 150, 180, 210, 240)
#' cohort <- do.call(rbind, lapply(seq_along(gaps), function(i) {
#'   data.frame(
#'     speaker = sprintf("s%02d", i),
#'     vowel = rep(c("ih", "eh"), each = 50),
#'     f1 = c(rnorm(50, 500, 55), rnorm(50, 500 + gaps[i], 55)),
#'     f2 = c(rnorm(50, 1900, 120), rnorm(50, 1900 - gaps[i], 120))
#'   )
#' }))
#' ranking <- rank_contrasts(cohort, c("f1", "f2"), "vowel", "speaker",
#'                           bw_check = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   inspect_contrast(ranking, "s04")
#' }
#' @export
inspect_contrast <- function(ranking,
                             group,
                             features = NULL,
                             bw_scales = c(0.5, 1, 2),
                             levels = c(0.5, 0.8, 0.95),
                             points = TRUE,
                             overlap = TRUE,
                             grid_n = NULL,
                             point_alpha = 0.55,
                             point_size = 1.6,
                             reverse_x = FALSE,
                             reverse_y = FALSE) {
  .require_ggplot2()
  .check_ranking(ranking)
  pa <- attr(ranking, "protocol")
  if (!is.character(group) || length(group) != 1L || !group %in% ranking$group) {
    stop(
      "`group` must be one of the ranked groups: ",
      paste(utils::head(ranking$group, 10L), collapse = ", "),
      if (nrow(ranking) > 10L) ", ..." else "", ".",
      call. = FALSE
    )
  }
  if (is.null(features)) {
    features <- pa$features
  }
  if (!is.character(features) || !length(features) || length(features) > 2L ||
      !all(features %in% pa$features)) {
    stop(
      "`features` must name one or two of the ranking's features (",
      paste(pa$features, collapse = ", "), ") to display.",
      call. = FALSE
    )
  }
  if (!is.numeric(bw_scales) || !length(bw_scales) || any(!is.finite(bw_scales)) ||
      any(bw_scales <= 0) || anyDuplicated(bw_scales)) {
    stop("`bw_scales` must be distinct positive finite numbers.", call. = FALSE)
  }
  .check_bool(points, "points")
  .check_bool(overlap, "overlap")
  .check_bool(reverse_x, "reverse_x")
  .check_bool(reverse_y, "reverse_y")
  .check_plot_number(point_alpha, "point_alpha", lower = 0, upper = 1)
  .check_plot_number(point_size, "point_size", lower = 0)
  .check_contrast_levels(levels)
  d <- length(features)
  grid_n <- .check_contrast_grid_n(grid_n, d)

  data <- pa$data
  category_col <- pa$category_col
  levs <- .two_levels(data[[category_col]], "category_col")
  df_g <- .split_groups(data, pa$group_col)[[group]]
  row <- as.data.frame(ranking)[ranking$group == group, , drop = FALSE]
  est <- pa$estimator

  # ---- one panel per bandwidth, drawn and scored at that bandwidth ----
  layer <- list(curves = list(), regions = list(), ov = list(), tokens = list())
  metrics <- vector("list", length(bw_scales))
  for (i in seq_along(bw_scales)) {
    s <- bw_scales[i]
    panel <- sprintf("bandwidth x%s", format(s))
    comp <- .contrast_panel_data(
      df_g = df_g, features = features, category_col = category_col,
      levs = levs, density = "kde", bw = "scott.diag", levels = levels,
      grid_n = grid_n, label = panel, bw_scale = s
    )
    layer$curves[[panel]] <- comp$curves
    layer$regions[[panel]] <- comp$regions
    layer$ov[[panel]] <- comp$ov
    layer$tokens[[panel]] <- comp$tokens
    metrics[[i]] <- .inspect_panel_metrics(
      df_g = df_g, features = pa$features, category_col = category_col,
      est = est, bw_scale = s, eval_seed = pa$eval_seed, panel = panel
    )
  }
  metrics <- do.call(rbind, metrics)
  metrics$pillai <- row$pillai
  rownames(metrics) <- NULL

  curves <- do.call(rbind, layer$curves)
  regions <- do.call(rbind, layer$regions)
  ov <- do.call(rbind, layer$ov)
  tokens <- do.call(rbind, layer$tokens)
  rownames(ov) <- rownames(tokens) <- NULL
  if (!is.null(curves)) rownames(curves) <- NULL
  if (!is.null(regions)) rownames(regions) <- NULL
  if (d == 2L && !is.null(ov) && nrow(ov)) {
    ov_max <- max(ov$ov, na.rm = TRUE)
    if (is.finite(ov_max) && ov_max > 0) ov$ov <- ov$ov / ov_max
    ov <- ov[ov$ov > 0.004, , drop = FALSE]
  }

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

  fmt <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "NA")
  lab_df <- data.frame(
    .group = metrics$panel,
    .lab = paste0(
      "sqrt(JSD) ", fmt(metrics$sqrt_jsd), "\n",
      "shared mass ", fmt(metrics$shared_mass), "\n",
      "Pillai ", fmt(metrics$pillai)
    ),
    stringsAsFactors = FALSE
  )
  lab_df$.x <- if (isTRUE(reverse_x)) Inf else -Inf
  lab_df$.y <- if (isTRUE(reverse_y) && d == 2L) -Inf else Inf
  p <- p + ggplot2::geom_label(
    data = lab_df,
    ggplot2::aes(x = .data[[".x"]], y = .data[[".y"]], label = .data[[".lab"]]),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.9, lineheight = 1.1,
    label.size = 0, fill = "white", alpha = 0.75, color = .phontrast_ink
  )

  if (isTRUE(reverse_x)) p <- p + ggplot2::scale_x_reverse()
  if (isTRUE(reverse_y) && d == 2L) p <- p + ggplot2::scale_y_reverse()
  p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[".group"]]))

  counts <- table(factor(df_g[[category_col]], levels = levs))
  p <- p + ggplot2::labs(
    title = group,
    subtitle = .inspect_subtitle(row, pa),
    caption = paste0(
      "density: kde (bw = scott.diag x",
      paste(vapply(bw_scales, format, character(1)), collapse = "/"), ")",
      if (d == 2L) paste0("; regions: ", paste0(round(100 * levels), "%", collapse = "/"),
                          " highest-density") else "",
      "; n: ", paste(levs, as.integer(counts), collapse = ", "),
      "; metrics on ", pa$d, " feature", if (pa$d == 1L) "" else "s",
      " (", paste(pa$features, collapse = ", "), ")"
    )
  )
  p <- .phontrast_style(p, fill = (d == 1L))
  attr(p, "inspect_metrics") <- metrics
  p
}

# ---- helpers -----------------------------------------------------------------

.check_ranking <- function(ranking) {
  if (!inherits(ranking, "phontrast_ranking") || is.null(attr(ranking, "protocol"))) {
    stop("`ranking` must be the result of rank_contrasts().", call. = FALSE)
  }
  invisible(ranking)
}

.ranking_status_levels <- c("agrees", "flagged", "set aside", "flag not readable")

.ranking_status <- function(df) {
  status <- rep("agrees", nrow(df))
  status[is.na(df$flag)] <- "flag not readable"
  status[df$flag %in% TRUE] <- "flagged"
  if ("set_aside" %in% names(df)) {
    status[df$set_aside %in% TRUE] <- "set aside"
  }
  factor(status, levels = .ranking_status_levels)
}

.ranking_status_colours <- function() {
  stats::setNames(
    unname(.phontrast_colors[c("blue", "vermillion", "orange", "grey")]),
    .ranking_status_levels
  )
}

.ranking_status_shapes <- function() {
  stats::setNames(c(16, 16, 4, 1), .ranking_status_levels)
}

.ranking_caption <- function(df, pa) {
  fmt_floor <- function(f) if (is.finite(f)) format(f) else "none"
  at_ceiling <- df$group[df$at_ceiling %in% TRUE]
  parts <- c(
    sprintf("k = %d ranked", pa$k),
    if (length(at_ceiling)) {
      sprintf("at ceiling (sqrt(JSD) >= %s, not drawn): %s", format(pa$ceiling),
              paste(at_ceiling, collapse = ", "))
    },
    sprintf("margin %s", format(pa$margin)),
    sprintf("d = %d: rank floor %s, flag floor %s tokens per category", pa$d,
            fmt_floor(pa$floors$rank_floor), fmt_floor(pa$floors$flag_floor)),
    sprintf("kernel %s/%s", pa$estimator$bw, pa$estimator$engine),
    if (isTRUE(pa$bw_check)) {
      sprintf("bandwidth check x0.5/x2: %d set aside", sum(df$set_aside %in% TRUE))
    } else {
      "bandwidth check not run"
    }
  )
  paste(parts, collapse = "; ")
}

.inspect_panel_metrics <- function(df_g, features, category_col, est, bw_scale,
                                   eval_seed, panel) {
  mc <- tryCatch(
    .kde_mc_pair(
      data = df_g, features = features, category_col = category_col,
      bw = "scott.diag", eval_n = est$eval_n, eval_seed = eval_seed,
      engine = est$engine, chunk_size = 1000L, metric = "inspect_contrast()",
      bw_scale = bw_scale
    ),
    error = function(e) NULL
  )
  data.frame(
    panel = panel,
    bw_scale = bw_scale,
    sqrt_jsd = if (is.null(mc)) NA_real_ else sqrt(.jsd_mc(mc, loo = est$loo)),
    shared_mass = if (is.null(mc)) NA_real_ else .overlap_mc(mc),
    stringsAsFactors = FALSE
  )
}

.inspect_subtitle <- function(row, pa) {
  fmt <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "NA")
  if (isTRUE(row$at_ceiling)) {
    first <- sprintf(
      "reported sqrt(JSD) %s: at ceiling, excluded from rank comparisons; Pillai %s; shared mass %s",
      fmt(row$sqrt_jsd), fmt(row$pillai), fmt(row$shared_mass)
    )
  } else {
    verdict <- if (is.na(row$flag)) {
      "flag not readable"
    } else if (isTRUE(row$flag)) {
      "flagged"
    } else {
      "agrees"
    }
    first <- sprintf(
      "reported sqrt(JSD) %s (rank %s), Pillai %s (rank %s): rank_diff %+.2f, %s",
      fmt(row$sqrt_jsd), fmt(row$pr_jsd), fmt(row$pillai), fmt(row$pr_pillai),
      row$rank_diff, verdict
    )
  }
  if (!isTRUE(pa$bw_check) || !"bw_shift" %in% names(row)) {
    return(first)
  }
  second <- if (is.na(row$bw_shift)) {
    "bandwidth check: no rank to shift"
  } else {
    sprintf(
      "bandwidth check: rank shift %+.2f%s, %s",
      row$bw_shift,
      if (isTRUE(row$sign_change)) ", rank difference changes sign" else "",
      if (isTRUE(row$set_aside)) "set aside" else "not set aside"
    )
  }
  paste0(first, "\n", second)
}
