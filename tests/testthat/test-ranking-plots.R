# Tests for the plotting side of the ranking protocol added in phontrast
# 2.5.0: plot_rank_agreement() (and the plot()/autoplot() methods on
# rank_contrasts() results), inspect_contrast(), and bw_scale in
# plot_contrast().

.ranking_plot_cohort <- function(n = 100, seed = 5) {
  set.seed(seed)
  gaps <- seq(25, 200, length.out = 8)
  make <- function(id, gap) data.frame(
    speaker = id,
    vowel = rep(c("ih", "eh"), each = n),
    f1 = c(rnorm(n, 500, 55), rnorm(n, 500 + gap, 55)),
    f2 = c(rnorm(n, 1900, 120), rnorm(n, 1900 - gap, 120)),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, Map(make, sprintf("s%02d", seq_along(gaps)), gaps))
  side <- sample(c(-1, 1), n, replace = TRUE)
  rbind(out, data.frame(
    speaker = "plant",
    vowel = rep(c("ih", "eh"), each = n),
    f1 = c(rnorm(n, 500, 30), 500 + 95 * side + rnorm(n, 0, 30)),
    f2 = c(rnorm(n, 1900, 70), 1900 - 220 * side + rnorm(n, 0, 70)),
    stringsAsFactors = FALSE
  ))
}

.ranking_layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1], character(1))
}

cohort_plot <- .ranking_plot_cohort()
ranking_plot <- rank_contrasts(cohort_plot, c("f1", "f2"), "vowel", "speaker")

test_that("plot_rank_agreement() draws band, identity, points, and labels", {
  skip_if_not_installed("ggplot2")
  p <- plot_rank_agreement(ranking_plot)
  expect_s3_class(p, "ggplot")
  geoms <- .ranking_layer_geoms(p)
  expect_true(all(c("GeomRibbon", "GeomAbline", "GeomPoint", "GeomText") %in% geoms))
  built <- ggplot2::ggplot_build(p)
  pts <- built$data[[which(geoms == "GeomPoint")]]
  expect_equal(nrow(pts), 9L)
  labs <- built$data[[which(geoms == "GeomText")]]
  expect_true("plant" %in% labs$label)
  expect_false("s01" %in% labs$label)
  expect_match(p$labels$caption, "k = 9 ranked")
  expect_match(p$labels$caption, "margin 0.25")
  expect_match(p$labels$x, "Pillai")

  p_all <- plot_rank_agreement(ranking_plot, label = "all")
  built_all <- ggplot2::ggplot_build(p_all)
  g_all <- .ranking_layer_geoms(p_all)
  expect_equal(nrow(built_all$data[[which(g_all == "GeomText")]]), 9L)

  p_none <- plot_rank_agreement(ranking_plot, label = "none")
  expect_false("GeomText" %in% .ranking_layer_geoms(p_none))
})

test_that("status encodes flagged, set-aside, and unreadable speakers", {
  df <- as.data.frame(ranking_plot)
  status <- phontrast:::.ranking_status(df)
  expect_s3_class(status, "factor")
  # the planted speaker is flagged, and the bandwidth check sets it aside:
  # the doubled bandwidth smooths away the bimodality its flag rests on, and
  # "set aside" outranks "flagged"
  expect_true(df$flag[df$group == "plant"])
  expect_true(df$set_aside[df$group == "plant"])
  expect_identical(as.character(status[df$group == "plant"]), "set aside")
  df_kept <- df
  df_kept$set_aside[df_kept$group == "plant"] <- FALSE
  expect_identical(
    as.character(phontrast:::.ranking_status(df_kept)[df_kept$group == "plant"]),
    "flagged"
  )
  expect_true(all(as.character(status[grepl("^s", df$group)]) == "agrees"))
  df$flag[1] <- NA
  expect_identical(as.character(phontrast:::.ranking_status(df)[1]), "flag not readable")
  df$set_aside[2] <- TRUE
  expect_identical(as.character(phontrast:::.ranking_status(df)[2]), "set aside")
  # without the bandwidth-check columns the status still resolves
  df$set_aside <- NULL
  expect_equal(nlevels(phontrast:::.ranking_status(df)), 4L)
})

test_that("plot() and autoplot() dispatch on phontrast_ranking", {
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(ranking_plot)
  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::autoplot(ranking_plot, label = "none"), "ggplot")
  expect_error(plot_rank_agreement(data.frame(a = 1)), "rank_contrasts")
})

test_that("inspect_contrast() draws one panel per bandwidth with that bandwidth's metrics", {
  skip_if_not_installed("ggplot2")
  p <- inspect_contrast(ranking_plot, "plant")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 3L)
  m <- attr(p, "inspect_metrics")
  expect_equal(m$bw_scale, c(0.5, 1, 2))
  expect_true(all(is.finite(m$sqrt_jsd)))
  expect_true(all(is.finite(m$shared_mass)))
  row <- ranking_plot[ranking_plot$group == "plant", ]
  expect_equal(m$sqrt_jsd[m$bw_scale == 0.5], row$sqrt_jsd_half)
  expect_equal(m$sqrt_jsd[m$bw_scale == 2], row$sqrt_jsd_double)
  expect_true(all(m$pillai == row$pillai))
  expect_identical(p$labels$title, "plant")
  expect_match(p$labels$subtitle, "flagged")
  expect_match(p$labels$subtitle, "bandwidth check")
  expect_match(p$labels$caption, "pooled scott.diag x0.5/1/2")
  geoms <- .ranking_layer_geoms(p)
  expect_true(all(c("GeomTile", "GeomPoint", "GeomPath", "GeomLabel") %in% geoms))
})

test_that("inspect_contrast() honours bw_scales, display options, and validates", {
  skip_if_not_installed("ggplot2")
  p <- inspect_contrast(ranking_plot, "s03", bw_scales = c(1, 3), points = FALSE,
                        overlap = FALSE, reverse_x = TRUE, reverse_y = TRUE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 2L)
  geoms <- .ranking_layer_geoms(p)
  expect_false(any(c("GeomTile", "GeomPoint") %in% geoms))
  expect_match(p$labels$subtitle, "agrees")
  expect_error(inspect_contrast(ranking_plot, "nobody"), "ranked groups")
  expect_error(inspect_contrast(ranking_plot, c("s01", "s02")), "ranked groups")
  expect_error(inspect_contrast(ranking_plot, "s01", bw_scales = c(1, 1)), "bw_scales")
  expect_error(inspect_contrast(ranking_plot, "s01", bw_scales = -1), "bw_scales")
  expect_error(inspect_contrast(ranking_plot, "s01", features = "f9"), "features")
  expect_error(inspect_contrast(ranking_plot, "s01", features = c("f1", "f2", "f1")), "features")
})

test_that("inspect_contrast() works for one-feature rankings and without the bandwidth check", {
  skip_if_not_installed("ggplot2")
  small <- cohort_plot[cohort_plot$speaker != "plant", ]
  r1 <- rank_contrasts(small, "f1", "vowel", "speaker", bw_check = FALSE)
  p <- inspect_contrast(r1, "s05")
  expect_s3_class(p, "ggplot")
  expect_true("GeomLine" %in% .ranking_layer_geoms(p))
  expect_false(grepl("bandwidth check", p$labels$subtitle))
  expect_equal(nrow(attr(p, "inspect_metrics")), 3L)
  # the agreement plot also works without set_aside columns
  expect_s3_class(plot_rank_agreement(r1), "ggplot")
})

test_that("annotation labels use whichever border argument ggplot2 supports", {
  skip_if_not_installed("ggplot2")
  # Newer ggplot2 deprecates geom_label(label.size =) for the `linewidth`
  # aesthetic; older releases only know `label.size`.
  expect_identical(phontrast:::.label_borderless(c("colour", "fill")), list(label.size = 0))
  expect_identical(phontrast:::.label_borderless(c("colour", "linewidth")), list(linewidth = 0))
  # Whatever is installed, building the annotated plots must not warn.
  d <- cohort_plot[cohort_plot$speaker == "s04", ]
  expect_no_warning(p <- plot_contrast(d, c("f1", "f2"), "vowel", bw = "scott.diag"))
  expect_no_warning(inspect_contrast(ranking_plot, "s03"))
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLabel" %in% geoms)
})

test_that("plot_contrast() accepts bw_scale and records it", {
  skip_if_not_installed("ggplot2")
  d <- cohort_plot[cohort_plot$speaker == "s04", ]
  p <- plot_contrast(d, c("f1", "f2"), "vowel", bw = "scott.diag", bw_scale = 2)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "scott.diag x2")
  ann <- attr(p, "contrast_metrics")
  expect_equal(ann$jsd, jsd_kde_nd(d, c("f1", "f2"), "vowel", bw = "scott.diag", bw_scale = 2))
  p1 <- plot_contrast(d, "f1", "vowel", bw_scale = 0.5)
  expect_match(p1$labels$caption, "x0.5")
  expect_error(plot_contrast(d, "f1", "vowel", bw_scale = 0), "bw_scale")
  expect_false(grepl(" x1", plot_contrast(d, "f1", "vowel")$labels$caption))
})
