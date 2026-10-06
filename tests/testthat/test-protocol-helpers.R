# Tests for the building blocks of the Sec. VII.A protocol added in phontrast
# 2.5.0: the `bw_scale` bandwidth multiplier on the kernel path,
# percentile_rank(), and recommended_estimator().

protocol_fixture <- function(n = 60, seed = 2026) {
  set.seed(seed)
  data.frame(
    speaker = rep(c("s01", "s02"), each = 2 * n),
    vowel = rep(rep(c("ih", "eh"), each = n), 2),
    f1 = c(rnorm(n, 500, 55), rnorm(n, 600, 60),
           rnorm(n, 510, 60), rnorm(n, 640, 65)),
    f2 = c(rnorm(n, 1980, 150), rnorm(n, 1800, 155),
           rnorm(n, 1960, 160), rnorm(n, 1760, 165)),
    stringsAsFactors = FALSE
  )
}

# ---- percentile_rank() -------------------------------------------------------

test_that("percentile_rank() is (average rank - 1/2) / k with ties averaged", {
  x <- c(0.2, 0.5, 0.5, 0.9)
  # ranks 1, 2.5, 2.5, 4 over k = 4
  expect_equal(percentile_rank(x), c(0.125, 0.5, 0.5, 0.875))
  # a single measurement sits at the middle of its own ordering
  expect_equal(percentile_rank(3), 0.5)
  # larger values get larger percentile ranks
  expect_equal(order(percentile_rank(c(3, 1, 2))), c(2, 3, 1))
})

test_that("percentile_rank() drops excluded and missing values from k", {
  x <- c(0.2, 0.5, 0.5, 0.995)
  pr <- percentile_rank(x, exclude = x >= 0.99)
  expect_equal(pr, c((1 - 0.5) / 3, (2.5 - 0.5) / 3, (2.5 - 0.5) / 3, NA))
  expect_equal(percentile_rank(c(NA, 1, 2)), c(NA, 0.25, 0.75))
  # NA in `exclude` counts as not excluded
  expect_equal(percentile_rank(c(1, 2), exclude = c(NA, FALSE)), c(0.25, 0.75))
  expect_true(all(is.na(percentile_rank(c(NA_real_, NA_real_)))))
  expect_true(all(is.na(percentile_rank(c(1, 2), exclude = c(TRUE, TRUE)))))
})

test_that("percentile_rank() validates its inputs", {
  expect_error(percentile_rank("a"), "numeric")
  expect_error(percentile_rank(1:3, exclude = TRUE), "same length")
  expect_error(percentile_rank(1:3, exclude = c(1, 0, 1)), "logical")
})

# ---- recommended_estimator() -------------------------------------------------

test_that("recommended_estimator() follows the three Table III tiers", {
  low <- recommended_estimator(2)
  expect_identical(low$bw, "Hpi")
  expect_identical(low$engine, "ks")
  expect_null(low$eval_n)
  expect_true(low$loo)
  expect_identical(recommended_estimator(4)$tier, low$tier)

  mid <- recommended_estimator(8)
  expect_identical(mid$bw, "scott.diag")
  expect_identical(mid$engine, "fast_diag")
  expect_identical(mid$eval_n, 200L)
  expect_true(mid$loo)
  # untabulated dimensionalities take the next higher calibrated row
  expect_identical(recommended_estimator(5)$tier, mid$tier)
  expect_identical(recommended_estimator(13)$tier, mid$tier)

  high <- recommended_estimator(32)
  expect_identical(high$bw, "scott.diag")
  expect_identical(high$engine, "fast_diag")
  expect_false(high$loo)
  expect_identical(recommended_estimator(14)$tier, high$tier)
  expect_identical(recommended_estimator(64)$d, 64L)
})

test_that("recommended_estimator() settings are accepted by jsd_kde_nd()", {
  d <- protocol_fixture(n = 30)
  est <- recommended_estimator(2)
  expect_no_error(
    jsd_kde_nd(d, c("f1", "f2"), "vowel", bw = est$bw, engine = est$engine,
               eval_n = est$eval_n, loo = est$loo)
  )
  est <- recommended_estimator(8)
  expect_no_error(
    jsd_kde_nd(d, c("f1", "f2"), "vowel", bw = est$bw, engine = est$engine,
               eval_n = est$eval_n, loo = est$loo)
  )
})

test_that("recommended_estimator() validates `d`", {
  expect_error(recommended_estimator(0), "positive integer")
  expect_error(recommended_estimator(2.5), "positive integer")
  expect_error(recommended_estimator(c(2, 3)), "positive integer")
  expect_error(recommended_estimator("2"), "positive integer")
})

# ---- bw_scale ----------------------------------------------------------------

test_that(".scale_bandwidth() scales h linearly and H quadratically", {
  expect_equal(phontrast:::.scale_bandwidth(2, 0.5), 1)
  H <- diag(c(1, 4))
  expect_equal(phontrast:::.scale_bandwidth(H, 2), 4 * H)
  expect_identical(phontrast:::.scale_bandwidth(H, 1), H)
})

test_that("bw_scale is validated and bw_scale = 1 reproduces the default", {
  d <- protocol_fixture(n = 30)
  base <- jsd_kde_nd(d, c("f1", "f2"), "vowel")
  expect_equal(jsd_kde_nd(d, c("f1", "f2"), "vowel", bw_scale = 1), base)
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "2")) {
    expect_error(jsd_kde_nd(d, c("f1", "f2"), "vowel", bw_scale = bad), "bw_scale")
  }
  expect_error(percent_overlap_kde(d, "f1", "vowel", bw_scale = 0), "bw_scale")
  expect_error(estimate_jsd(d, "f1", "vowel", bw_scale = -2), "bw_scale")
  expect_error(estimate_overlap(d, "f1", "vowel", bw_scale = NA), "bw_scale")
  expect_error(
    phontrast(d, "f1", "vowel", metrics = "jsd", bw_scale = 0),
    "bw_scale"
  )
})

test_that("bw_scale moves the kernel JSD and doubling the bandwidth lowers it", {
  d <- protocol_fixture()
  feats <- c("f1", "f2")
  half <- jsd_kde_nd(d, feats, "vowel", bw = "scott.diag", bw_scale = 0.5)
  one  <- jsd_kde_nd(d, feats, "vowel", bw = "scott.diag")
  dbl  <- jsd_kde_nd(d, feats, "vowel", bw = "scott.diag", bw_scale = 2)
  # Oversmoothing flattens both densities toward each other, so the doubled
  # bandwidth must read lower; the halved bandwidth moves the estimate but its
  # direction at two dimensions depends on the sample (the leave-one-out
  # correction offsets the undersmoothing), so only a change is asserted.
  expect_gt(one, dbl)
  expect_false(isTRUE(all.equal(half, one)))

  # the one-dimensional path scales h directly
  half1 <- jsd_kde_nd(d, "f1", "vowel", bw = "scott.diag", bw_scale = 0.5)
  dbl1  <- jsd_kde_nd(d, "f1", "vowel", bw = "scott.diag", bw_scale = 2)
  expect_gt(half1, dbl1)

  # the fast_diag engine and the legacy estimator honour it too
  expect_gt(
    jsd_kde_nd(d, feats, "vowel", bw = "scott.diag", engine = "fast_diag", bw_scale = 0.5),
    jsd_kde_nd(d, feats, "vowel", bw = "scott.diag", engine = "fast_diag", bw_scale = 2)
  )
  expect_false(isTRUE(all.equal(
    jsd_kde_nd(d, feats, "vowel", method = "legacy", bw_scale = 2),
    jsd_kde_nd(d, feats, "vowel", method = "legacy")
  )))
})

test_that("bw_scale threads through estimate_jsd(), estimate_overlap(), phontrast()", {
  d <- protocol_fixture()
  feats <- c("f1", "f2")
  j2 <- estimate_jsd(d, feats, "vowel", bw = "scott.diag", bw_scale = 2)$jsd_point
  expect_equal(j2, jsd_kde_nd(d, feats, "vowel", bw = "scott.diag", bw_scale = 2))

  o2 <- estimate_overlap(d, feats, "vowel", bw = "scott.diag", bw_scale = 2)$overlap
  expect_equal(o2, percent_overlap_kde(d, feats, "vowel", bw = "scott.diag", bw_scale = 2))
  expect_false(isTRUE(all.equal(
    o2, percent_overlap_kde(d, feats, "vowel", bw = "scott.diag")
  )))

  p <- phontrast(d, feats, "vowel", metrics = c("jsd", "overlap"),
                 bw = "scott.diag", bw_scale = 2)
  expect_equal(p$jsd, j2)
  expect_equal(p$percent_overlap, o2)

  # grouped path (speaker_jsd() forwards it to jsd_kde_nd())
  g <- estimate_jsd(d, feats, "vowel", group_col = "speaker",
                    bw = "scott.diag", bw_scale = 2)
  expect_equal(
    g$jsd_point[g$group == "s01"],
    jsd_kde_nd(d[d$speaker == "s01", ], feats, "vowel", bw = "scott.diag", bw_scale = 2)
  )
})

test_that("bw_scale is ignored under density = 'mvnorm'", {
  d <- protocol_fixture(n = 30)
  a <- jsd_kde_nd(d, c("f1", "f2"), "vowel", density = "mvnorm",
                  mc_n = 2000, eval_seed = 1)
  b <- jsd_kde_nd(d, c("f1", "f2"), "vowel", density = "mvnorm",
                  mc_n = 2000, eval_seed = 1, bw_scale = 3)
  expect_equal(a, b)
})
