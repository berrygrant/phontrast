# Tests for the shared-density kernel family behind phontrast() (2.5.0): the
# refactor that reads Jensen-Shannon, overlap, total variation, and the
# matched-kernel Bhattacharyya / Hellinger distances off one density estimate,
# plus the Euclidean distance of standardized means.

.family_fixture <- function(n = 80, seed = 2026) {
  set.seed(seed)
  data.frame(
    speaker = rep(c("s01", "s02"), each = 2 * n),
    vowel = rep(rep(c("ih", "eh"), each = n), 2),
    f1 = c(rnorm(n, 500, 55), rnorm(n, 590, 60),
           rnorm(n, 510, 60), rnorm(n, 640, 65)),
    f2 = c(rnorm(n, 1980, 150), rnorm(n, 1820, 155),
           rnorm(n, 1960, 160), rnorm(n, 1760, 165)),
    stringsAsFactors = FALSE
  )
}

test_that("the shared pass reproduces the single-metric estimators", {
  d <- .family_fixture()
  feats <- c("f1", "f2")
  p <- phontrast(d, feats, "vowel", metrics = c("jsd", "js_distance", "overlap"))
  expect_equal(p$jsd, jsd_kde_nd(d, feats, "vowel"))
  expect_equal(p$js_distance, sqrt(p$jsd))
  expect_equal(p$percent_overlap, percent_overlap_kde(d, feats, "vowel"))

  g <- phontrast(d, feats, "vowel", group_col = "speaker",
                 metrics = c("jsd", "overlap", "tv"))
  s02 <- d[d$speaker == "s02", ]
  expect_equal(g$jsd[g$group == "s02"], jsd_kde_nd(s02, feats, "vowel"))
  expect_equal(g$percent_overlap[g$group == "s02"], percent_overlap_kde(s02, feats, "vowel"))
  expect_equal(g$total_variation, 1 - g$percent_overlap)

  # bw_scale reaches the shared pass
  p2 <- phontrast(d, feats, "vowel", metrics = "jsd", bw = "scott.diag", bw_scale = 2)
  expect_equal(p2$jsd, jsd_kde_nd(d, feats, "vowel", bw = "scott.diag", bw_scale = 2))
})

test_that("default metrics and column order are unchanged", {
  d <- .family_fixture(n = 40)
  p <- phontrast(d, c("f1", "f2"), "vowel")
  expect_named(p, c(
    "scope", "n_tokens", "pillai", "pillai_p_value", "bhatt_dist",
    "bhatt_affinity", "jsd", "js_distance", "mahalanobis_dist", "percent_overlap"
  ))
  expect_false(any(c("total_variation", "hellinger", "euclidean_dist") %in% names(p)))
})

test_that(".bhatt_mc() returns 1 for identical densities and matches the vector helper", {
  mc <- list(
    logp1 = log(c(0.2, 0.5, 0.3)), logq1 = log(c(0.2, 0.5, 0.3)),
    logp2 = log(c(0.1, 0.4)), logq2 = log(c(0.1, 0.4)),
    n1 = 3L, n2 = 2L, kh0_1 = 1, kh0_2 = 1
  )
  expect_equal(phontrast:::.bhatt_mc(mc, loo = FALSE), 1)
  v <- phontrast:::.kernel_family_vector(jsd = 0.25, overlap = 0.4, bc = 0.8)
  expect_equal(unname(v["js_distance"]), 0.5)
  expect_equal(unname(v["total_variation"]), 0.6)
  expect_equal(unname(v["bhatt_kde_dist"]), -log(0.8))
  expect_equal(unname(v["hellinger"]), sqrt(0.2))
  expect_identical(names(v), phontrast:::.kernel_family_columns())
})

test_that("matched-kernel Bhattacharyya tracks the closed form on Gaussian data", {
  d <- .family_fixture(n = 300, seed = 3)
  d <- d[d$speaker == "s01", ]
  feats <- c("f1", "f2")
  p <- phontrast(d, feats, "vowel", metrics = c("bhattacharyya", "bhattacharyya_kde"))
  expect_true(p$bhatt_kde_affinity > 0 && p$bhatt_kde_affinity < 1)
  expect_equal(p$bhatt_kde_dist, -log(p$bhatt_kde_affinity))
  expect_equal(p$hellinger, sqrt(1 - p$bhatt_kde_affinity))
  # same quantity, two estimators: they agree to within kernel bias
  expect_lt(abs(p$bhatt_kde_affinity - p$bhatt_affinity), 0.08)

  # under the mvnorm backend the kernel-family columns are Monte-Carlo
  # estimates between the fitted Gaussians, so they approach the closed form
  pm <- phontrast(d, feats, "vowel", metrics = c("bhattacharyya", "bhattacharyya_kde"),
                  density = "mvnorm", mc_n = 40000, eval_seed = 9)
  expect_lt(abs(pm$bhatt_kde_affinity - pm$bhatt_affinity), 0.02)

  # the legacy estimator still returns a full, finite family
  pl <- phontrast(d, feats, "vowel", method = "legacy",
                  metrics = c("jsd", "overlap", "tv", "bhattacharyya_kde"))
  expect_true(all(is.finite(unlist(pl[, c("jsd", "percent_overlap", "total_variation",
                                           "bhatt_kde_affinity", "hellinger")]))))
})

test_that("euclidean standardizes by the pooled two-category SD", {
  x <- c(0, 0, 0, 0, 2, 2, 2, 2)
  d <- data.frame(cat = rep(c("a", "b"), each = 4), f1 = x, f2 = 2 * x,
                  stringsAsFactors = FALSE)
  got <- phontrast:::.euclidean_means(d, c("f1", "f2"), "cat")
  expect_equal(got, sqrt(2) * (2 / stats::sd(x)))
  # bounded by 2 * sqrt(d) for equal category sizes
  expect_lt(got, 2 * sqrt(2))
  # scale invariance of the standardized distance
  d2 <- d
  d2$f1 <- 1000 * d2$f1
  expect_equal(phontrast:::.euclidean_means(d2, c("f1", "f2"), "cat"), got)
  expect_error(
    phontrast:::.euclidean_means(data.frame(cat = c("a", "a", "b", "b"), f1 = c(1, 1, 1, 1)),
                                 "f1", "cat"),
    "spread"
  )

  dd <- .family_fixture(n = 40)
  p <- phontrast(dd, c("f1", "f2"), "vowel", group_col = "speaker",
                 metrics = c("euclidean", "mahalanobis"))
  expect_true(all(c("euclidean_dist", "mahalanobis_dist") %in% names(p)))
  s01 <- dd[dd$speaker == "s01", ]
  expect_equal(p$euclidean_dist[p$group == "s01"],
               phontrast:::.euclidean_means(s01, c("f1", "f2"), "vowel"))
})

test_that("long output and bootstrap handle the new metrics", {
  d <- .family_fixture(n = 40)
  long <- phontrast(d, "f1", "vowel", output = "long",
                    metrics = c("overlap", "tv", "bhattacharyya_kde", "euclidean"))
  expect_true(all(c("Total variation", "Hellinger distance",
                    "Bhattacharyya affinity (kernel)",
                    "Euclidean distance of standardized means") %in% long$metric))
  tv <- long[long$metric == "Total variation", ]
  ov <- long[long$metric == "Percent overlap", ]
  expect_equal(tv$separation_value, ov$separation_value)
  expect_identical(tv$orientation, "separation")
  expect_true(all(is.finite(long$separation_rank)))

  b <- phontrast(d, "f1", "vowel", metrics = c("tv", "euclidean"),
                 do_boot = TRUE, n_boot = 4, progress = FALSE)
  expect_true(all(c("total_variation_mean", "total_variation_ci_lower",
                    "euclidean_dist_mean", "euclidean_dist_n_boot") %in% names(b)))
  expect_equal(b$euclidean_dist_n_boot, 4L)
})
