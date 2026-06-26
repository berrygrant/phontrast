test_that("discrete JSD handles zeros correctly", {
  expect_equal(jsd(c(1, 0), c(1, 0)), 0)
  expect_equal(jsd(c(1, 0), c(0, 1)), 1)
  expect_error(jsd(c(1, -1), c(1, 1)), "finite, non-negative")
})

test_that("classical grouped metrics drop incomplete rows consistently", {
  set.seed(1)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 20),
    category = rep(rep(c("a", "b"), each = 10), 2),
    f1 = rnorm(40),
    f2 = rnorm(40)
  )
  data$f1[1] <- NA_real_

  pillai <- speaker_pillai(
    data = data,
    group_col = "speaker",
    category_col = "category",
    features = c("f1", "f2"),
    min_tokens = 10
  )
  bhatt <- estimate_bhatt(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 10
  )

  expect_equal(nrow(pillai), 2)
  expect_equal(nrow(bhatt), 2)
  expect_equal(pillai$n_tokens[pillai$group == "s1"], 19)
  expect_equal(bhatt$n_tokens[bhatt$group == "s1"], 19)
})

test_that("KDE JSD wrappers return valid global and grouped estimates", {
  set.seed(2)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )
  data$f2[1] <- NA_real_

  point <- jsd_kde_nd(data, c("f1", "f2"), "category")
  global <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    min_tokens = 20
  )
  grouped <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20
  )

  expect_true(is.finite(point))
  expect_true(point >= 0 && point <= 1)
  expect_equal(global$scope, "global")
  expect_equal(global$n_tokens, 159)
  expect_equal(nrow(grouped), 2)
  expect_true(all(is.finite(grouped$jsd_point)))
})

test_that("bootstrap JSD summaries preserve successful replicate counts", {
  set.seed(3)
  data <- data.frame(
    speaker = rep("s1", 80),
    category = rep(c("a", "b"), each = 40),
    f1 = c(rnorm(40, 0), rnorm(40, 1)),
    f2 = c(rnorm(40, 0), rnorm(40, 1))
  )

  out <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 3,
    min_tokens = 20
  )

  expect_equal(out$n_boot, 3)
  expect_true(is.finite(out$jsd_mean))
  expect_true(is.finite(out$jsd_low))
  expect_true(is.finite(out$jsd_high))
})

test_that("one-dimensional KDE metrics are supported", {
  set.seed(4)
  data <- data.frame(
    category = rep(c("a", "b"), each = 40),
    f1 = c(rnorm(40, 0), rnorm(40, 1))
  )

  jsd_val <- jsd_kde_nd(data, "f1", "category")
  overlap <- percent_overlap_kde(data, "f1", "category")

  expect_true(is.finite(jsd_val))
  expect_true(jsd_val >= 0 && jsd_val <= 1)
  expect_true(is.finite(overlap))
  expect_true(overlap >= 0 && overlap <= 1)
})

test_that("global JSD wrapper forwards KDE controls", {
  set.seed(5)
  data <- data.frame(
    category = rep(c("a", "b"), each = 60),
    f1 = c(rnorm(60, 0), rnorm(60, 1)),
    f2 = c(rnorm(60, 0), rnorm(60, 1))
  )

  direct <- jsd_kde_nd(data, c("f1", "f2"), "category", bw = "Hpi.diag")
  wrapped <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    min_tokens = 10,
    bw = "Hpi.diag"
  )

  expect_equal(wrapped$jsd_point, direct)
})

test_that("grouped bootstrap JSD respects confidence level", {
  set.seed(6)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = rnorm(160),
    f2 = rnorm(160)
  )

  set.seed(7)
  ci50 <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 12,
    min_tokens = 20,
    conf_level = 0.50
  )
  set.seed(7)
  ci95 <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 12,
    min_tokens = 20,
    conf_level = 0.95
  )

  width50 <- ci50$jsd_high - ci50$jsd_low
  width95 <- ci95$jsd_high - ci95$jsd_low
  expect_true(all(width50 <= width95))
})

test_that("grouped metrics return shaped empty outputs when all groups are filtered", {
  set.seed(8)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 20),
    category = rep(rep(c("a", "b"), each = 10), 2),
    f1 = rnorm(40),
    f2 = rnorm(40)
  )

  pillai <- estimate_pillai(data, c("f1", "f2"), "category", "speaker", min_tokens = 100)
  bhatt <- estimate_bhatt(data, c("f1", "f2"), "category", "speaker", min_tokens = 100)
  overlap <- estimate_overlap(data, c("f1", "f2"), "category", "speaker", min_tokens = 100)

  expect_equal(nrow(pillai), 0)
  expect_equal(names(pillai), c("scope", "group", "n_tokens", "pillai", "p_value"))
  expect_equal(nrow(bhatt), 0)
  expect_equal(names(bhatt), c("scope", "group", "n_tokens", "bhatt_dist", "bhatt_affinity"))
  expect_equal(nrow(overlap), 0)
  expect_equal(names(overlap), c("scope", "group", "n_tokens", "overlap"))
})

test_that("prepare_jsd_beta rejects invalid JSD values", {
  ok <- prepare_jsd_beta(data.frame(jsd_mean = c(0, 0.5, 1, NA)), eps = 0.01)
  expect_equal(ok$jsd_beta, c(0.01, 0.5, 0.99, NA))
  expect_error(
    prepare_jsd_beta(data.frame(jsd_mean = c(-0.1, 0.2, Inf))),
    "values in \\[0, 1\\]"
  )
})

test_that("small KDE samples fail with a clear data-size message", {
  set.seed(9)
  data <- data.frame(
    category = rep(c("a", "b"), each = 2),
    f1 = rnorm(4),
    f2 = rnorm(4)
  )

  expect_error(
    jsd_kde_nd(data, c("f1", "f2"), "category"),
    "at least 3 finite observations"
  )
})

test_that("compare_overlap_metrics returns wide and long comparisons", {
  set.seed(10)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )

  wide <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20
  )
  long <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20,
    output = "long"
  )

  expect_equal(nrow(wide), 2)
  expect_true(all(c(
    "pillai", "bhatt_dist", "bhatt_affinity", "jsd", "js_distance",
    "mahalanobis_dist", "percent_overlap"
  ) %in% names(wide)))
  expect_equal(nrow(long), 14)
  expect_true(all(c(
    "metric", "estimate", "orientation", "separation_value", "separation_rank"
  ) %in% names(long)))
})

test_that("hierarchical bootstrap preserves repeated group draws", {
  set.seed(11)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )
  rows_seen <- integer()
  fit_fun <- function(formula, data, ...) {
    rows_seen <<- c(rows_seen, nrow(data))
    stats::lm(formula = formula, data = data, ...)
  }

  set.seed(12)
  out <- hier_boot_jsd_model(
    data = data,
    group_col = "speaker",
    category_col = "category",
    features = c("f1", "f2"),
    formula = jsd_beta ~ 1,
    fit_fun = fit_fun,
    n_outer = 4,
    min_tokens = 20,
    progress = FALSE
  )

  expect_true(nrow(out) > 0)
  expect_true(all(rows_seen == 2))
})

test_that("extract_mfcc validates arguments before optional audio dependencies", {
  expect_error(
    extract_mfcc(data.frame(file = "missing.wav"), "file", numcep = 0),
    "positive integer"
  )
  expect_error(
    extract_mfcc(data.frame(file = "missing.wav"), "file", fs = -1),
    "positive finite"
  )
  expect_error(
    extract_mfcc(data.frame(file = "missing.wav", start = "0", end = 1), "file",
                 start_col = "start", end_col = "end"),
    "numeric column"
  )
})
