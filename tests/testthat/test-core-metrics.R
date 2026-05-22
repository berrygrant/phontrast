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
