test_that("Warn if observed_counts start at 0", {
  set.seed(12345)
  y <- c(0, 0, rpois(18, 3))
  expect_warning(o0 <- estimate_rt(y, nsol = 20))
})


test_that("Error if observed_counts has neg values", {
  y <- c(0, -2, -3, rpois(20, 3))
  expect_error(estimate_rt(observed_counts = y))
})

test_that("estimate_rt works with uneven spacing", {
  y <- c(1:10)
  x <- c(1, 2, 4, 6, 7, 9, 10, 15, 18, 20)
  expect_no_error(estimate_rt(observed_counts = y, x = x, korder = 2, lambda = 1))
})
