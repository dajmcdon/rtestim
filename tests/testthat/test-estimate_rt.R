test_that("Report error if observed_counts start at 0", {
  y <- c(0, rpois(20, 3))
  expect_error(estimate_rt(y))
})


test_that("Report error if observed_counts has non-neg values 0", {
  y <- c(0, -2, -3, rpois(20, 3))
  expect_error(estimate_rt(y))
})

test_that("estimate_rt works with uneven spacing", {
  y <- c(1:10)
  x <- c(1,2,4,6,7,9,10,15,18,20)
  expect_no_error(estimate_rt(y, x = x, degree = 2, lambda = 1))
})
