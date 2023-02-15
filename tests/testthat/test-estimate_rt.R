test_that("works if observed_counts contain 0", {
  y <- c(0,0,rnorm(20, 3))
  expect_no_error(estimate_rt(y))
})
