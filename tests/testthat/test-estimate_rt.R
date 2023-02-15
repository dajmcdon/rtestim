test_that("works if observed_counts contain 0", {
  y <- c(rep(0, 5), rpois(20, 3))
  expect_no_error(estimate_rt(y))
})
