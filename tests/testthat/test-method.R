test_that("plot error for unseen lambda", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
  expect_error(plot(out, which_lambda = 2))
})
