test_that("plot error for unseen lambda", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
  expect_error(plot(out, which_lambda = 2))
})


test_that("plot warning if plot_Rt FALSE but use_lambda not NULL", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  cv <- cv_estimate_rt(y, degree = 3, fold = 2, nsol=30)
  expect_warning(plot(cv, plot_Rt = FALSE, which_lambda = 2))
})
