test_that("plot error for unseen lambda", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
  expect_error(plot(out, which_lambda = 2))
})


test_that("plot warning if plot_Rt FALSE but use_lambda not NULL", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  cv <- cv_estimate_rt(y, degree = 3, nfold = 2, nsol=30)
  expect_warning(plot(cv, plot_Rt = FALSE, which_lambda = 2))
})

test_that("match_lambda return abort and warn", {
  expect_warning(w <- match_lambda(which_lambda = c(1,2,3), lambda = c(2,7,3)))
  expect_equal(w, c(1,3))
  expect_error(match_lambda(which_lambda = c(4,5,6), lambda = c(1,2)))

  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  cv <- cv_estimate_rt(y, degree = 3, nfold = 2, nsol=30)
  expect_warning(plot(cv, plot_Rt = TRUE, which_lambda = c(cv$lambda[1], 3)))
})
