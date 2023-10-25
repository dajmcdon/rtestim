test_that("cv passes correct parameters to inner solvers (when lamdba is not
          specified)", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  cv <- cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 30)
  mod <- estimate_rt(y, korder = 3, nsol = 30)
  expect_identical(cv$lambda, mod$lambda)
})

test_that("test CV returns a warning message when max iteration is insufficient", {
  set.seed(1001)
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  expect_warning(cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 5, maxiter = 1e2L))
  expect_no_warning(cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 5, maxiter = 1e4L))
})
