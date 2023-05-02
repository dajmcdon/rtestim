test_that("run linearized admm with RcppEigen integration", {
  set.seed(12345)
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  n <- length(y)
  x <- as.double(1:n)
  ord <- 1L
  w <- delay_calculator(y)

  mod_linear <- rtestim_path(
    algo = 1,
    y = y,
    x = x,
    w = w,
    korder = ord,
    maxiter = 1e4L,
    lambda = .2,
    tolerance = 1e-5
  )

  mod_prox <- rtestim_path(
    algo = 2,
    y = y,
    x = x,
    w = w,
    korder = ord,
    maxiter = 1e4L,
    lambda = .2,
    tolerance = 1e-5
  )

  expect_true(mod_linear$niter < 1e4L)
  expect_true(mod_prox$niter < 1e4L)

  #expect_equal(drop(mod_linear$Rt), drop(mod_prox$Rt), tolerance = 1e-2)
  expect_true(mean(abs(drop(mod_linear$Rt) - drop(mod_prox$Rt))) <= .1)
})
