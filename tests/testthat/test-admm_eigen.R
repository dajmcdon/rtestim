test_that("run linearized admm with RcppEigen integration", {
  set.seed(12345)
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  n <- length(y)
  x <- as.double(1:n)
  tt_eigen <- admm_eigen_testing(
    M = 1e4L,
    y = y,
    x = x,
    w = delay_calculator(y),
    n = n,
    ord = 4L,
    theta = double(n),
    z = double(n - 4L),
    u = double(n - 4L),
    lambda = 2.0,
    rho = 2.0,
    mu = 2 * 4^4,
    tol = 1e-4,
    iter = 10
  )
})
