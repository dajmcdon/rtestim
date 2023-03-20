set.seed(1116)
n <- 100
# settings
M <- 1e3L
x <- 1:n
w <- rep(1, n)
theta <- double(n)
lambda <- .1
rho <- .1
M_inner <- 3L
iter <- integer(1)
alpha <- .5
gamma <- .9
tol <- 1e-4


# k=0
ord <- 0
D <- buildDx_tilde(n, ord, x)
z <- double(n - ord)
u <- double(n - ord)
mu <- 2 * 4^ord
for (i in 1:3) {
  y <- rpois(n, 10)
  lin_mod <- linear_admm_testing(
    M, y, x, w, n, ord, theta, z, u, lambda, rho, mu, tol, iter
  )
  pro_mod <- prox_newton_testing(
    M, M_inner, n, ord, y, x, w, theta, z,
    u, lambda, rho, alpha, gamma, D, tol, iter
  )
  # print(cbind(head(exp(lin_mod$theta)),head(exp(pro_mod$theta))))
  test_that("test two algorithms equal for k=0", {
    expect_equal(exp(lin_mod$theta[, 1]), exp(lin_mod$theta[, 1]), tolerance = .1)
  })
}

# k=1
ord <- 1
D <- buildDx_tilde(n, ord, x)
z <- double(n - ord)
u <- double(n - ord)
mu <- 2 * 4^ord
for (i in 1:3) {
  y <- rpois(n, 10)
  lin_mod <- linear_admm_testing(
    M, y, x, w, n, ord, theta, z, u, lambda, rho, mu, tol, iter
  )
  pro_mod <- prox_newton_testing(
    M, M_inner, n, ord, y, x, w, theta, z,
    u, lambda, rho, alpha, gamma, D, tol, iter
  )
  test_that("test two algorithms equal for k=1", {
    expect_equal(exp(lin_mod$theta[, 1]), exp(lin_mod$theta[, 1]), tolerance = .1)
  })
}

# k=2
ord <- 2
D <- buildDx_tilde(n, ord, x)
z <- double(n - ord)
u <- double(n - ord)
mu <- 2 * 4^ord
for (i in 1:3) {
  y <- rpois(n, 10)
  lin_mod <- linear_admm_testing(
    M, y, x, w, n, ord, theta, z, u, lambda, rho, mu, tol, iter
  )
  pro_mod <- prox_newton_testing(
    M, M_inner, n, ord, y, x, w, theta, z, u, lambda, rho,
    alpha, gamma, D, tol, iter
  )
  test_that("test two algorithms equal for k=2", {
    expect_equal(exp(lin_mod$theta[, 1]), exp(lin_mod$theta[, 1]), tolerance = .1)
  })
}
