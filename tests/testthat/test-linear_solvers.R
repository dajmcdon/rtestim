test_that("warn if linear solver is not 1 or 2", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  expect_warning(estimate_rt(y, linear_solver = 3))
})

test_that("same results using two solvers", {
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  rt_kf <- estimate_rt(y, linear_solver = 1)$Rt
  rt_qr <- estimate_rt(y, linear_solver = 2)$Rt
  expect_equal(rt_kf, rt_qr)
})
