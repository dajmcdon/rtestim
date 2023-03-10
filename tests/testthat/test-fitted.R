test_that("fitted method for poisson_rt works", {
  set.seed(12345)
  y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
  out <- estimate_rt(y, nsol = 10)
  expect_identical(out$Rt, fitted(out))
  expect_identical(out$Rt[,1], fitted(out, out$lambda[1]))
  expect_warning(ff <- fitted(out, out$lambda[1] - .001))
  expect_identical(out$Rt[,1], ff)

  fff <- suppressWarnings(fitted(out, c(out$lambda[1] - .001, out$lambda[1])))
  expect_identical(out$Rt[,c(1,1)], fff)

  l <- mean(out$lambda[1:2])
  expect_equal(exp(rowMeans(log(out$Rt[,1:2]))), fitted(out, l))
})
