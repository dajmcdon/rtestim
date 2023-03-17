test_that("lambda interpolation works", {
  lam <- 1:10
  res_exact <- list(left = 2, right = 2, frac = 1)
  expect_identical(interpolate_lambda(lam, 2), res_exact)
  res_mid <- list(left = 1, right = 2, frac = 0.5)
  expect_identical(interpolate_lambda(lam, 1.5), res_mid)
  res_two <- list(left = c(1, 3), right = c(2, 4), frac = c(0.5, 0.25))
  expect_equal(interpolate_lambda(lam, c(1.5, 3.75)), res_two)

  expect_warning(interpolate_lambda(lam, 0))
  expect_warning(rtb <- interpolate_lambda(lam, 11))

  res_too_big <- list(left = 10, right = 10, frac = 1)
  expect_identical(rtb, res_too_big)
})
