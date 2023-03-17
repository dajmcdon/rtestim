test_that("logarithmic matrix interpolation works", {
  lam1 <- 5.5
  lam_list <- interpolate_lambda(1:10, lam1)

  res1 <- exp(log(21:25) * .5 + log(26:30) * .5)
  expect_equal(res1, interpolate_mat(matrix(1:50, ncol = 10), lam_list, lam1))

  lam2 <- c(5.5, 6.25, 7)
  lam_list <- interpolate_lambda(1:10, lam2)
  res2 <- cbind(unname(res1), exp(log(26:30) * .75 + log(31:35) * .25), 31:35)
  expect_equal(res2, interpolate_mat(matrix(1:50, ncol = 10), lam_list, lam2))
})
