test_that("gcd works", {
  expect_identical(gcd(c(1, NA)), NA)
  expect_identical(gcd(c(1, NA), TRUE), 1L)
  expect_identical(gcd(c(1, 2, 4)), 1L)
  expect_error(gcd(1.3))
  expect_identical(gcd(c(2, 4, 6)), 2L)
  expect_identical(gcd(c(2, 3, 7)), 1L)
})
