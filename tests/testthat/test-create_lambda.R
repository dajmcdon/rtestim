test_that("create_lambda works", {
  # should be an increasing sequence
  expect_equal(
    create_lambda_test(double(10), .1, 10, 1e-4, 10),
    10^(seq(log10(10), log10(.1), length.out = 10))
  )
  expect_equal(
    create_lambda_test(10:1, .1, 10, 1e-4, 10),
    as.double(10:1)
  )
  expect_equal(
    create_lambda_test(double(10), -1, 10, 1e-4, 10),
    10^(seq(log10(10), log10(1e-3), length.out = 10))
  )
})
