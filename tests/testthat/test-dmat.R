test_that("building the D matrix works", {
  n <- 10
  x <- 1:n
  Id <- get_D(-1, x)
  D1 <- get_D(0, x)
  D2 <- get_D(1, x)
  D3 <- get_D(2, x)

  # for even spacing, get_D(k-1, x) == get_Dtil(k, x)
  Idm <- sparseMatrix(1:10, 1:10, x = 1, dims = c(10, 10))
  expect_identical(Id, Idm)
  expect_identical(Id, get_Dtil(0, x))
  expect_identical(D1, get_Dtil(1, x))
  expect_identical(D2, get_Dtil(2, x))

  # note that get_Dtil() does NOT include the diagonal weighting,
  # because we don't need it in the ADMM. So for uneven spacing,
  # D1 %*% Dtil(2, x) != D3(x)

  z <- c(1:7, 9:10)
  expect_identical(get_D(2, z), dspline::d_mat(3, z, TRUE)) # note 2 vs 3
  expect_identical(get_Dtil(2, z), dspline::d_mat(2, z, FALSE))
})
