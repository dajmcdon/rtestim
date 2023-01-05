test_that("building the D matrix works", {
  n <- 10
  x <- 1:n
  D1 <- buildD(10, 0)
  D2 <- buildD(10, 1)
  D3 <- buildD(10, 2)
  expect_identical(D1, buildDx(n, 0, x))
  expect_identical(D2, buildDx(n, 1, x))
  expect_identical(D3, buildDx(n, 2, x))
  expect_identical(D1, buildDx_tilde(n, 1, x))
  expect_identical(D2, buildDx_tilde(n, 2, x))
  expect_identical(D3, buildDx_tilde(n, 3, x))

  z <- c(1:7, 9:10)
  n <- length(z)
  D1 <- buildD(n, 0)
  expect_identical(D1, buildDx(n, 0, z))
  delx <- Matrix::Diagonal(n - 1, 1 / diff(z, 1))
  Dx2_tilde <- delx %*% D1
  Dx2 <- D1[1:(n - 2), 1:(n - 1)] %*% Dx2_tilde
  expect_identical(Dx2_tilde, buildDx_tilde(n, 1, z))
  expect_identical(Dx2, buildDx(n, 1, z))
  delx <- Matrix::Diagonal(n - 2, 2 / diff(z, 2))
  Dx3_tilde <- delx %*% Dx2
  Dx3 <- D1[1:(n - 3), 1:(n - 2)] %*% Dx3_tilde
  expect_identical(Dx3_tilde, buildDx_tilde(n, 2, z))
  expect_identical(Dx3, buildDx(n, 2, z))

})

test_that("building D works with empty x", {
  x = double(0) # no
  expect_identical(buildD(10, 2), buildDx(10, 2, x))
  expect_identical(buildDx_tilde(10, 2, x), buildDx_tilde(10, 2, 1:10))
})
