

test_that("fold index helper function works", {
  n <- 10
  fold <- 4
  expected_idx <- c(0, 1, 2, 3, 4, 1, 2, 3, 4, 0)
  expect_equal(expected_idx, fold_calculator(n, fold))
})
