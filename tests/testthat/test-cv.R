test_that("even spaced observation", {
  observed_counts <- c(1:10)
  x <- c(1:10) # even observation
  lambda <- c(1,2,3)
  cv_result <- cv_estimate_rt(observed_counts = observed_counts,
                              x = x,
                              lambda = lambda)
  expect_equal(length(cv_result$cv_scores),  3)
})


test_that("uneven spaced observation", {
  observed_counts <- c(1, 3, 4, 6, 8, 10, 11, 12, 15, 16)
  x <- observed_counts # spacing "happens" to corresponds to the observed_counts
  lambda <- c(1,2,3)
  cv_result <- cv_estimate_rt(observed_counts = observed_counts,
                              x = x,
                              lambda = lambda)
  expect_equal(length(cv_result$cv_scores) == 0)
})


test_that("fold index helper function works", {
  n <- 10
  fold <- 4
  expected_idx <- c(0, 1, 2, 3, 4, 1, 2, 3, 4, 0)
  expect_equal(expected_idx, fold_calculator(n, fold))
})


test_that("correctly predict Rt at hold out fold for evenly spaced observation", {
  n <- 10
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) # not needed, for visualization only
  foldid <- c(0, 1, 2, 1, 2, 1, 2, 1, 2, 0)

  train_idx <- c(1, 3, 5, 7, 9, 10)
  test_idx <- c(2, 4, 6, 8)

  train_x <- c(1, 3, 5, 7, 9, 10)
  test_x <- c(2, 4, 6, 8)

  r1 <- matrix(nrow=6, ncol= 2)
  r1[,1] <- c(1.1, 1.3, 1.5, 1.7, 1.9, 2.0)
  r1[,2] <- c(1.5, 1.9, 2.3, 2.7, 3.1, 3.3)

  r1_true <- matrix(nrow=4, ncol = 2)
  r1_true[ ,1] <- c(1.2, 1.4, 1.6, 1.8)
  r1_true[ ,2 ] <- c(1.7, 2.1, 2.5, 2.9)

  pred_r1 <- pred_kth_rt(r1, n, train_idx, test_idx, train_x, test_x)

  expect_equal(pred_r1, r1_true)
})
