test_that("CV 'runs' with even spaced observation", {
  observed_counts <- c(1:10)
  x <- c(1:10) # even observation
  lambda <- c(1,2,3)
  cv_result <- cv_estimate_rt(observed_counts = observed_counts,
                              x = x,
                              lambda = lambda)
  expect_equal(length(cv_result$cv_scores),  3)
})


test_that("CV 'runs' on uneven spaced observation", {
  observed_counts <- c(1, 3, 4, 6, 8, 10, 11, 12, 15, 16)
  x <- observed_counts # spacing "happens" to corresponds to the observed_counts
  lambda <- c(1,2,3)
  cv_result <- cv_estimate_rt(observed_counts = observed_counts,
                              x = x,
                              lambda = lambda)
  expect_equal(length(cv_result$cv_scores) == 0)
})
