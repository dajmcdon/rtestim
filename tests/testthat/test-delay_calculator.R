test_that("delay calculator with even spacing works", {
  n <- 5
  counts <- 1:n
  w <- discretize_gamma(1:n)
  delay_calculator_output <- delay_calculator(observed_counts = counts)
  manual_output <- rep(0, (n - 1))
  w <- discretize_gamma(1:n)
  cw <- cumsum(w)
  for (j in 2:n) {
    manual_output[j - 1] <- sum(w[1:(j - 1)] * counts[(j - 1):1] / cw[(j - 1)])
  }
  manual_output <- c(manual_output[1], manual_output)
  expect_equal(manual_output, delay_calculator_output)
})

test_that("weighted_past_counts has same length as observed_counts", {
  y <- c(1, 3, 4, 6, 7, 13)
  weighted_past_counts <- delay_calculator(y, x = c(1:6))
  expect_true(length(y) == length(weighted_past_counts))
})
