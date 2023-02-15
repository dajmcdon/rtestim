test_that("delay calculator with even spacing works", {
  n <- 5
  counts <- 1:n
  w <- discretize_gamma(1:n)
  delay_calculator_output <- delay_calculator(observed_counts = counts)
  manual_output <- rep(0, (n-1))
  w <- discretize_gamma(1:n)
  cw <- cumsum(w)
  for (j in 2:n){
    manual_output[j-1] <- sum(w[1:(j-1)] * counts[(j-1):1]/cw[(j-1)])
  }
  manual_output <- c(manual_output[1], manual_output)
  expect_equal(manual_output, delay_calculator_output)
})


test_that("fill_case_counts work as intended", {
  # Uneven integer spacing, at most 1 gap
  x1 <- c(1,3,4,5,7,9)
  o1 <- 2*x1
  filled_o1 <- fill_case_counts(x1, o1)
  o1_true <- 2*c(1:9)
  expect_equal(filled_o1, o1_true)

  # Uneven integer spacing, works with any number of gaps
  x2 <- c(1,3,4,5,7,9,16)
  o2 <- 2*x2
  filled_o2 <- fill_case_counts(x2, o2)
  o2_true <- 2*c(1:16)
  expect_equal(filled_o2, o2_true)

  # Uneven double spacing, works with any number of gaps
  x3 <- c(1, 1.6, 2.2, 3.4, 4, 5.2, 7, 8.2)
  o3 <- 4*x3
  filled_o3 <- fill_case_counts(x3, o3)
  o3_true <- 4*seq(1, 8.2, 0.6)
  expect_equal(filled_o3, o3_true)
})


test_that("weighted_past_counts[1] is reasonable when observed_counts[1] is 0", {
  y <- c(0, 0, rpois(5, lambda = 3))
  delay_calculator_output <- delay_calculator(y)
  expect_false(delay_calculator_output[1] < 1e-3)
})
