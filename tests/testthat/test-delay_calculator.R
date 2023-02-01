test_that("delay calculator with even spacing works", {
  n <- 5
  counts <- 1:n
  w <- discretize_gamma(1:n)
  delay_calculator_output <- delay_calculator(observed_counts = counts)
  manual_output <- rep(0, (n-1))
  for (j in 2:n){
    manual_output[j-1] <- sum(discretize_gamma(1:(j-1)) * counts[(j-1):1])
  }
  manual_output <- c(manual_output[1], manual_output)
  expect_equal(manual_output, delay_calculator_output)

})
