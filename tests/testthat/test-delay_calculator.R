test_that("delay calculator with even spacing works", {
  counts <- c(1:5, 5:1)
  n <- length(counts)
  delay_calculator_output <- delay_calculator(observed_counts = counts)
  manual_output <- double(n)
  w <- discretize_gamma(0:(n - 1))
  cw <- cumsum(w)
  for (j in 1:n) manual_output[j] <- sum(w[1:j] * counts[j:1] / cw[j])
  manual_output[1] <- manual_output[2] # first value is 0 / 0
  expect_equal(manual_output, delay_calculator_output)

  w <- c(1, 2, 3, 2, 1, 1)
  delay_calculator_output <- delay_calculator(counts, delay_distn = w)
  w <- c(w, rep(0, n - length(w)))
  cw <- cumsum(w)
  manual_output <- double(n)
  for (j in 1:n) manual_output[j] <- sum(w[1:(j)] * counts[j:1]) / cw[j]
  expect_equal(manual_output, delay_calculator_output)
})

test_that("weighted_past_counts has same length as observed_counts", {
  y <- c(1, 3, 4, 6, 7, 13)
  weighted_past_counts <- delay_calculator(y, x = 1:6)
  expect_true(length(y) == length(weighted_past_counts))
})

test_that("delay calculator errors out when necessary", {
  expect_error(delay_calculator(1:10, dist_gamma = 1))
  expect_error(delay_calculator(1:10, dist_gamma = c(1, -1)))
  expect_error(delay_calculator(1:10, delay_distn = -1:10 / 10))
  expect_error(delay_calculator(1:10, x = 0:10))
  expect_error(delay_calculator(1:10, x = c(1, NA, 3:10)))
  expect_error(delay_calculator(1:10, x = c(1, 1, 3:10)))
  expect_error(delay_calculator(1:10, x = c(3, 2, 4:11)))
  expect_error(delay_calculator(1:10, xout = 0:10))
  expect_error(delay_calculator(1:10, xout = c(1, NA, 3:10)))
  expect_error(delay_calculator(1:10, xout = c(1, 1, 3:10)))
  expect_error(delay_calculator(1:10, xout = c(3, 2, 4:11)))
  expect_error(delay_calculator(1:10, xout = 0:11))
  expect_error(delay_calculator(1:10, x = 1:10, xout = 2:20 / 2))
})

test_that("delay calculator correctly handles periodicity", {
  xd <- seq(as.Date("2020-01-01"), as.Date("2020-05-01"), by = 1)
  x <- seq_along(xd)
  xw <- xd[seq(1, max(x), by = 7)]
  yd <- rep(c(1:7, 8:2), length.out = length(x))
  yw <- yd[seq(1, max(x), by = 7)]

  expect_error(delay_calculator(yw, xw, xout = xd))
  expect_error(delay_calculator(yw, xw, xout = xd[xd <= max(xw)]))
  expect_error(delay_calculator(
    yw, xw,
    delay_distn_periodicity = "1 week", xout = xd[xd <= max(xw)]
  ))
  dweekly <- delay_calculator(
    yw, xw,
    delay_distn_periodicity = 1, xout = xd[xd <= max(xw)]
  )
  dweekly_text <- delay_calculator(
    yw, xw,
    delay_distn_periodicity = "1 day", xout = xd[xd <= max(xw)]
  )
  expect_identical(dweekly, dweekly_text)
  ddaily <- delay_calculator(yd[xd <= max(xw)], xd[xd <= max(xw)])
  expect_identical(dweekly, ddaily)

  delay_distn <- 7:1
  ddn <- c(delay_distn, rep(0, length(yd) - length(delay_distn)))
  dwn <- c(delay_distn, rep(0, length(yw) - length(delay_distn)))
  expect_equal(
    delay_calculator(yd, xd, delay_distn = delay_distn),
    stats::convolve(yd, rev(ddn), type = "open")[seq_along(yd)] / cumsum(ddn)
  )
  expect_equal(
    delay_calculator(yw, xw, delay_distn = delay_distn),
    stats::convolve(yw, rev(dwn), type = "open")[seq_along(yw)] / cumsum(dwn)
  )
  expect_error(delay_calculator(
    yw, xw,
    delay_distn = delay_distn, delay_distn_periodicity = 7,
    xout = xd[xd <= max(xw)]
  ))
  s <- xd <= max(xw)
  expect_equal(
    delay_calculator(
      yw, xw,
      delay_distn = delay_distn, delay_distn_periodicity = 1,
      xout = xd[s]
    ),
    stats::convolve(yd[s], rev(ddn[s]), type = "open")[seq(sum(s))] / cumsum(ddn[s])
  )

  cw <- stats::convolve(yd[s], rev(ddn[s]), type = "open")[seq(sum(s))] / cumsum(ddn[s])
  expect_equal(
    delay_calculator(yw, xw, delay_distn = delay_distn, delay_distn_periodicity = 1),
    cw[xw - min(xw) + 1]
  )
})

test_that("delay calculator accommodates alternative delays", {
  y <- cancovid$incident_cases[1:100]
  wpc <- delay_calculator(y)

  dist_gamma <- discretize_gamma(0:99)

  expect_equal(wpc, delay_calculator(y, delay_distn = dist_gamma))

  # wrong dimensions
  d_mat <- matrix(1, nrow = 200, ncol = 200)
  d_mat[upper.tri(d_mat)] <- 0
  expect_error(delay_calculator(y, delay_distn = d_mat))

  # not lower triangular
  d_mat <- matrix(1, nrow = 100, ncol = 100)
  expect_error(delay_calculator(y, delay_distn = d_mat))

  d_mat[upper.tri(d_mat, diag = TRUE)] <- 0
  expect_error(delay_calculator(y, delay_distn = d_mat))

  d_mat[1,1] <- 1
  roll_avg <- cumsum(y) / seq_along(y)
  roll_avg <- c(roll_avg[1], roll_avg[-100])
  expect_equal(roll_avg, delay_calculator(y, delay_distn = d_mat))
  expect_equal(
    roll_avg,
    delay_calculator(y, delay_distn = drop0(as(d_mat, "CsparseMatrix")))
  )

  d_mat <- matrix(0, 100, 100)
  d_mat[1,1] <- 1
  for (i in 2:100) d_mat[i,1:i] <- rev(dist_gamma[1:i])
  d_mat <- drop0(as(d_mat, "CsparseMatrix"))
  d_mat <- d_mat / rowSums(d_mat)
  expect_equal(wpc, delay_calculator(y, delay_distn = d_mat))
})
