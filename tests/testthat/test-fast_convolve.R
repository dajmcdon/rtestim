test_that("fast convolution matches existing version", {
  set.seed(12345)
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  delay_distn <- discretize_gamma(0:100)

  ## old version
  cw <- cumsum(delay_distn)
  convolved_seq <- stats::convolve(y, rev(delay_distn), type = "open")
  convolved_seq <- convolved_seq[seq_along(y)] / cw
  convolved_seq[1] <- 0 # there's an Inf here

  exact_convolved_seq <- function(y, delay) {
    n <- length(delay)
    cmat <- t(toeplitz2(c(rev(delay), rep(0, n)), n, n))
    drop(cmat %*% y) / rowSums(cmat)
  }

  ecs <- exact_convolved_seq(y, delay_distn)
  ecs[1] <- 0 # there's an NaN here
  expect_equal(ecs, convolved_seq)

  fcs <- fast_convolve(y, delay_distn)
  expect_equal(ecs, fcs)
})

test_that("negatives no longer appear in convolution", {
  length <- 300
  ## piecewise constant Rt:
  Rt <- c(rep(2, floor(2 * length / 5)), rep(0.8, length - floor(2 * length / 5)))
  Mu <- double(length)
  Mu[1] <- 2
  incidence <- double(length)
  w <- double(length)
  ## gamma parameters of measles:
  gamma_pars <- c(14.5963, 1.0208)

  set.seed(317)
  incidence[1] <- Mu[1]
  for (i in 2:length) {
    w[i] <- delay_calculator(incidence[1:i], dist_gamma = gamma_pars)[i]
    Mu[i] <- Rt[i] * w[i]
    incidence[i] <- rpois(1, Mu[i])
  }
  dd <- discretize_gamma(0:(length - 1), gamma_pars[1], gamma_pars[2])
  ## old version
  cw <- cumsum(dd)
  convolved_seq <- stats::convolve(incidence, rev(dd), type = "open")
  convolved_seq <- convolved_seq[seq_along(incidence)] / cw


  bcs <- delay_calculator(incidence, dist_gamma = gamma_pars)
  expect_true(all(bcs >= 0))

  expect_equal(bcs[10:length], convolved_seq[10:length])

  skip_on_ci() # seeds not necessarily reproducible
  expect_true(all(convolved_seq[1:2] < 0))
})

test_that("difficult case of potential negatives still works", {
  set.seed(12345)
  n <- 99
  y <- c(0, 0, 1, 2, 4, 4, 2, 0, 0, 1, rep(0, 20), rpois(69, 4))
  delay <- c(discretize_gamma(0:9), rep(0, n - 10))

  o <- fast_convolve(y, delay) # should be exact
  expect_true(all(o >= 0))

  y <- c(y, 1, 2, 3)
  delay <- c(delay, 0, 0, 0) # now uses fft
  o2 <- fast_convolve(y, delay)
  expect_true(all(o2 >= 0))
  expect_equal(o[1:n], o2[1:n])
})
