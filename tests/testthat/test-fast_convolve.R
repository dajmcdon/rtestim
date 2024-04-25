test_that("fast convolution matches existing version", {
  set.seed(12345)
  y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
  delay_distn <- discretize_gamma(0:100)

  ## old version
  cw <- cumsum(delay_distn)
  convolved_seq <- stats::convolve(y, rev(delay_distn), type = "open")
  convolved_seq <- convolved_seq[seq_along(y)] / cw

  exact_convolved_seq <- function(y, delay) {
    n <- length(delay)
    cmat <- t(toeplitz2(c(rev(delay), rep(0, n)), n, n))
    drop(cmat %*% y) / rowSums(cmat)
  }

  ecs <- exact_convolved_seq(y, delay_distn)
  expect_equal(ecs[2:100], convolved_seq[2:100])

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
  expect_true(all(convolved_seq[1:2] < 0))

  bcs <- delay_calculator(incidence, dist_gamma = gamma_pars)
  expect_true(all(bcs >= 0))

  expect_equal(bcs[10:length], convolved_seq[10:length])

})
