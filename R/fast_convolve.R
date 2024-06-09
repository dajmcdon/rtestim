fast_convolve <- function(y, delay) {
  n <- length(y)
  arg_is_length(n, delay)
  if (n < 100) return(calc_delays(y, delay)) # this is exact, normalizes

  # Has potential errors due to fft / ifft
  convolved_seq <- stats::convolve(y, rev(delay), type = "open")
  convolved_seq <- convolved_seq[seq_along(y)]

  # now we find the support and 0 out the rest
  eps <- .Machine$double.eps
  if (any(delay < eps) || any(y < eps)) {
    idy <- as.double(y > eps)
    idd <- as.double(delay > eps)
    id0 <- stats::convolve(idy, rev(idd), type = "open")[seq_along(y)]
    convolved_seq[id0 < 0.5] <- 0
    convolved_seq[convolved_seq < eps] <- 0
  }
  # normalize
  cw <- cumsum(delay)
  wcs <- convolved_seq / cw
  wcs[cw < eps] <- 0
  wcs
}
