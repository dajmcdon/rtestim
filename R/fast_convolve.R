fast_convolve <- function(y, delay) {
  n <- length(y)
  arg_is_length(n, delay)
  if (n < 20) return(calc_delays(y, delay))
  cw <- cumsum(delay)
  convolved_seq <- stats::convolve(y, rev(delay), type = "open")
  convolved_seq <- convolved_seq[seq_along(y)] / cw
  convolved_seq[1:20] <- calc_delays(y[1:20], delay[1:20])
  convolved_seq
}
