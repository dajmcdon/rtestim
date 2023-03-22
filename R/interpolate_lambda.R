interpolate_from_ref_sequence <- function(ref_sequence, sequence) {
  r <- range(ref_sequence)
  rr <- range(sequence)
  if (rr[1] < r[1]) {
    cli::cli_warn(
      c("Supplied `ref_sequence` has values smaller than those used to estimate Rt originally.",
        i = "You may want to refit at these values.",
        i = "Using the smallest available `ref_sequence`'s available."))
  }
  if (rr[2] > r[2]) {
    cli::cli_warn(
      c("Supplied `ref_sequence` has values larger than those used to estimate Rt originally.",
        i = "You may want to refit at these values.",
        i = "Using the largest available `ref_sequence`'s available."))
  }
  if (length(ref_sequence) == 1) {# degenerate case of only one ref_sequence
    nums <- length(sequence)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    k <- length(ref_sequence)
    sfrac <- (ref_sequence[1] - sequence) / (ref_sequence[1] - ref_sequence[k])
    ref_sequence <- (ref_sequence[1] - ref_sequence) / (ref_sequence[1] - ref_sequence[k])
    sfrac[sfrac < min(ref_sequence)] <- min(ref_sequence)
    sfrac[sfrac > max(ref_sequence)] <- max(ref_sequence)
    coord <- stats::approx(ref_sequence, seq(ref_sequence), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - ref_sequence[right]) / (ref_sequence[left] - ref_sequence[right])
    sfrac[left == right] <- 1
    sfrac[abs(ref_sequence[left] - ref_sequence[right]) < .Machine$double.eps] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}
