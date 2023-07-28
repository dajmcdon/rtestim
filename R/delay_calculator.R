#' Calculate the total infectiousness at each observed time point.
#'
#' The total infectiousness at each observed time point is calculated
#' by \eqn{\sum_{a=1}^t I_{t-s}w_s}, where \eqn{I} denotes the vector containing
#' all observed case counts, and \eqn{w} denotes the serial interval
#' distribution. The serial interval distribution expresses the probability
#' of the symptom onset of a secondary infection occurred a given
#' number of days after the primary infection
#'
#' @inheritParams estimate_rt
#' @param output_partial_seq Logical. By default, this function returns the
#'   convolved weight vector with observed cases at the original `x` values.
#'   However, when `x` is irregular, setting this to `FALSE` will result
#'   in a result at the interpolated `x` sequence.
#'
#' @return A vector containing the total infectiousness at each
#'   observed time point
#' @export
#'
#' @examples
#' delay_calculator(c(3,2,5,3,1), dist_gamma = c(2.5, 2.5))
delay_calculator <- function(
    observed_counts,
    x = NULL,
    dist_gamma = c(2.5, 2.5),
    delay_distn = NULL,
    output_partial_seq = TRUE) {

  arg_is_length(2, dist_gamma)
  arg_is_lgl_scalar(output_partial_seq)
  arg_is_positive(dist_gamma)
  arg_is_positive(delay_distn, allow_null = TRUE)
  n <- length(observed_counts)
  arg_is_length(n, x, allow_null = TRUE)
  if (is.null(x)) x <- 1:n
  else {
    if (any(is.na(x))) cli_abort("x may not contain missing values.")
    if (is.unsorted(x, strictly = TRUE))
      cli_abort("x must be sorted and contain no duplicates.")
  }

  if (inherits(x, "Date")) x <- as.numeric(x)
  arg_is_numeric(x)
  if (!is.null(delay_distn)) delay_distn <- delay_distn / sum(delay_distn)
  regular <- vctrs::vec_unique_count(diff(x)) == 1L
  if (regular) xout <- x
  else xout <- seq(from = min(x), to = max(x), by = min(diff(x)))

  if (is.null(delay_distn)) {
    delay_distn <- discretize_gamma(xout, dist_gamma[1], dist_gamma[2])
  } else {
    if (length(delay_distn) > length(xout)) {
      cli_abort(
        "User specified `delay_distn` must have no more than {length(xout)} elements."
      )
    }
    # pad the tail with zero if too short
    delay_distn <- c(delay_distn, rep(0, length(xout) - length(delay_distn)))
  }

  y <- stats::approx(x, observed_counts, xout = xout)$y
  cw <- cumsum(delay_distn)

  convolved_seq <- stats::convolve(y, rev(delay_distn), type = "open")
  convolved_seq <- convolved_seq[seq_along(xout)] / cw
  convolved_seq <- c(convolved_seq[1], convolved_seq[-length(convolved_seq)])
  if (!regular && output_partial_seq)
    convolved_seq <- convolved_seq[xout %in% x]
  convolved_seq
}


