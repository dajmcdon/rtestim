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
#'
#' @return A vector of size n-1, containing the total infectiousness at each
#'   observed time point
#' @export
#'
#' @examples
#' delay_calculator(c(3,2,5,3,1), c(2.5, 2.5))
delay_calculator <- function(observed_counts, x = NULL,
                             dist_gamma = c(2.5, 2.5)) {
  if (length(dist_gamma) != 2L)
    cli:cli_abort("dist_gamma must have length 2.")
  if (any(dist_gamma) <= 0)
    cli::cli_abort("dist_gamma must be positive.")
  n <- length(observed_counts)
  if (!is.null(x)) {
    if (any(is.na(x)))
      cli::cli_abort("x may not contain missing values.")
    if (length(x) != n)
      cli::cli_abort("x must be the same length as observed_counts {n}.")
    if (is.unsorted(x, strictly = TRUE))
      cli::cli_abort("x must be sorted and contain no duplicates.")
  } else {
    x <- 1:n
  }

  w <- discretize_gamma(x, dist_gamma[1], dist_gamma[2])
  cw <- cumsum(w)
  regular <- vctrs::vec_unique_count(diff(x)) == 1L
  if (regular) {
    convolved_seq <- convolve(current_counts, rev(w))[1:n] / cw
    return(c(convolved_seq[1], convolved_seq[1:(n - 1)]))
  } else {
    cli::cli_abort("Uh oh. We don't support irregular x yet...")
  }
}
