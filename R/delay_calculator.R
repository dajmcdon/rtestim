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
#' delay_calculator(c(3,2,5,3,1), dist_gamma = c(2.5, 2.5))
delay_calculator <- function(observed_counts, x = NULL,
                             dist_gamma = c(2.5, 2.5)) {
  arg_is_length(2, dist_gamma)
  arg_is_positive(dist_gamma)
  n <- length(observed_counts)
  arg_is_length(n, x, allow_null = TRUE)
  if (!is.null(x)) {
    if (any(is.na(x)))
      cli::cli_abort("x may not contain missing values.")
    if (is.unsorted(x, strictly = TRUE))
      cli::cli_abort("x must be sorted and contain no duplicates.")
  } else {
    x <- 1:n
  }

  regular <- vctrs::vec_unique_count(diff(x)) == 1L

  if (!regular) {
    observed_counts <- fill_case_counts(x, observed_counts)
    n <- length(observed_counts)
    x <- 1:n
  }
  w <- discretize_gamma(x, dist_gamma[1], dist_gamma[2])
  cw <- cumsum(w)
  convolved_seq <- stats::convolve(observed_counts, rev(w), type = "open")[1:n] / cw
  c(convolved_seq[1], convolved_seq[1:(n - 1)])
}



#' Interpolate case counts for uneven spaced time points
#'
#' Given observation time points `x` and observed case counts `observed_counts`,
#' this function finds the minimal difference `m` from consecutive `x` and
#' construct a full and even `x` with difference `m`. This function then find
#' the missing index and interpolate the observed case counts using
#' `na.fill` method with argument `fill = "extend"`
#'
#' @inheritParams delay_calculator
#'
#' @return interpolated `observed_counts`
#' @export
#'
#' @examples
#' x1 <- c(1,3,4,5,7,9)
#' o1 <- 2*x1
#' filled_o1 <- fill_case_counts(x1, o1)
#' o1_true <- 2*c(1:9)
#' # o1_true should equal filled_o1
fill_case_counts <- function(x, observed_counts) {
  min_diff <- min(diff(x))
  full_x <- seq(1, max(x), min_diff)
  n_full <- length(full_x)
  full_x <- round(full_x, 5) # avoid numerical issue
  full_case_counts <- rep(NA, n_full)

  missing_idx <- !full_x %in% x
  non_missing_idx <- full_x %in% x
  full_case_counts[non_missing_idx] <- observed_counts

  zoo::na.fill(full_case_counts, fill = "extend")
}
