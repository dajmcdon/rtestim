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
    delay_distn = NULL) {

  arg_is_length(2, dist_gamma)
  arg_is_positive(dist_gamma)
  arg_is_positive(delay_distn, allow_null = TRUE)
  n <- length(observed_counts)
  arg_is_length(n, x, allow_null = TRUE)
  if (is.null(x)) x <- 1:n
  else {
    if (any(is.na(x)))
      cli::cli_abort("x may not contain missing values.")
    if (is.unsorted(x, strictly = TRUE))
      cli::cli_abort("x must be sorted and contain no duplicates.")
  }

  if (!is.null(delay_distn)) delay_distn <- delay_distn / sum(delay_distn)
  regular <- vctrs::vec_unique_count(diff(x)) == 1L
  if (regular) xout <- x
  else xout <- seq(from = min(x), to = max(x), by = min(diff(x)))

  if (is.null(delay_distn)) {
    delay_distn <- discretize_gamma(xout, dist_gamma[1], dist_gamma[2])
  } else {
    if (length(delay_distn) > length(xout)) {
      cli::cli_abort(
        "User specified `w` must have no more than {length(xout)} elements."
      )
    }
    # pad the tail with zero if too short
    delay_distn <- c(delay_distn, rep(0, length(xout) - length(delay_distn)))
  }

  y <- approx(x, observed_counts, xout = xout)$y
  cw <- cumsum(delay_distn)

  convolved_seq <- stats::convolve(
    y, rev(delay_distn), type = "open")[seq_along(xout)] / cw
  if (!regular) convolved_seq <- convolved_seq[xout %in% x]
  return(c(0, convolved_seq[1:(n - 1)]))
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

  output <- list()

  output$full_counts <- zoo::na.fill(full_case_counts, fill = "extend")
  output$missing_idx <- missing_idx

  return(output)
}
