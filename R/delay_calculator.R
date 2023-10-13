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
#' @param xout a vector of positions at for which the results should be returned.
#'   By default, this will be the same as `x`, but in the case that observations
#'   are unequally spaced, alternatives may be desired. Note that `xout` must
#'   satisfy `min(x) <= min(xout)` and `max(x) >= max(xout)`.
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
    xout = x) {

  arg_is_length(2, dist_gamma)
  arg_is_positive(dist_gamma)
  arg_is_positive(delay_distn, allow_null = TRUE)
  n <- length(observed_counts)
  arg_is_length(n, x, allow_null = TRUE)
  if (is.null(x)) x <- 1:n
  else {
    if (any(is.na(x))) cli_abort("`x` may not contain missing values.")
    if (is.unsorted(x, strictly = TRUE))
      cli_abort("`x` must be sorted and contain no duplicates.")
  }

  if (inherits(x, "Date")) x <- as.numeric(x)
  arg_is_numeric(x)
  if (!is.null(delay_distn)) delay_distn <- delay_distn / sum(delay_distn)

  if (any(is.na(xout))) cli_abort("`xout` may not contain missing values.")
  if (is.unsorted(xout, strictly = TRUE))
    cli_abort("`xout` must be sorted and contain no duplicates.")
  if (inherits(xout, "Date")) xout <- as.numeric(xout)
  arg_is_numeric(xout)
  if (min(xout) < min(x)) cli_abort("`min(xout)` man not be less than `min(x)`.")
  if (max(xout) > max(x)) cli_abort("`max(xout)` man not exceed `max(x)`.")

  allx <- union(x, xout)
  dallx <- diff(allx)

  ## TODO: handle weekly / monthly incidence automatically
  # regular <- vctrs::vec_unique_count(dallx) == 1L
  # if (!regular) {
  min_spacing <- 1L  #gcd(unique(dallx))
  allx <- seq(from = min(x), to = max(x), by = min_spacing)
  # }

  if (is.null(delay_distn)) {
    delay_distn <- discretize_gamma(allx, dist_gamma[1], dist_gamma[2])
  } else {
    if (length(delay_distn) > length(allx)) {
      cli_abort(
        "User specified `delay_distn` must have no more than {length(allx)} elements."
      )
    }
    # pad the tail with zero if too short
    delay_distn <- c(delay_distn, rep(0, length(allx) - length(delay_distn)))
  }

  y <- stats::approx(x, observed_counts, xout = allx)$y
  cw <- cumsum(delay_distn)
  convolved_seq <- stats::convolve(y, rev(delay_distn), type = "open")
  convolved_seq <- convolved_seq[seq_along(allx)] / cw
  convolved_seq <- c(convolved_seq[1], convolved_seq[-length(convolved_seq)])
  convolved_seq[allx %in% xout]
}


gcd <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (anyNA(x)) return(NA)
  stopifnot(is.numeric(x))
  if (length(x) < 2L) return(x)
  if (!rlang::is_integerish(x)) cli_abort("`x` must contain only integers.")
  x <- x[x != 0]
  compute_gcd(x)
}
