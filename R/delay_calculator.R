#' Calculate the total infectiousness at each observed time point.
#'
#' The total infectiousness at each observed time point is calculated
#' by \eqn{\sum_{s=1}^t I_{t-s}w_s}, where \eqn{I} denotes the vector containing
#' observed incidence, and \eqn{w} denotes the generation interval
#' distribution. Typically, the generation interval is challenging to estimate
#' from data, so the serial interval is used instead. The serial interval
#' distribution expresses the probability
#' of a secondary infection caused by a primary infection which occurred \eqn{s}
#' days earlier.
#'
#' @inheritParams estimate_rt
#' @param xout a vector of positions at which the results should be returned.
#'   By default, this will be the same as `x`, but in the case that observations
#'   are unequally spaced, alternatives may be desired. Note that `xout` must
#'   satisfy `min(x) <= min(xout)` and `max(x) >= max(xout)`.
#'
#'
#' @return A vector containing the total infectiousness at each
#'  point `xout`.
#' @export
#'
#' @examples
#' delay_calculator(c(3, 2, 5, 3, 1), dist_gamma = c(2.5, 2.5))
delay_calculator <- function(
    observed_counts,
    x = NULL,
    dist_gamma = c(2.5, 2.5),
    delay_distn = NULL,
    delay_distn_periodicity = NULL,
    xout = x) {
  assert_numeric(dist_gamma, lower = 0, len = 2)
  if (inherits(delay_distn, "Matrix")) {
    assert_numeric(delay_distn@x, lower = 0)
  } else {
    assert_numeric(delay_distn, lower = 0, null.ok = TRUE)
  }

  n <- length(observed_counts)
  if (is.null(x)) {
    x <- 1:n
  } else {
    if (inherits(x, "Date")) x <- as.numeric(x)
    assert_numeric(x, len = n, any.missing = FALSE)
    if (is.unsorted(x, strictly = TRUE)) {
      cli_abort("`x` must be strictly sorted and contain no duplicates.")
    }
  }

  if (is.null(xout)) {
    xout <- x
  } else {
    if (inherits(xout, "Date")) xout <- as.numeric(xout)
    assert_numeric(xout, any.missing = FALSE)
    if (is.unsorted(xout, strictly = TRUE)) {
      cli_abort("`xout` must be strictly sorted and contain no duplicates.")
    }
  }
  if (min(xout) < min(x)) cli_abort("`min(xout)` may not be less than `min(x)`.")
  if (max(xout) > max(x)) cli_abort("`max(xout)` may not exceed `max(x)`.")

  ddx <- gcd(unique(diff(x)))
  if (is.null(delay_distn_periodicity)) {
    ddp <- ddx
  } else if (is.character(delay_distn_periodicity)) {
    ddp <- new_period(delay_distn_periodicity)
    ddp <- vctrs::field(ddp, "day")
  } else if (is.numeric(delay_distn_periodicity)) {
    ddp <- delay_distn_periodicity
  } else {
    cli_abort("`delay_distn_periodicity` must be a scalar, character or numeric.")
  }
  if (ddx %% ddp != 0) {
    cli_abort(c(
      "`delay_distn_periodicity` may be at most the minimum spacing in `x`,",
      "!" = "and must divide the minimum spacing evenly.",
      i = "`delay_distn_periodicity` = {.val {ddp}} compared to {.val {ddx}} for `x`."
    ))
  }
  allx <- seq(from = min(x), to = max(x), by = ddp)
  ddxout <- gcd(unique(diff(xout)))
  if (ddxout < ddp) {
    cli_abort(c(
      "The minimum spacing in `xout` must be at least `delay_distn_periodicity`.",
      i = "`delay_distn_periodicity` = {.val {ddp}} compared to {.val {ddxout}} for `xout`."
    ))
  }

  y <- stats::approx(x, observed_counts, xout = allx)$y

  if (is.null(delay_distn)) {
    delay_distn <- discretize_gamma(allx - min(x), dist_gamma[1], dist_gamma[2])
  } else if (is.matrix(delay_distn) || inherits(delay_distn, "Matrix")) {
    if (!all(dim(delay_distn) == length(allx))) {
      cli_abort(c(
        "User specified `delay_distn` has dimensions {.val {dim(delay_distn)}},",
        "!" = "but it must have both dimensions {.val {length(allx)}}."
      ))
    }
    if (!Matrix::isTriangular(delay_distn, upper = FALSE)) {
      cli_abort(
        "User specified `delay_distn` must be square and lower triangular."
      )
    }
    if (delay_distn[1,1] < .Machine$double.eps) {
      cli_abort(c(
        "The upper-left most entry of `delay_distn` must be strictly positive.",
        i = "Otherwise, we will eventually try to calculate 0/0 and produce {.val NaN}. "
      ))
    }
    delay_distn <- delay_distn / Matrix::rowSums(delay_distn)
    convolved_seq <- drop(delay_distn %*% y)
    return(convolved_seq[allx %in% xout])
  } else {
    if (length(delay_distn) > length(allx)) {
      cli_warn(c(
        "User specified `delay_distn` has {.val {length(delay_distn)}} when",
        "only {.val {length(allx)}} are necessary.",
        i = "Truncating to match."
      ))
      delay_distn <- delay_distn[seq_along(allx)]
    } else {
      delay_distn <- c(delay_distn, rep(0, length(allx) - length(delay_distn)))
    }
    delay_distn <- delay_distn / sum(delay_distn)
  }

  convolved_seq <- fast_convolve(y, delay_distn)
  # (polish up the beginning of the delay calculation)
  # when delay_distn[1] == 0, we're putting no weight on today.
  # This is typical and results in division by zero for the first observation
  # in convolved_seq
  if (abs(delay_distn[1]) < sqrt(.Machine$double.eps)) {
    convolved_seq <- c(convolved_seq[2], convolved_seq[-1])
  }
  convolved_seq[allx %in% xout]
}
