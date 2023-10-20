new_period <- function(x) {
  arg_is_scalar(x)

  n <- as.numeric(strextract("^[0-9]+", x))
  if (length(n) == 0L) n <- 1L
  names_in <- strextract("[a-zA-Z]+$", x)
  if (length(names_in) == 0L) names_in <- "days"
  names_allowed <- paste0(rlang::fn_fmls_names(default_period), "s")
  if (is.na(pmatch(names_in, names_allowed))) {
    cli_abort(c(
      "Requested periodicity {.var {names_in}} is not available.",
      i = "Must be one of {.val {names_allowed}}."
    ))
  }
  names_in <- gsub("s$", "", names_in)
  l <- rlang::list2(!!names_in := n)
  res <- eval(rlang::call2("default_period", !!!l))
  vctrs::new_rcrd(res, class = "period")
}

default_period <- function(year = 0, quarter = 0, month = 0, week = 0, day = 0) {
  enlist(year = year, month = month + 3 * quarter, day = day + 7 * week)
}

#' @method format period
#' @export
format.period <- function(x, ...) {
  nms <- c("Y", "M", "D")
  val <- vctrs::vec_c(!!!vctrs::vec_data(x))
  paste0(val[val != 0], nms[val != 0])
}

strextract <- function(pattern, x) {
  m <- regexec(pattern, x)
  unlist(regmatches(x, m))
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

