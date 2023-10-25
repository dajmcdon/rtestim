new_period <- function(x) {
  arg_is_scalar(x)
  arg_is_chr(x)

  n <- as.integer(strextract("^[0-9]+", x))
  names_in <- tolower(strextract("[a-zA-Z]+$", x))
  names_allowed <- paste0(rlang::fn_fmls_names(default_period), "s")
  if (length(n) == 0L || length(names_in) == 0L ||
      is.na(pmatch(names_in, names_allowed))) {
    cli_abort(c(
      "Requested periodicity {.var {names_in}} is not available.",
      i = "Input must be a positive integer followed by one of {.val {names_allowed}}."
    ))
  }
  names_in <- gsub("s$", "", names_in)
  l <- rlang::list2(!!names_in := n)
  res <- eval(rlang::call2("default_period", !!!l))
  vctrs::new_rcrd(res, class = "period")
}

default_period <- function(year = 0L, quarter = 0L, month = 0L, week = 0L,
                           day = 0L) {
  enlist(year = year, month = month + 3L * quarter, day = day + 7L * week)
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
  if (anyNA(x)) {
    return(NA)
  }
  stopifnot(is.numeric(x))
  if (!rlang::is_integerish(x)) cli_abort("`x` must contain only integers.")
  if (length(x) == 1L) {
    return(as.integer(x))
  }
  x <- x[x != 0]
  compute_gcd(x)
}
