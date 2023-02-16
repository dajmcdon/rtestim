# See https://github.com/r-lib/rlang/blob/main/R/compat-purrr.R


map <- function(.x, .f, ...) {
  .f <- rlang::as_function(.f, env = rlang::global_env())
  lapply(.x, .f, ...)
}

walk <- function(.x, .f, ...) {
  map(.x, .f, ...)
  invisible(.x)
}

walk2 <- function(.x, .y, .f, ...) {
  map2(.x, .y, .f, ...)
  invisible(.x)
}

map_lgl <- function(.x, .f, ...) {
  .rlang_purrr_map_mold(.x, .f, logical(1), ...)
}

map_int <- function(.x, .f, ...) {
  .rlang_purrr_map_mold(.x, .f, integer(1), ...)
}

map_dbl <- function(.x, .f, ...) {
  .rlang_purrr_map_mold(.x, .f, double(1), ...)
}

map_chr <- function(.x, .f, ...) {
  .rlang_purrr_map_mold(.x, .f, character(1), ...)
}

.check_list_of_data_frames <- function (x, error_call = rlang::caller_env()) {
  vctrs::vec_check_list(x, call = error_call)
  is_df_or_null <- map_lgl(x, function(x) is.data.frame(x) || is.null(x))
  if (all(is_df_or_null)) return()
  bad <- which(!is_df_or_null)
  cli::cli_abort(
    c("Each element of {.arg x} must be either a data frame or {.code NULL}.",
      i = "Elements {bad} are not."), arg = "x", call = error_call
  )
}

map_dfr <- function(.x, .f, ..., .id = NULL) {
  .f <- rlang::as_function(.f, env = rlang::global_env())
  res <- map(.x, .f, ...)
  .check_list_of_data_frames(res)
  vctrs::vec_rbind(!!!res, .names_to = .id)
}

map2_dfr <- function(.x, .y, .f, ..., .id = NULL) {
  .f <- rlang::as_function(.f, env = rlang::global_env())
  res <- map2(.x, .y, .f, ...)
  .check_list_of_data_frames(res)
  vctrs::vec_rbind(!!!res, .names_to = .id)
}

.rlang_purrr_map_mold <- function(.x, .f, .mold, ...) {
  .f <- rlang::as_function(.f, env = rlang::global_env())
  out <- vapply(.x, .f, .mold, ..., USE.NAMES = FALSE)
  names(out) <- names(.x)
  out
}

.rlang_purrr_args_recycle <- function(args) {
  lengths <- map_int(args, length)
  n <- max(lengths)

  stopifnot(all(lengths == 1L | lengths == n))
  to_recycle <- lengths == 1L
  args[to_recycle] <- map(args[to_recycle], function(x) rep.int(x, n))

  args
}

map2 <- function(.x, .y, .f, ...) {
  .f <- rlang::as_function(.f, env = rlang::global_env())
  out <- mapply(.f, .x, .y, MoreArgs = list(...), SIMPLIFY = FALSE)
  if (length(out) == length(.x)) {
    rlang::set_names(out, names(.x))
  } else {
    rlang::set_names(out, NULL)
  }
}
map2_lgl <- function(.x, .y, .f, ...) {
  as.vector(map2(.x, .y, .f, ...), "logical")
}
map2_int <- function(.x, .y, .f, ...) {
  as.vector(map2(.x, .y, .f, ...), "integer")
}
map2_dbl <- function(.x, .y, .f, ...) {
  as.vector(map2(.x, .y, .f, ...), "double")
}
map2_chr <- function(.x, .y, .f, ...) {
  as.vector(map2(.x, .y, .f, ...), "character")
}
imap <- function(.x, .f, ...) {
  map2(.x, names(.x) %||% seq_along(.x), .f, ...)
}

pmap <- function(.l, .f, ...) {
  .f <- as.function(.f)
  args <- .rlang_purrr_args_recycle(.l)
  do.call("mapply", c(
    FUN = list(quote(.f)),
    args, MoreArgs = quote(list(...)),
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  ))
}

reduce <- function(.x, .f, ..., .init) {
  f <- function(x, y) .f(x, y, ...)
  Reduce(f, .x, init = .init)
}
