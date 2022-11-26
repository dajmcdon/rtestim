# package level documentation

utils::globalVariables(c("R_rate", "Time", "pois_para", "Signal"))
Est <- Type <- NULL
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang abort is_integerish
#' @importFrom Matrix bandSparse diff
#' @importFrom RSpectra eigs
#' @importFrom spam crossprod.spam
#' @importFrom data.table data.table
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by
#' @importFrom ggplot2 ggplot geom_point geom_line labs scale_colour_manual aes theme_bw
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils globalVariables
#' @useDynLib rtestim, .registration = TRUE
NULL
