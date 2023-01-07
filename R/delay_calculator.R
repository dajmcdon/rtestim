#' Delay Calculator
#'
#' @param current_counts the current daily infection counts
#' @param dist_gamma scale and shape parameter of the discretized Gamma distribution
#'
#' @return list of counts weighted by the discretized Gamma distribution
#' @export
#'
#' @examples delay_calculator(c(3,2,5,3,1), c(2,2))
delay_calculator <- function(current_counts, dist_gamma){
  n = length(current_counts)
  output = rep(0, n)
  output[1] = current_counts[1]

  for(idx in 2:n){
    w  = discretize_gamma(idx, dist_gamma[1], dist_gamma[2])
    output[idx] = current_counts[idx:1]*w
  }
  return(output)
}
