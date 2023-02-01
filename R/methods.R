#' Summarize results of a poisson_rt model
#'
#' @method summary poisson_rt
#'
#' @param object a fitted model of class `poisson_rt`
#' @param ... .
#'
#' @return a data table of estimates and a logical value of convergence.
#' The data table columns include the current observed daily counts (Signal),
#' the estimated reproduction rate (R_rate), and the estimated Poisson mean
#' parameter (pois_mean)
#'
#' @export
#'
#' @examples
#' TODO: Need to change example
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod = admm_solver(
#'   current_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(current_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
#' )
#' summary(mod)
summary.poisson_rt <- function(object, ...){
  n = length(object$observed_counts)
  res <- data.frame(
    Time = 1:n,
    Signal = object$observed_counts,
    R_rate = object$Rt,
    pois_mean = object$Rt * object$weighted_past_counts
  )

  lst = list(Results = res)
  class(lst) = "summary.poisson_rt"
  return(lst)
}


#' Plot predicted observed_count and estimated Rt from `summary(poisson_rt)` object
#'
#' @method plot summary.poisson_rt
#' @param x summary of `poisson_rt` models
#' @param ... .
#'
#' @return Figure with two panels. Top panel shows the predicted observed_counts
#' calculated from \eqn{weighted_count * Rt_estim}
#' @export
#'
#' @examples
#' TODO: change this example
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- estimate_rt(
#'   observed_counts = y, degree = 1, lambda = 1)
#' plot(summary(mod))
plot.summary.poisson_rt <- function(summary, ...){
  fig_cases <- summary$Results %>%
    ggplot(aes(x = .data$Time)) +
    geom_point(aes(y = .data$Signal)) +
    geom_line(aes(y = .data$pois_mean), col = "#08519C") +
    labs(x = "Time", y = "Daily infection counts (on dots)",
         title = "The estimated piecewise polynomial curve (in line)") +
    theme_bw()

  fig_rt <- summary$Results %>%
    ggplot(aes(x = .data$Time)) +
    geom_line(aes(y = .data$R_rate), col = "#14754C") +
    labs(x = "Time", y = "Estimated Rt",
         title = "The estimated Rt") +
    theme_bw()

  fig <- ggpubr::ggarrange(fig_cases, fig_rt, ncol=2)

  print(fig)
}

#' Summary of `cv_result` object
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.cv_result <- function(object, ...) {
  res <- list(
    lambdas = object$lambda,
    cv_scores = object$cv_scores,
    optimal_Rt = object$optimal_Rt,
    optimal_lambda = object$optimal_lambda,
    x = object$x,
    weighted_past_counts = object$weighted_past_counts,
    observed_counts = object$observed_counts,
    pois_mean = object$optimal_Rt*object$weighted_past_counts
  )

  class(res) = "summary.cv_result"

  return(res)
}

#' Plot of `summary.cv_result` object
#'
#' @param cv_result
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.summary.cv_result <- function(cv_result, ...) {
  lambda <- cv_result$lambdas
  cv_scores <- cv_result$cv_scores
  Rt <- cv_result$optimal_Rt
  x <- cv_result$x
  weighted_past_counts <- cv_result$weighted_past_counts
  observed_counts <- cv_result$observed_counts
  pois_mean <- cv_result$pois_mean

  # Score plot
  cv_scores_plot <- ggplot(data.frame(lambdas =lambda, scores = cv_scores),
                           aes(x=lambdas))+
    geom_line(aes(x = lambda, y = scores),col = "#08519C")+
    labs(x = "Lambdas", y = "Cross Validation Scores",
         title = "Cross validation scores") +
    theme_bw()

  # Rt plot
  optimal_rt_plot <- ggplot(data.frame(x=x, Rt = Rt), aes(x=x, y = Rt))+
    geom_line(col = "#08519C")+
    labs(x = "Observation time point", y = "Rt",
         title = "Optimal Rt from Cross Validation") +
    geom_hline(yintercept=1, linetype="dashed", color = "red")+
    theme_bw()

  # True vs pred
  truevspred<- ggplot(data.frame(x = x, true = observed_counts, pred=pois_mean),
                      aes(x = x)) +
    geom_point(aes(y = observed_counts), shape = 1) +
    geom_line(aes(y = pois_mean), col = "#14754C") +
    labs(x = "Time", y = "Daily infection counts (on dots)",
         title = "The estimated piecewise polynomial curve (in line)") +
    theme_bw()

  fig <- ggpubr::ggarrange(cv_scores_plot,
                           optimal_rt_plot,
                           truevspred,
                           ncol=3)

  print(fig)
}

