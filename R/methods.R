#' Summarize results of a poisson_rt model
#'
#' @method summary poisson_rt
#'
#' @param object output of `estimate_rt` function with class `poisson_rt`
#' @param ... .
#'
#' @return a data table of estimates and a logical value of convergence.
#' The data table columns include the current observed daily counts (Signal),
#' the estimated reproduction rate (R_rate), and the estimated Poisson mean
#' parameter (pois_mean)
#'
#' @exportS3Method
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod = estiamte_rt(observed_counts = y, degree = 1, lambda = 0.001)
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
#' @param summary summary of `poisson_rt` models, of class `summary.poisson_rt`
#' @param ... .
#'
#' @return Panel with two figures. Left figure shows the predicted observed_counts
#' calculated from \eqn{weighted_past_count * R} against the true observed_counts.
#' Right figure shows the estimated Rt across observation time
#' @exportS3Method
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- estimate_rt(observed_counts = y, degree = 1, lambda = 0.001)
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
    ylim(0, 10)+
    geom_hline(yintercept=1, linetype="dashed", color = "red")+
    labs(x = "Time", y = "Estimated Rt",
         title = "The estimated Rt") +
    theme_bw()

  fig <- ggpubr::ggarrange(fig_cases, fig_rt, ncol=2)

  print(fig)
}

#' Summary of `cv_result` object
#'
#' @param object output of `cv_estimate_rt` function with class `cv_result`
#' @param ... .
#'
#' @return a data table containing `weighted_past_counts`, `observed_counts`,
#' `lambda` and `x` which are used to find the optimal lambda from cross
#' validation; `cv_scores`, `optimal_Rt`, `optimal_lambda`, `pois_mean` are the
#' results of the cross validation
#' @exportS3Method
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- cv_estimate_rt(observed_counts = y, degree = 1,
#' lambda = log(seq(1.01, 5, 0.05)))
#' summary(mod)
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
#' @param cv_result summary of the function `cv_estimate_rt` with class
#' `summary.cv_result`
#' @param ...
#'
#' @return Generates three figures. First figure shows the cross validation scores
#' of each lambda. Second figure shows the optimal Rt, estimated with
#' the optimal lambda from cross validation, across observed time point. Third
#' figure shows the predicted observed case counts calculated by
#' \eqn{weighted_past_count * optimal_Rt}, against the true observed_counts
#' @exportS3Method
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- cv_estimate_rt(observed_counts = y, degree = 1,
#' lambda = log(seq(1.01, 5, 0.05)))
#' plot(summary(mod))
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

  fig <- ggpubr::ggarrange(optimal_rt_plot,
                           cv_scores_plot,
                           truevspred,
                           ncol=2, nrow=2)

  print(fig)
}

