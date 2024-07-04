globalVariables(c("T_val", "Y_val"))
#' Transformed Data to variable from extreme value distribution
#'
#' @param T_vec_dt a data table which contains value of T statistic and information of observation's censoring
#' @param theta_hat a estimator of parameter theta (shape parameter)
#' @param lambda_hat a estimator of parameter lambda (scale parameter)
#'
#' @import data.table
#' @keywords internal
transform_to_Y <- function(T_vec_dt, theta_hat, lambda_hat){
  T_vec_dt[, .(Y_val = theta_hat * (log(T_val) - log(lambda_hat)) , delta)]
}

#' Sorting transformed data
#'
#' @param Y_vec_dt a data table which contains value of transformed data statistic Y and
#' information of observation's censoring
#'
#' @import data.table
#' @keywords internal
sort_Y_vals <- function(Y_vec_dt){
  Y_values <- Y_vec_dt[, Y_val]
  names(Y_values) <- Y_vec_dt[, delta]
  sort_Y_values <- sort(Y_values)
  if(sum(Y_vec_dt[, delta]) != length(Y_values)){
    sort_delta <- as.logical(names(sort_Y_values))
  }
  else{
    sort_delta <- rep(TRUE, length(Y_values))
  }
  dt <- data.table(sort_Y_val = sort_Y_values, delta = sort_delta)
  dt
}

#' Kaplan-Meier Estimator of distribution of sample
#'
#' @param t numeric, the value, where the estimator will be counted
#' @param n integer, number of observations in data
#' @param sort_vec_dt a data table which contains value of sorted statistic and
#' information of observation's censoring
#'
#' @import data.table
#' @keywords internal
Est_kaplan_meier <- function(t, n, sort_vec_dt){
  Vec <- sort_vec_dt[, get(names(sort_vec_dt)[1])]
  delta_Vec <- sort_vec_dt[, get(names(sort_vec_dt)[2])]
  if(Vec[1] > t) return(0)
  if(Vec[n] < t) return(1)
  else{
    k <- which.max(Vec >= t)
    k_before <- k
    expr1 <- prod(((n - 1:k_before) / (n + 1 - 1:k_before)) ^ delta_Vec[1:k_before])
    return(1 - expr1)
  }
}


#' Function
#'
#' @param x a numeric vector with observations of full sample
#' @param perc_of_cens numeric, percentage of censorship (value from (0, 1])
#'
#' @importFrom stats punif
#' @keywords internal
m_find_function <- function(x, perc_of_cens){
  function(m){
    mean(punif(x, min = 0, max = m)) - perc_of_cens
  }
}

#' Function to find
#'
#' @param x a numeric vector with observations of full sample
#' @param perc_of_cens numeric, percentage of censorship (value from (0, 1])
#'
#' @export
m_find <- function(x, perc_of_cens){
  UseMethod("m_find")
}

#' @importFrom pracma bisect
#' @export
m_find.default <- function(x, perc_of_cens){
  pracma::bisect(m_find_function(x, perc_of_cens), 0, 10000)$root
}
