#' Statistics of new tests
#'
#' @param T_vec_dt a data table which contains value of T statistic and information of observation's censoring
#' @param theta_hat a estimator of parameter theta (shape parameter)
#' @param lambda_hat a estimator of parameter lambda (scale parameter)
#'
#' @rdname statisticTests
#' @keywords internal
TransformedDataToY <- function(T_vec_dt, theta_hat, lambda_hat) {
  UseMethod("TransformedDataToY")
}

#' @import data.table
#'
#' @keywords internal
TransformedDataToY.default <- function(T_vec_dt, theta_hat, lambda_hat){
  T_vec_dt[, .(Y_val = theta_hat * (log(T_val) - log(lambda_hat)) , delta)]
}
