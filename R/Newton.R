globalVariables(c("T_vals"))
#' Newton-Raphson Algorithm - k-th step of algorithm
#'
#' @param theta_before_k (k-1)-th estimatation of parameter theta (shape parameter)
#' @param lambda_before_k (k-1)-th estimatation of parameter theta (scale parameter)
#' @param Hess_before_k ..
#' @param grad_before_k ..
#'
#' @keywords internal
Newton_k_step <- function(theta_before_k, lambda_before_k, Hess_before_k, grad_before_k) {
  UseMethod("Newton_k_step")
}

#'
#' @keywords internal
Newton_k_step.numeric <- function(theta_before_k, lambda_before_k, Hess_before_k, grad_before_k){
  x_before_k <- c(theta_before_k, lambda_before_k)
  x_before_k - as.vector(solve(Hess_before_k) %*% grad_before_k)
}

#' Starting Points - theta
#'
#' @param n integer, number of observations in data
#' @param T_vals vector of observations
#'
#' @keywords internal
theta_hat_approx <- function(n, T_vals) {
  UseMethod("theta_hat_approx")
}

#'
#' @keywords internal
theta_hat_approx.numeric <- function(n, T_vals){
  function(x) n / x + sum(log(T_vals)) - (n / sum(T_vals ^ x)) * sum(log(T_vals) * (T_vals) ^ x)
}

#' Starting Points - theta, full sample
#'
#' @param n integer, number of observations in data
#' @param T_vals vector of observations
#'
#' @keywords internal
mle_theta_full_sample <- function(n, T_vals) {
  UseMethod("mle_theta_full_sample")
}

#' @importFrom pracma bisect
#' @keywords internal
mle_theta_full_sample.numeric <- function(n, T_vals) {
  pracma::bisect(theta_hat_approx(n, T_vals), 0, 10)$root
  }


#' Starting Points - lambda, full sample
#'
#' @param n integer, number of observations in data
#' @param T_vals vector of observations
#' @param theta_hat_est estimation of parametr theta (shape parametr); result of
#' function \code{\link{mle_theta_full_sample}}
#'
#' @keywords internal
mle_lambda_full_sample <- function(n, T_vals, theta_hat_est = mle_theta_full_sample(n, T_vals)){
  UseMethod("mle_lambda_full_sample")
}

#' @keywords internal
mle_lambda_full_sample.numeric <- function(n, T_vals, theta_hat_est = mle_theta_full_sample(n, T_vals)){
  (1 / n * (sum(T_vals ^ theta_hat_est))) ^ (1 / theta_hat_est)
}

#' Newton-Raphson Algorithm
#'
#' @param T_dt data table with one column with values of observation, second one contains information if the observation
#' is censured or not
#'
#' @export
Newton_result <- function(T_dt){
  UseMethod("Newton_result")
}

#' @import data.table
#' @keywords internal
Newton_result.default <- function(T_dt){
  theta_0 <- mle_theta_full_sample(T_dt[, .N], T_dt[, get(names(T_dt)[1])])
  lambda_0 <- mle_lambda_full_sample(T_dt[, .N], T_dt[, get(names(T_dt)[1])], theta_0)
  T_v <- T_dt[, get(names(T_dt)[1])]
  delta <-  T_dt[, get(names(T_dt)[2])]
  x_k <- c(theta_0, lambda_0)
  Hess_k <- Hessian_k(T_v, delta, x_k[1], x_k[2])
  grd_k <- grad_k(T_v, delta, x_k[1], x_k[2])
  x_k2 <- Newton_k_step(x_k[1], x_k[2], Hess_k, grd_k)

  while((abs(x_k2[1] - x_k[1]) > 0.00005 | abs(x_k2[2] - x_k[2]) > 0.00005)){
    x_k <- c(x_k2[1], x_k2[2])
    Hess_k <- Hessian_k(T_v, delta, x_k[1], x_k[2])
    if(NaN %in% Hess_k) break
    grd_k <- grad_k(T_v, delta, x_k[1], x_k[2])
    if(NaN %in% grd_k) break
    x_k2 <- try(Newton_k_step(x_k[1], x_k[2], Hess_k, grd_k), silent = TRUE)
    if(isTRUE(is.character(x_k2[1]))) break
  }

  if(isTRUE(is.character(x_k2[1])) | isTRUE(is.character(x_k[1])) ){
    theta_k2 <- theta_0
    lambda_k2 <- lambda_0
  }
  else{
    theta_k2 <- x_k2[1]
    lambda_k2 <- x_k2[2]
  }

  c(theta_k2, lambda_k2)
}






