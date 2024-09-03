globalVariables(c("T_vals"))
#' Newton-Raphson Algorithm - k-th step of algorithm
#'
#' @param theta_before_k (k-1)-th estimatation of parameter theta (shape parameter)
#' @param lambda_before_k (k-1)-th estimatation of parameter theta (scale parameter)
#' @param Hess_before_k a value of Hessian of function \eqn{l(\theta, \lambda)} (more in \code{\link{derivatives}}) in point (\code{theta_before_k}, \code{lambda_before_k})
#' @param grad_before_k a value of gradient  of function \eqn{l(\theta, \lambda)} in point (\code{theta_before_k}, \code{lambda_before_k})
#'
#' @keywords internal
Newton_k_step <- function(theta_before_k, lambda_before_k, Hess_before_k, grad_before_k){
  x_before_k <- c(theta_before_k, lambda_before_k)
  x_before_k - as.vector(solve(Hess_before_k) %*% grad_before_k)
}

#' Starting points of Newton-Raphson Algorithm (censured case) - parameter \eqn{\theta}
#'
#' @param n integer, number of observations in data
#' @param T_vals vector of observations
#'
#' @keywords internal
theta_hat_approx <- function(n, T_vals){
  function(x) n / x + sum(log(T_vals)) - (n / sum(T_vals ^ x)) * sum(log(T_vals) * (T_vals) ^ x)
}

#' Starting points of Newton-Raphson Algorithm (full sample case) - parameter \eqn{\theta}
#'
#' @param n integer, number of observations in data
#' @param T_vals vector of observations
#'
#' @importFrom pracma bisect
#' @keywords internal
mle_theta_full_sample<- function(n, T_vals) {
  pracma::bisect(theta_hat_approx(n, T_vals), 0, 10)$root
}


#' Starting points of Newton-Raphson Algorithm (full sample case) - parameter \eqn{\lambda}
#'
#' @param n integer, number of observations in data
#' @param T_vals vector of observations
#' @param theta_hat_est estimation of parameter \eqn{\theta} (shape parametr); result of
#' function \code{\link{mle_theta_full_sample}}
#'
#' @keywords internal
mle_lambda_full_sample <- function(n, T_vals, theta_hat_est = mle_theta_full_sample(n, T_vals)){
  (1 / n * (sum(T_vals ^ theta_hat_est))) ^ (1 / theta_hat_est)
}

#' Maximum Likelihood Estimation of two-parameter Weibull distribution using Newton-Raphson algorithm
#'
#' Newton-Raphson algorithm is numerical method applied to find estimation of roots of function. In this case method was used
#' to find the root of loglikelihood function for Weibull distribution which is described by the formula:
#' \deqn{l(\lambda, \theta|\textbf{X}) =
#' \log{(\theta)}\sum^n_{i=1}\delta_i -
#' \theta\log{(\lambda)}\sum^n_{i=1}\delta_i +
#' (\theta - 1)\sum^n_{i=1}\delta_i\log{(X_i)} -
#' \lambda^{-\theta}\sum^n_{i=1}X_i^{\theta},}
#' where \eqn{\delta_i} ...
#' The starting points of algorithm are solution of ...
#'
#' @param T_dt data table with two columns. First one contains values of observation and second one - information if the observation
#' is censured or not
#'
#' @return a list which contains estimation of \eqn{\theta} and \eqn{\lambda} parameters
#'
#' @import data.table
#' @importFrom stats rweibull runif
#'
#' @examples
#' X <- rweibull(100, 2, 1)
#' C <- runif(100, 0, m_find(X, 0.1))
#' T_dt <- data.table("T_val" = pmin(X, C), "delta" = pmin(X, C) == X)
#' Newton_result(T_dt)
#'
#' @export
Newton_result <- function(T_dt){
  UseMethod("Newton_result")
}
#' @export
Newton_result.data.table <- function(T_dt){
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

  list("theta" = theta_k2, "lambda" = lambda_k2)
}






