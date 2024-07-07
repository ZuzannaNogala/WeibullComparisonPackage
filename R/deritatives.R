globalVariables(c("Y_val"))
#' Derivatives of loglikehood function of sample from Weibull distribution
#' (full and censored samples)
#'
#' The loglikehood function of sample from from Weibull distribution is described by from:
#' \deqn{l(\lambda, \theta|\textbf{X}) = \log{(\theta)}\sum^n_{i=1}\delta_i
#'  - \theta\log{(\lambda)}\sum^n_{i=1}\delta_i + (\theta - 1)\sum^n_{i=1}
#'  \delta_i\log{(X_i)} - \lambda^{-\theta}\sum^n_{i=1}X_i^{\theta}.}
#' If the sample is full \eqn{\sum^n_{i=1}\delta_i = n}. The functions below compute values of
#' first-order and second-order partial derivatives of \eqn{l(\lambda, \theta|\textbf{X})} with respect
#' to \eqn{\theta}, \eqn{\lambda} at \code{T_k}. Additionally, they compute mixed partial derivative,
#' Hessian and gradient of \eqn{l(\lambda, \theta|\textbf{X})}.
#'
#' @param T_k numeric, a value of k-th observation, a point where the function are computed
#' @param delta vector with 0 and 1 values, information of censorship all observations
#' in sample
#' @param theta_k numeric, k-th estimation of parameter \eqn{\theta} (shape parameter)
#' @param lambda_k numeric, k-th estimation of parameter \eqn{\lambda} (scale parameter)
#'
#' @rdname derivatives
#' @keywords internal
first_deritative_theta_k<- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  d_k * 1 / theta_k - d_k * log(lambda_k) +
    sum(delta * log(T_k)) +
    log(lambda_k) * lambda_k ^ (-theta_k) * sum(T_k ^ theta_k) -
    lambda_k ^ (-theta_k) * sum(log(T_k) * T_k ^ theta_k)
}

#' @rdname derivatives
#' @keywords internal
first_deritative_lambda_k <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  - d_k * theta_k / lambda_k + theta_k * lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k)
}

#' @rdname derivatives
#' @keywords internal
second_deritative_theta_k <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  - d_k / theta_k ^ 2 -
    log(lambda_k) ^ 2 * lambda_k ^ (-theta_k) * sum(T_k ^ theta_k) +
    2 * log(lambda_k) * lambda_k ^ (-theta_k) * sum(log(T_k) * T_k ^ theta_k) -
    lambda_k ^ (-theta_k) * sum(log(T_k) ^ 2 * T_k ^ theta_k)
}

#' @rdname derivatives
#' @keywords internal
second_deritative_lambda_k <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  d_k * theta_k / lambda_k ^ 2 - theta_k * (theta_k + 1) * lambda_k ^ (-theta_k - 2) * sum(T_k ^ theta_k)
}

#' @rdname derivatives
#' @keywords internal
mix_deritative_k<- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  -d_k / lambda_k +
    lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k) -
    log(lambda_k) * theta_k * lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k) +
    theta_k * lambda_k ^ (-theta_k - 1) * sum(log(T_k) * T_k ^ theta_k)
}

#' @rdname derivatives
#' @keywords internal
grad_k <- function(T_k, delta, theta_k, lambda_k){
  dx_theta_k <- first_deritative_theta_k(T_k, delta, theta_k, lambda_k)
  dx_lambda_k <- first_deritative_lambda_k(T_k, delta, theta_k, lambda_k)

  c(dx_theta_k, dx_lambda_k)
}

#' @rdname derivatives
#' @keywords internal
Hessian_k<- function(T_k, delta, theta_k, lambda_k){
  H_k <- matrix(0, ncol = 2, nrow = 2)
  H_k[1,1] <- second_deritative_theta_k(T_k, delta, theta_k, lambda_k)
  H_k[2,2] <- second_deritative_lambda_k(T_k, delta, theta_k, lambda_k)
  H_k[1,2] <- mix_deritative_k(T_k, delta, theta_k, lambda_k)
  H_k[2,1] <- H_k[1,2]

  H_k
}
