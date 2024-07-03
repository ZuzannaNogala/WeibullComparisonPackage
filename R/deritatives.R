globalVariables(c("Y_val"))
#' Deritatives used in Newton-Raphson Algorithm
#'
#' @param T_k numeric, k-th observation's value
#' @param delta vector with 0 and 1 values, information of censorship each observation
#' in sample
#' @param theta_k numeric, k-th estimatation of parameter theta (shape parameter)
#' @param lambda_k numeric, k-th estimatation of parameter lambda (scale parameter)
#'
#' @rdname deritatives
#' @keywords internal
first_deritative_theta_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("first_deritative_theta_k")
}

#' @keywords internal
first_deritative_theta_k.default <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  d_k * 1 / theta_k - d_k * log(lambda_k) +
    sum(delta * log(T_k)) +
    log(lambda_k) * lambda_k ^ (-theta_k) * sum(T_k ^ theta_k) -
    lambda_k ^ (-theta_k) * sum(log(T_k) * T_k ^ theta_k)
}

#' @rdname deritatives
#' @keywords internal
first_deritative_lambda_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("first_deritative_lambda_k")
}

#' @import data.table
#'
#' @keywords internal
first_deritative_lambda_k.default <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  - d_k * theta_k / lambda_k + theta_k * lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k)
}

#' @rdname deritatives
#' @keywords internal
second_deritative_theta_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("second_deritative_theta_k")
}

#' @import data.table
#'
#' @keywords internal
second_deritative_theta_k.default <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  - d_k / theta_k ^ 2 -
    log(lambda_k) ^ 2 * lambda_k ^ (-theta_k) * sum(T_k ^ theta_k) +
    2 * log(lambda_k) * lambda_k ^ (-theta_k) * sum(log(T_k) * T_k ^ theta_k) -
    lambda_k ^ (-theta_k) * sum(log(T_k) ^ 2 * T_k ^ theta_k)
}

#' @rdname deritatives
#' @keywords internal
second_deritative_lambda_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("second_deritative_lambda_k")
}

#' @import data.table
#'
#' @keywords internal
second_deritative_lambda_k.default <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  d_k * theta_k / lambda_k ^ 2 - theta_k * (theta_k + 1) * lambda_k ^ (-theta_k - 2) * sum(T_k ^ theta_k)
}

#' @rdname deritatives
#' @keywords internal
mix_deritative_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("mix_deritative_k")
}

#' @import data.table
#'
#' @keywords internal
mix_deritative_k.default <- function(T_k, delta, theta_k, lambda_k){
  d_k <- sum(delta)

  -d_k / lambda_k +
    lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k) -
    log(lambda_k) * theta_k * lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k) +
    theta_k * lambda_k ^ (-theta_k - 1) * sum(log(T_k) * T_k ^ theta_k)
}

#' @rdname deritatives
#' @keywords internal
grad_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("grad_k")
}

#' @import data.table
#'
#' @keywords internal
grad_k.default <- function(T_k, delta, theta_k, lambda_k){
  dx_theta_k <- first_deritative_theta_k(T_k, delta, theta_k, lambda_k)
  dx_lambda_k <- first_deritative_lambda_k(T_k, delta, theta_k, lambda_k)

  c(dx_theta_k, dx_lambda_k)
}

#' @rdname deritatives
#' @keywords internal
Hessian_k <- function(T_k, delta, theta_k, lambda_k) {
  UseMethod("Hessian_k")
}

#' @import data.table
#'
#' @keywords internal
Hessian_k.default <- function(T_k, delta, theta_k, lambda_k){
  H_k <- matrix(0, ncol = 2, nrow = 2)
  H_k[1,1] <- second_deritative_theta_k(T_k, delta, theta_k, lambda_k)
  H_k[2,2] <- second_deritative_lambda_k(T_k, delta, theta_k, lambda_k)
  H_k[1,2] <- mix_deritative_k(T_k, delta, theta_k, lambda_k)
  H_k[2,1] <- H_k[1,2]

  H_k
}
