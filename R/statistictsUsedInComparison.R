globalVariables(c("sort_Y_val", "delta", "."))
#' Komogorov-Smirnov test - Full Sample Case
#'
#' @param X_vec numeric, the values from uncensured sample
#' @param theta_hat a estimator of parameter theta (shape parameter)
#' @param lambda_hat a estimator of parameter lambda (scale parameter)
#' @param n integer, number of observations in data
#'
#' @keywords internal
KS_test_full <- function(X_vec, theta_hat, lambda_hat, n) {
  UseMethod("KS_test_full")
}

#' @import data.table
#'
#' @keywords internal
KS_test_full.default <- function(X_vec, theta_hat, lambda_hat, n){
  U <- 1 - exp( -exp(theta_hat * (log(X_vec) - log(lambda_hat))))
  U_sort <- sort(U)
  expr1 <- max(1:n / n - U_sort)
  expr2 <- max(U_sort - 0:(n - 1) / n)

  sqrt(n) * max(expr1, expr2)
}

#' Cramer-von Mises test - Full Sample Case
#'
#' @param X_vec numeric, the values from uncensured sample
#' @param theta_hat a estimator of parameter theta (shape parameter)
#' @param lambda_hat a estimator of parameter lambda (scale parameter)
#' @param n integer, number of observations in data
#'
#' @keywords internal
#' @keywords internal
CM_test_full <- function(X_vec, theta_hat, lambda_hat, n) {
  UseMethod("CM_test_full")
}

#' @import data.table
#'
#' @keywords internal
CM_test_full.default <- function(X_vec, theta_hat, lambda_hat, n){
  U <- 1 - exp( -exp(theta_hat * (log(X_vec) - log(lambda_hat))))
  U_sort <- sort(U)

  sum((U_sort - (2 * (1:n) - 1) / (2 * n))^2) + 1 / (12 * n)
}

#' Liao and Shimokawa test - Full Sample Case
#'
#' @param X_vec numeric, the values from uncensured sample
#' @param n integer, number of observations in data
#'
#' @keywords internal
#' @keywords internal
LS_test_full <- function(X_vec, n) {
  UseMethod("LS_test_full")
}

#' @import data.table
#' @importFrom EWGoF LSEst
#'
#' @keywords internal
LS_test_full.default <- function(X_vec, n){
  lambda_hat_LS <- EWGoF::LSEst(X_vec)$eta
  theta_hat_LS <- EWGoF::LSEst(X_vec)$beta

  U <- 1 - exp( -exp(theta_hat_LS * (log(X_vec) - log(lambda_hat_LS))))
  U_sort <- sort(U)

  expr_1 <- (1:n) / n - U_sort
  expr_2 <- U_sort - (0:(n - 1)) / n

  1 / sqrt(n) * sum(pmax(expr_1, expr_2) / sqrt(U_sort * (1 - U_sort)))
}

#' Kolmogorov-Smirnov test - Right-Censored Sample Case
#'
#' @param Y_vec_dt a data table which contains value of transformed data statistic Y and
#' information of observation's censoring
#' @param n integer, number of observations in data
#'
#' @keywords internal
KS_test <- function(Y_vec_dt, n) {
  UseMethod("KS_test")
}

#' @import data.table
#'
#' @keywords internal
KS_test.default <- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[delta == 1, sort_Y_val]
  Est_KM_sort_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  expr1 <- 1 - exp(-exp(sort_Y_values))

  S1 <- max(Est_KM_sort_Y - expr1)
  S2 <- max(expr1 - Est_KM_sort_Y)

  max(S1, S2)
}

#' Cramer-von Mises test - Right_censored Sample Case
#'
#' @param Y_vec_dt a data table which contains value of transformed data statistic Y and
#' information of observation's censoring
#' @param n integer, number of observations in data
#'
#' @keywords internal
#' @keywords internal
CM_test <- function(Y_vec_dt, n) {
  UseMethod("CM_test")
}

#' @import data.table
#'
#' @keywords internal
CM_test.default <- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[, sort_Y_val]
  num_uncens <- sort_Y_vec_dt[, which(delta == 1)]
  d <- sort_Y_vec_dt[delta == 1, .N]

  Est_KM_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))

  expr1 <- sort_Y_values[num_uncens][2:d] - sort_Y_values[num_uncens][1:(d-1)]
  expr2 <- sort_Y_values[num_uncens][2:d]^2 - sort_Y_values[num_uncens][1:(d-1)]^2

  expr_full <- Est_KM_Y[num_uncens][2:d]^2 * expr1 - Est_KM_Y[num_uncens][2:d] * expr2

  sum(expr_full) * n + n / 3
}

#' Liao and Shimokawa test - Right-Censored Sample Case
#'
#' @param Y_vec_dt a data table which contains value of transformed data statistic Y and
#' information of observation's censoring
#' @param n integer, number of observations in data
#'
#' @keywords internal
#' @keywords internal
LS_test <- function(Y_vec_dt, n) {
  UseMethod("LS_test")
}

#' @import data.table
#'
#' @keywords internal
LS_test.default <- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[, sort_Y_val]
  Est_KM_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))

  d <- length(sort_Y_values)

  max_S1 <- sort_Y_vec_dt[, .("sort_Y_val" = 1:n / n * delta, delta)][delta == 1, sort_Y_val] - Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]]
  max_S2 <- Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]] - sort_Y_vec_dt[, .("sort_Y_val" = 0:(n-1) / n * delta, delta)][delta == 1, sort_Y_val]
  max_S <- pmax(max_S1, max_S2)

  expr1 <- sqrt(Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]] * (1 - Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]]))

  1 / sqrt(n) * sum(ifelse(max_S / expr1 == Inf, 0, max_S / expr1))
}
