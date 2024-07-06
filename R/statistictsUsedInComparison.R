globalVariables(c("sort_Y_val", "delta", "."))
#' Computing value of Komogorov-Smirnov test (\eqn{KS_n}) for full samples
#'
#' Classic one-sample Kolmogorov-Smirnov test (\eqn{KS_n}) applied to data without censoring.
#'
#' @param X_vec numeric, the values of observations from uncensored sample
#' @param theta_hat a estimator of parameter theta (shape parameter)
#' @param lambda_hat a estimator of parameter lambda (scale parameter)
#' @param n integer, amount of observations in data
#'
#' @return numeric,a value of test statistic \eqn{KS_n}
#'
#' @source R. B. D’Agostino, M. A. Stephens. "Goodness-of-Fit Techniques". New York: Marcel Dekker (1986)
#'
#' @keywords internal
KS_test_full <- function(X_vec, theta_hat, lambda_hat, n){
  U <- 1 - exp( -exp(theta_hat * (log(X_vec) - log(lambda_hat))))
  U_sort <- sort(U)
  expr1 <- max(1:n / n - U_sort)
  expr2 <- max(U_sort - 0:(n - 1) / n)

  sqrt(n) * max(expr1, expr2)
}

#' Computing value of Cramér von Mises test (\eqn{CM_n}) for full samples
#'
#' Classic Cramér von Mises test (\eqn{CM_n}) applied to data without censoring.
#'
#' @param X_vec numeric, the values of observations from uncensored sample
#' @param theta_hat a estimator of parameter theta (shape parameter)
#' @param lambda_hat a estimator of parameter lambda (scale parameter)
#' @param n integer, amount of observations in data
#'
#' @return numeric,a value of test statistic \eqn{CM_n}
#'
#' @source R. B. D’Agostino, M. A. Stephens. "Goodness-of-Fit Techniques". New York: Marcel Dekker (1986)
#'
#' @keywords internal
#' @keywords internal
CM_test_full<- function(X_vec, theta_hat, lambda_hat, n){
  U <- 1 - exp( -exp(theta_hat * (log(X_vec) - log(lambda_hat))))
  U_sort <- sort(U)

  sum((U_sort - (2 * (1:n) - 1) / (2 * n))^2) + 1 / (12 * n)
}

#' Computing value of Liao and Shimokawa test (\eqn{LS_n}) for full samples
#'
#' Classic Liao and Shimokawa test (\eqn{LS_n}) applied to data without censoring.
#'
#' @param X_vec numeric, the values of observations from uncensored sample
#' @param n integer, amount of observations in data
#'
#' @return numeric,a value of test statistic \eqn{LS_n}
#'
#' @details By contrast to another tests, Liao and Shimokawa's statistic is constructed base on
#' the least squares estimators of parameters instead of maximum likelihood ones.
#'
#' @source M. Liao, T. Shimokawa. “A new goodness-of-fit test for type-I extreme-value and
#' 2-parameter weibull distributions with estimated parameters”.
#' Journal of Statistical Computation and Simulation 64 (1999), pp. 23–48
#'
#' @importFrom EWGoF LSEst
#'
#' @keywords internal
LS_test_full <- function(X_vec, n){
  lambda_hat_LS <- EWGoF::LSEst(X_vec)$eta
  theta_hat_LS <- EWGoF::LSEst(X_vec)$beta

  U <- 1 - exp( -exp(theta_hat_LS * (log(X_vec) - log(lambda_hat_LS))))
  U_sort <- sort(U)

  expr_1 <- (1:n) / n - U_sort
  expr_2 <- U_sort - (0:(n - 1)) / n

  1 / sqrt(n) * sum(pmax(expr_1, expr_2) / sqrt(U_sort * (1 - U_sort)))
}

#' Computing value of Komogorov-Smirnov test (\eqn{KS_n^*}) for right-censored samples
#'
#' Classic one-sample Kolmogorov-Smirnov test (\eqn{KS_n^*}) which was modified to usage for data with
#' presence of censoring.
#'
#' @details
#' The test can be applied for full samples as well because of a Kaplan-Meier
#' estimator's universality. However, it can be lead too much simplification of algorithm to computing quantiles.
#' Hence, in that case it is recommended use \code{\link{KS_test_full}}.
#'
#' @source J.A. Koziol, S.B. Green. “A Cram ́er-von Mises statistic for randomly censored data”. Biometrika 63 (1976), pp. 465–474
#'
#' @param Y_vec_dt a data.table object, two-columns table which contains value of statistic \eqn{Y}
#' and information of observation's censoring (1 - uncensored, 0 - censored).
#' @param n integer, amount of observations in data
#'
#' @return numeric,a value of test statistic \eqn{KS_n^*}
#'
#' @import data.table
#'
#' @keywords internal
KS_test<- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[delta == 1, sort_Y_val]
  Est_KM_sort_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  expr1 <- 1 - exp(-exp(sort_Y_values))

  S1 <- max(Est_KM_sort_Y - expr1)
  S2 <- max(expr1 - Est_KM_sort_Y)

  max(S1, S2)
}

#' Computing value of Cramér von Mises test (\eqn{CM_n^*}) for right-censored samples
#'
#' Classic Cramér von Mises test (\eqn{CM_n^*}) which was modified by to usage for data with presence of
#' censoring.
#'
#' @details
#' The test can be applied for full samples as well because of a Kaplan-Meier
#' estimator's universality. However, it can be lead too much simplification of algorithm to computing quantiles.
#' Hence, in that case it is recommended use \code{\link{CM_test_full}}.
#'
#' @source J.A. Koziol, S.B. Green. “A Cram ́er-von Mises statistic for randomly censored data”. Biometrika 63 (1976), pp. 465–476
#'
#' @param Y_vec_dt a data.table object, two-columns table which contains value of statistic \eqn{Y}
#' and information of observation's censoring (1 - uncensored, 0 - censored).
#' @param n integer, amount of observations in data
#'
#' @return numeric,a value of test statistic \eqn{CM_n^*}
#'
#' @import data.table
#' @keywords internal
CM_test <- function(Y_vec_dt, n){
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

#' Computing value of Liao and Shimokawa test (\eqn{LS_n^*}) for right-censored samples
#'
#' Classic Liao and Shimokawa test (\eqn{LS_n^*}) which was modified by to usage for data with presence of
#' censoring.
#'
#' @details
#' The test can be applied for full samples as well because of a Kaplan-Meier
#' estimator's universality. However, it can be lead too much simplification of algorithm to computing quantiles.
#' Hence, in that case it is recommended use \code{\link{LS_test_full}}.
#'
#' @source N. Kim. “Goodness-of-fit tests for randomly censored Weibull distributions with estimated parameters”.
#' Communications for Statistical Applications and Methods 24 (2017), pp. 519–531
#'
#' @param Y_vec_dt a data.table object, two-columns table which contains value of statistic \eqn{Y}
#' and information of observation's censoring (1 - uncensored, 0 - censored).
#' @param n integer, amount of observations in data
#'
#' @return numeric,a value of test statistic \eqn{LS_n^*}
#'
#' @import data.table
#' @keywords internal
LS_test <- function(Y_vec_dt, n){
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
