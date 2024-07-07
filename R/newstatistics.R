#' @import data.table
#'
#' @keywords internal
Delta_k <- function(k, n, sort_Y_vec_dt){
  delta_vec <- sort_Y_vec_dt[, delta]
  if(k == 1) delta_vec[1] / n
  if(k == n) prod(((n - 1 : (n - 1)) / (n - 1 : (n - 1) + 1)) ^ delta_vec[1 : (n - 1)])
  else (delta_vec[k] / (n - k + 1)) * prod(((n - 1 : (k - 1)) / (n - 1 : (k - 1) + 1)) ^ delta_vec[1 : (k - 1)])
}

#' Computing value of new statistics - \eqn{S_{n,a}^{(1)}} and \eqn{S_{n,a}^{(2)}}
#'
#' The test statistic \eqn{S_{n,a}} is described by form:
#' \deqn{S_{n,a} = n\int_{-\infty}^{\infty} \left \lvert \sum ^n_{j=1} \Delta_j
#' \left[ite^{itY_j} + (1 - e^{Y_j})e^{itY_j}\right] \right \rvert ^2 w_a(t) \textrm{d}t,}
#' where \eqn{w_a(t)} is a weight function which depends on a tuning parameter \eqn{a > 0}.
#' If \eqn{w_a(t) = e^{-at^2}} the statistic simplifies to the form:
#' \deqn{S_{n,a}^{(1)} = \sqrt{\frac{\pi}{a}}\sum^n_{j=1}\sum^n_{k=1}\Delta_j
#' \Delta_k e^{\frac{-(Y_j-Y_k)^2}{4a}}\Bigl(-\frac{1}{4a^2} \left((Y_j-Y_k)^2
#' - 2a \right) + \frac{1 - e^{Y_j}}{a}(Y_j - Y_k) + (1 - e^{Y_j})(1 - e^{Y_k}) \Bigr)}
#' If \eqn{w_a(t) = e^{-a|t|}}: \deqn{S_{n,a}^{(2)} = n \sum^n_{j=1}\sum^n_{k=1}\Delta_j
#' \Delta_k \Bigl( \frac{-4a(3(Y_j - Y_k)^2) - a^2)}{((Y_j - Y_k)^2) + a^2)^3} +
#' \frac{8a(Y_j - Y_k)(1 - e^{Y_j})}{((Y_j - Y_k)^2) + a^2)^2} +
#' \frac{2a(1 - e^{Y_j})(1 - e^{Y_k})}{(Y_j - Y_k)^2) + a^2}\Bigr).}
#'
#' @param a numeric, tuning parameter > 0
#' @param n integer, number of observations in data
#' @param Y_vec_dt a data table which contains value of transformed statistic T (Y) and
#' information of observation's censoring
#'
#' @return numeric, a value of test statistics \eqn{S_{n,a}^{(1)}} or \eqn{S_{n,a}^{(2)}}
#'
#' @source E. Bothma, J. S Allison, I. J. H Vigagie. “New classes of tests for the Weibull
#' distribution using Stein’s method in the presence of random right censoring”.
#' Computational Statistics 37 (2022), pp. 1751–1770
#'
#' @rdname newStatistics
#' @keywords internal
S_n_1 <- function(a, n, Y_vec_dt){
  Y_val <- sort(Y_vec_dt[, Y_val])

  Delta <- sapply(1:n, function(k) Delta_k(k, n, sort_Y_vals(Y_vec_dt)))

  S_j <- sapply(1:n, function(i){
    Y_j <- Y_val[i]
    S_1 <- exp((- (Y_j - Y_val) ^ 2) / (4 * a))
    S_2 <- - ((Y_j - Y_val) ^ 2 - 2 * a) / (4 * a ^ 2)
    S_3 <- 2 * (1 - exp(Y_j)) * (Y_j - Y_val) / (2 * a) + (1 - exp(Y_j)) * (1 - exp(Y_val))

    sum(Delta * (S_1 * (S_2 + S_3)))
  })

  n * sqrt(pi / a) * sum(Delta * S_j)
}

#' @rdname newStatistics
#'
#' @keywords internal
S_n_2 <- function(a, n, Y_dt){
  Y_val <- sort(Y_dt[, Y_val])

  Delta <- sapply(1:n, function(k) Delta_k(k, n, sort_Y_vals(Y_dt)))

  S_j <- sapply(1:n, function(i){
    Y_j <- Y_val[i]
    S_1 <- (3 * (Y_j - Y_val) ^ 2 - a ^ 2) * (-4 * a)
    S_2 <- (Y_j - Y_val) ^ 2 + a ^ 2
    S_3 <- 8 * a * (Y_j - Y_val) * (1 - exp(Y_j))
    S_4 <- 2 * a * (1 - exp(Y_j)) * (1 - exp(Y_val))

    sum(Delta * (S_1 / (S_2) ^ 3 + S_3 / (S_2) ^ 2 + S_4 / S_2))
  })

  n * sum(Delta * S_j)
}

