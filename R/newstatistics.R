#' Computing value of statistic \eqn{\Delta_k}
#'
#' Statistic represents the size of jump in \eqn{G_n(T_{(k)})} and can be described by a formula
#' \deqn{\Delta_k = G_n(T_{(k)}) - \lim_{t \rightarrow T_{(k)}} G_n(t), \ k = 1, \ldots n} This difference can be simplified to
#' calculable expression and the function \code{Delta_k} used it to computes value of statistic.
#'
#' @param k integer, k-th observation in sorted data
#' @param n integer, amount of observations in data
#' @param sort_Y_vec_dt data.table object, two-columns table which contains value of sorted statistic \eqn{Y}
#' and information of observation's censoring (1 - uncensored, 0 - censored).
#'
#' @return numeric, a value of statistic \eqn{\Delta_k}
#'
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
#'
#'
#' @param a numeric, the value of statistic's parameter
#' @param n integer, number of observations in data
#' @param Y_vec_dt a data table which contains value of transformed statistic T (Y) and
#' information of observation's censoring
#'
#' @return numeric, a value of test statistic \eqn{S_{n,a}^{(1)}} or \eqn{S_{n,a}^{(2)}}
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

