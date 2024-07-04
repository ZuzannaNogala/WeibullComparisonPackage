#' Computing Quantiles for tests - statistic tests
#'
#' @param n integer, number of observations in full sample
#'
#' @import data.table
#'
#' @keywords internal
compute_test_values_uncens <- function(n){
  X <- rweibull(n, 1.5, 1)

  T_v_dt <- data.table("T_val" = X, "delta" = 1)

  theta_hat <- mle_theta_full_sample(n, T_v_dt[, T_val])
  lambda_hat <- mle_lambda_full_sample(n, T_v_dt[, T_val], theta_hat)

  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)

  list("S_n_1_1" =  S_n_1(a = 1, n, Y_v_dt),
       "S_n_1_2" =  S_n_1(a = 2, n, Y_v_dt),
       "S_n_1_5" =  S_n_1(a = 5, n, Y_v_dt),
       "S_n_2_1" =  S_n_2(a = 1, n, Y_v_dt),
       "S_n_2_2" =  S_n_2(a = 2, n, Y_v_dt),
       "S_n_2_5" =  S_n_2(a = 5, n, Y_v_dt),
       "KS" = KS_test_full(X, theta_hat, lambda_hat, n),
       "CM" = CM_test_full(X, theta_hat, lambda_hat, n),
       "LS" = LS_test_full(X, n))
}

#' Generate list of data.frames with statistic values to compute Quantiles
#'
#' @param m integer, amount of iterations of Parallel loop
#' @param n1 integer, amount of observations in full sample
#' @param ... integer, another amounts of observations in full sample
#'
#' @import parallel
#' @import data.table
#'
#'
#' @keywords internal
generation_quantiles_uncens <- function(m, n1, ...){
  n_params_list <- list(n1, ...)
  df <- as.data.frame(replicate(m, unlist(n_params_list)))
  cols_df <- paste0("uncesured_", unlist(n_params_list))

  n_cores <- detectCores() - 1

  list_of_values_to_quantiles <- mclapply(1:ncol(df),
                                          function(i){
                                            df <- sapply(1:nrow(df), function(row) compute_test_values_uncens(df[row, i]))
                                            colnames(df) <- cols_df
                                            df
                                            },
                                          mc.cores = n_cores)

  return(list_of_values_to_quantiles)
}

#' Find statistic's value for choosen type of data
#'
#' @param list_of_vals list of objects from class data.frame; result of \link{\code{generation_quantiles_uncens}} function
#' @param i integer, a position of test statistic in rownames
#' @param strColName characteristic, type of data (uncensored/censored and amount of observations)
#'
#'
#'
#' @keywords internal
getStatisticTestValue <- function(list_of_vals, i, strColName){
  unlist(lapply(list_of_vals, function(df) df[i, strColName]))
}

#' Find all statistics' value for choosen type of data and get data.frame
#'
#' @param list_of_vals list of objects from class data.frame; result of \link{\code{generation_quantiles_uncens}} function
#' @param strColName characteristic, type of data (uncensored/censored and amount of observations)
#'
#'
#'
#' @keywords internal
getDataFrameForSample <- function(list_of_vals, strColName){
  dt <- as.data.frame(sapply(1:9, function(i) getStatisticTestValue(list_of_vals, i, strColName)))
  names(dt) <- c("S_n_1_1", "S_n_1_2", "S_n_1_5", "S_n_2_1", "S_n_2_2", "S_n_2_5", "KS", "CM", "LS")
  dt
}

#' Find the estimation of quantile for choosen test value and type of data
#'
#' @param df data frame object...
#' @param strColName characteristic, type of data (uncensored/censored and amount of observations)
#' @param m description
#' @param alpha description
#'
#'
#' @keywords internal
getQuantileForTest <- function(df, strTestName, m = 10000, alpha = 0.05){
  test_vals <- df[, strTestName]
  sort(test_vals)[floor(m * (1 - alpha))]
}

#' Find the estimation of quantile for all test value and one type of data
#'
#' @param df data frame object...
#' @param m description
#' @param alpha description
#'
#'
#' @keywords internal
getQuantilesForEachTest <- function(df, m, alpha = 0.05){
  sapply(names(df), function(name) getQuantileForTest(df, name, m, alpha))
}

#' Compute Quantiles for each test and each type of data
#'
#' @param alpha data frame object...
#' @param m description
#' @param n1 description
#' @param ... description
#'
#' @import parallel
#'
#' @export
getQuantiles_DF <- function(alpha, m, n1, ...){
  UseMethod("getQuantiles_DF")
}
#' @export
getQuantiles_DF <- function(alpha, m, n1, ...){
  list_of_comp_stats <- generation_quantiles_uncens(m, n1, ...)
  rownames_df <- paste0("uncesured_", c(n1, ...))
  n_cores <- detectCores() - 1

  list_of_df <- mclapply(rownames_df, function(namerow){
    getDataFrameForSample(list_of_comp_stats, namerow)
  }, mc.cores = n_cores)

  list_of_all_quantiles <- lapply(list_of_df, function(df) getQuantilesForEachTest(df, m, alpha))

  df <- t(do.call(cbind, lapply(list_of_all_quantiles, data.frame)))
  rownames(df) <- rownames_df
  df
}

