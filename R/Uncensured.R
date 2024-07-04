#' Computing Quantiles for tests - statistic tests
#'
#' @param n integer, number of observations in full sample
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

# p <- generation_quantiles_uncens(2, 10, 20)


