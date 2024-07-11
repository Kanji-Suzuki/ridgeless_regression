# ---------------------------#
# function for the kernel regression
# ---------------------------#

ridge_rgr <- function(y_t, S_t, singular_value, ridge_param){
  #' This function conducts the ridge regression
  #' 
  #' @param y_t The return vec
  #' @param S_t The explanatory vars mat
  #' @param singular_value The result of the singular value decomposition
  #' @param ridge_param Ridge Parameter
  #'
  #' @return The ridge estimates
  
  window_size <- y_t %>% 
    length()
  observable_feature_num <- S_t %>% 
    ncol()

  # simple regime
  if(length(singular_value$d) == 1){
    temp_mat <- as.matrix(singular_value$d / (1*window_size + singular_value$d**2))
  } else{
    temp_mat <- diag(singular_value$d / (ridge_param*window_size + singular_value$d**2))
  }
  beta_hat <- singular_value$v %*% temp_mat %*% t(singular_value$u) %*% y_t
  
  beta_hat <- beta_hat %>% 
    as.vector()
  
  return(beta_hat)
}

calc_kernel <- function(random_feature_method, y_t, S_t, feature_i) {
  #' This function computes the kernel distance between the mat S_t and the vec feature_i
  #' 
  #' @param random_feature_method The type of kernel
  #' @param y_t The return vec
  #' @param S_t The explanatory vars mat
  #' @param feature_i The vector to compute the distance
  #'
  #' @return The kernel distance between the mat S_t and the vec feature_i
  
  window_size <- length(y_t)
  gamma <- 2
  
  if(random_feature_method == "RFFNormal"){
    kernel_distance <- exp( - (1/gamma^2) * colSums( (t(S_t) - feature_i)**2) )
  } else if(random_feature_method == "RFFLaplace"){
    kernel_distance <- exp( - colSums(abs(t(S_t) - feature_i)) )
  } else if(random_feature_method == "RFFCaucy"){
    kernel_distance <- colProds( 2 / (1 + (t(S_t) - feature_i)**2 ) )
  } else if(random_feature_method == "ReLU"){
    S_t <- S_t / sqrt(window_size)
    feature_i <- feature_i / sqrt(window_size)
    kernel_distance <- 1/pi * ( (S_t %*% feature_i) * ( pi - acos(pmin((S_t %*% feature_i) / sqrt(rowSums(S_t**2)) / sqrt(sum(feature_i**2)), 1)) ) +
                                          + sqrt( pmax( rowSums(S_t**2) * sum(feature_i**2) - (S_t %*% feature_i)**2, 0 ) ))
  } else if(random_feature_method == "erf"){
    # Williams, C. 1997
    S_t <- S_t / sqrt(window_size)
    feature_i <- feature_i / sqrt(window_size)
    kernel_distance <- 2/pi * asin( 2 * (S_t %*% feature_i) / sqrt(1+2*rowSums(S_t**2)) / sqrt(1+2*sum(feature_i**2)) )
  }
  
  return(kernel_distance)
}

calc_kernel_gram <- function(random_feature_method, y_t, S_t) {
  #' This function computes the kernel gram matrix
  #' 
  #' @param random_feature_method The type of kernel
  #' @param y_t The return vec
  #' @param S_t The explanatory vars mat
  #'
  #' @return Gram matrix
  
  # Compute the gram matrix of the kernel
  window_size <- length(y_t)
  raw_feature_num <- ncol(S_t)
  
  kernel_gram_matrix <- matrix(
    rep(0, window_size**2), nrow=window_size
  )
  
  for(i in 1:window_size){
    feature_i <- S_t[i, ]
    kernel_gram_matrix[i, ] <- calc_kernel(random_feature_method, y_t, S_t, feature_i)
  }

  return( kernel_gram_matrix )
}

kernel_regression <- function(y, S, window_size, random_feature_method) {
  #' For time and window_size, loop the out-of-sample estimation of the kernel regression and record the result
  #' 
  #' @param y The return vec
  #' @param S The factor mat
  #' @param window_size The window size
  #' @param random_feature_method The type of kernel
  #'
  #' @return df including the stats for each ridge_param
  
  time <- length(y)
  observable_feature_num <- ncol(S)
  
  # loop the sample period
  df_return <- data.frame()
  for(i in window_size:(time-1)){
    # predict y at i+1 using the last window_size records
    y_t <- y[(i-(window_size-1)):i]
    S_t <- S[(i-(window_size-1)):i, ]
    S_t_std <- sqrt(diag(cov(S_t)))
    S_t_scaled <- t(t(S_t) / ifelse(S_t_std == 0, 1, S_t_std))
    
    # Gram
    X_t <- calc_kernel_gram(random_feature_method, y_t, S_t_scaled)
    
    # compute the singular value decomposition. 
    # simple regime
    singular_value <- fast.svd( X_t )
    
    # list of beta 
    beta_hat <- ridge_rgr(y_t, X_t, singular_value, 0)
    
    # list of beta norm
    beta_hat_norm <- sum(beta_hat**2)
    
    # list of predicted return
    S_t_test <- S[i+1, ] / ifelse(S_t_std == 0, 1, S_t_std)
    rt_hat <- sum(beta_hat * calc_kernel(random_feature_method, y_t, S_t_scaled, S_t_test) )
    
    df_return <- df_return %>%
      bind_rows(
        data.frame(
          beta_hat_norm = beta_hat_norm,
          rt_hat = rt_hat,
          timing_str_rt = rt_hat * y[i+1],
          ridge_param = 0
        ) %>%
          mutate(
            rt = y[i+1],
            time = i
          )
      )
  }
  
  # compute the stats for each ridge_param
  df_stats <- calc_stats(df_return, 0, observable_feature_num, window_size)

  return( df_stats )
}

calc_stats <- function(df, ridge_param, observable_feature_num, window_size){
  #' This function computes the summary stats from the prediction result
  #' 
  #' @param df Prediction data
  #' @param ridge_param Ridge Parameter
  #' @param observable_feature_num The number of features
  #' @param window_size The window size
  #'
  #' @return The summary stats from the prediction result
  
  beta_hat_norm_ <- pull(df, beta_hat_norm)
  rt_hat_ <- pull(df, rt_hat)
  timing_str_rt_ <- pull(df, timing_str_rt)
  rt_ <- pull(df, rt)
  
  average_beta_hat_norm <- mean(beta_hat_norm_)
  R2 <- 1 - sum((rt_-rt_hat_)**2) / sum((rt_-mean(rt_))**2)
  timing_str_rt_average <- mean(timing_str_rt_)
  timing_str_rt_vol <- sd(timing_str_rt_)
  timing_str_rt_SR <- timing_str_rt_average / timing_str_rt_vol
  
  # Infomation Ratio
  model <- lm(rt_ ~ rt_hat_)
  beta <- coef(model)
  residuals <- residuals(model)
  
  return(
    data.frame(
      list("average_beta_hat_norm"=average_beta_hat_norm, "R2"=R2,
           "timing_str_rt_average"=timing_str_rt_average, "timing_str_rt_vol"=timing_str_rt_vol,
           "timing_str_rt_SR"=timing_str_rt_SR, 
           "alpha" = beta[1], "alpha_t" = summary(model)$coefficients[1, "t value"],
           "IR" = beta[1] / sd(residuals, na.rm = TRUE),
           "observable_feature_num"=observable_feature_num,
           "window_size"=window_size, "c" = observable_feature_num / window_size)
    )
  )
}
