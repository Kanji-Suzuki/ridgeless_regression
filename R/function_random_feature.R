# ---------------------------#
# function for the random feature regression
# ---------------------------#


generate_RFF <- function( observable_feature_num_max, gamma, G_scaled, 
                         random_feature_method ){
  #' This function generates the Random Fourier Feature.
  #' 
  #' @param observable_feature_num_max The number of features to generate
  #' @param gamma Band width
  #' @param G_scaled The raw feature
  #' @param random_feature_method Type of the random weight
  #'
  #' @return The raw feature matrix
  
  
  raw_feature_num <- ncol(G_scaled)
  time <-  nrow(G_scaled)
  
  if(random_feature_method == "RFFNormal"){
    W <- rnorm( raw_feature_num * observable_feature_num_max / 2 , mean = 0, sd = 1) %>%
      matrix(nrow = observable_feature_num_max / 2)
  } else if(random_feature_method == "RFFLaplace"){
    W <- rcauchy( raw_feature_num * observable_feature_num_max / 2 , location = 0, scale = 1) %>%
      matrix(nrow = observable_feature_num_max / 2)
  } else if(random_feature_method == "RFFCaucy"){
    W <- rLaplace(raw_feature_num * observable_feature_num_max / 2 , mu = 0, b = 1) %>%
      matrix(nrow = observable_feature_num_max / 2)
  }
  
  # Construct factors by RFF
  # The first half cols are sin(), and the half left is cos().
  S <- matrix(rep(0, time * observable_feature_num_max ), nrow = time)
  S[, 1:(observable_feature_num_max / 2)*2 - 1 ] <- cos(gamma * (G_scaled %*% t(W))) 
  S[, 1:(observable_feature_num_max / 2)*2 ] <- sin(gamma * (G_scaled %*% t(W)))
  
  return(S)
}


ridge_rgr <- function(y_t, S_t, singular_value, ridge_param){
  #' This function conducts the ridge regression
  #' 
  #' @param y_t The return vec
  #' @param S_t The explanatory vars mat
  #' @param singular_value The result of the singular value decomposition
  #' @param ridge_param Ridge Parameter
  #'
  #' @return The ridge estimates
  
  window_size <- S_t %>% 
    nrow()
  observable_feature_num <- S_t %>% 
    ncol()
  
  if(observable_feature_num <= window_size){
    # simple regime
    if(length(singular_value$d) == 1){
      temp_mat <- as.matrix(singular_value$d / (1*window_size + singular_value$d**2))
    } else{
      temp_mat <- diag(singular_value$d / (ridge_param*window_size + singular_value$d**2))
    }
    beta_hat <- singular_value$v %*% temp_mat %*% t(singular_value$u) %*% y_t
  }else{
    # complex regime
    if(length(singular_value$d) == 1){
      temp_mat <- as.matrix(1 / sqrt( singular_value$d * window_size ))
      temp_mat2 <- as.matrix( 1 / (ridge_param + singular_value$d) )
    } else{
      temp_mat <- diag( 1 / sqrt( singular_value$d * window_size ) )
      temp_mat2 <- diag( 1 / (ridge_param + singular_value$d) )
    }
    W <- t(S_t) %*% singular_value$u %*% temp_mat
    beta_hat <- crossprod(t(W), crossprod(temp_mat2,  crossprod(W,  (t(S_t) %*% y_t / window_size))))
  }
  
  beta_hat <- beta_hat %>% 
    as.vector()
  return(beta_hat)
}



RF_regression <- function(y, S, observable_feature_num, ridge_param_list, window_size){
  #' For time and window_size, loop the out-of-sample estimation and record the result
  #' 
  #' @param y The return vec
  #' @param S The factor mat
  #' @param observable_feature_num The result of the singular value decomposition
  #' @param ridge_param_list Ridge parameter list
  #' @param window_size The window size
  #'
  #' @return df including the stats for each ridge_param
  
  time <- length(y)
  
  # loop the sample period
  df_return <- data.frame()
  for(i in window_size:(time-1)){
    # predict y at i+1 using the last window_size records
    y_t <- y[(i-(window_size-1)):i]
    S_t <- S[(i-(window_size-1)):i, 1:observable_feature_num]
    S_t_std <- sqrt(diag(cov(S_t)))
    S_t_scaled <- t(t(S_t) / ifelse(S_t_std == 0, 1, S_t_std))
    
    # compute the singular value decomposition. 
    # This is common across all elements of ridge_param_list.
    if(observable_feature_num <= window_size){
      # simple regime
      singular_value <- fast.svd( S_t_scaled )
    }else{
      # complex regime
      singular_value <- fast.svd( S_t_scaled %*% t(S_t_scaled) / window_size )
    }
    
    # list of beta for each ridge_param
    beta_hat_list <- map(ridge_param_list,
                              ~ridge_rgr(y_t, S_t_scaled, singular_value, .x))
    # print(beta_hat_list)
    # list of beta norm
    beta_hat_norm_list <- map(beta_hat_list,
                              ~{sum(.x**2)}) %>%
      unlist()
    # list of predicted return
    rt_hat_list <- map(beta_hat_list,
                       ~ sum(.x * (S[i+1, 1:observable_feature_num] / ifelse(S_t_std == 0, 1, S_t_std)))) %>%
      unlist()

    df_return <- df_return %>%
      bind_rows(
        data.frame(
          beta_hat_norm = beta_hat_norm_list,
          rt_hat = rt_hat_list,
          timing_str_rt = rt_hat_list * y[i+1],
          ridge_param = ridge_param_list
        ) %>%
          mutate(
            rt = y[i+1],
            time = i
          )
      )
  }
  
  # compute the stats for each ridge_param
  df_stats <- df_return %>% 
    group_nest(ridge_param) %>% 
    mutate(stats = map2(data, ridge_param,
                       ~calc_stats(.x, .y,
                                      observable_feature_num, window_size))) %>% 
    unnest(stats) %>% 
    select(-data)
  
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

aggreagate_stats <- function(df_stats, threshold){
  #' This function aggregates the summary stats across simulations
  #' 
  #' @param df_stats Statistics data
  #' @param threshold threshold to cut off the estreme values
  #'
  #' @return The summary stats from the prediction result
  
  df_stats <- df_stats %>% 
    drop_na()
  average_beta_hat_norm_v <- df_stats %>% pull(average_beta_hat_norm)
  R2_v <- df_stats %>% pull(R2)
  timing_str_rt_average_v <- df_stats %>% pull(timing_str_rt_average)
  timing_str_rt_vol_v <- df_stats %>% pull(timing_str_rt_vol)
  timing_str_rt_SR_v <- df_stats %>% pull(timing_str_rt_SR)
  
  df_return <- data.frame(
    average_beta_hat_norm = mean(average_beta_hat_norm_v[quantile(average_beta_hat_norm_v, threshold) <= average_beta_hat_norm_v & average_beta_hat_norm_v <= quantile(average_beta_hat_norm_v, 1-threshold)]),
    R2 = mean(R2_v[quantile(R2_v, threshold) <= R2_v & R2_v <= quantile(R2_v, 1-threshold)]),
    timing_str_rt_average = mean(timing_str_rt_average_v[quantile(timing_str_rt_average_v, threshold) <= timing_str_rt_average_v & timing_str_rt_average_v <= quantile(timing_str_rt_average_v, 1-threshold)]),
    timing_str_rt_vol = mean(timing_str_rt_vol_v[quantile(timing_str_rt_vol_v, threshold) <= timing_str_rt_vol_v & timing_str_rt_vol_v <= quantile(timing_str_rt_vol_v, 1-threshold)]),
    timing_str_rt_SR = mean(timing_str_rt_SR_v[quantile(timing_str_rt_SR_v, threshold) <= timing_str_rt_SR_v & timing_str_rt_SR_v <= quantile(timing_str_rt_SR_v, 1-threshold)])
  )
  
  return(df_return)
}


