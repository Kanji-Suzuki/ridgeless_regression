# ---------------------------#
# Replicate the Random Features Excercise
# ---------------------------#
rm(list=ls())
gc()
gc()
library(pacman)
pacman::p_load(tidyverse, stringr, stringi, csv, ggplot2, mvtnorm, readxl, ExtDist,
               zoo, xts, parallel, beepr, fs, corpcor, doParallel, doRNG, progress, doSNOW, docstring)
options(max.print=100000) 
options(digits=6) 
#----------------------------#
path_input <- "./data/intermediate"
path_output<- "./data/output"
source('R/function_random_feature.R')
#----------------------------#

### Read data. 
G_scaled <- readRDS( path_join(c(path_input, "G_scaled.rds")) )
y_scaled <- readRDS( path_join(c(path_input, "y_scaled.rds")) )
dates <- readRDS( path_join(c(path_input, "dates.rds")) )
time <- readRDS( path_join(c(path_input, "time.rds")) )

### set parameters
observable_feature_num_max <- 600
gamma <- 2
num_of_simul <- 1000
window_size <- 12 * 1
ridge_param_list <- c(10**-3, 10**-2, 10**-1, 1, 10**1, 10**2, 10**3)
observable_feature_num_list <- c( (1:15)*2, (1:(observable_feature_num_max/50-1))*50, observable_feature_num_max)
random_feature_method <- "erf"
raw_feature_num <- ncol(G_scaled)

# if directory doesn't exst, create
if(!dir.exists( path_join( c(path_output, random_feature_method )) )){
  dir.create( path_join( c(path_output, random_feature_method )) )
}

### Simulation and Estimation. Random Features.
start_time <- Sys.time()
# progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed",
  total = num_of_simul, 
  clear = FALSE
)
progress <- function(n){
  pb$tick()
} 
opts <- list(progress = progress)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
registerDoSNOW(cl)
# rv seed
registerDoRNG(seed = 1)
result <- foreach(i = 1:num_of_simul,
                  .packages = c('dplyr', 'corpcor', 'purrr', 'tidyr', 'fs', 'stringr', 'ExtDist', 'progress'),
                  .export = c('ridge_rgr', 'calc_stats'),
                  .options.snow = opts) %dopar% {

  # generate random linear converter to construct the factors
  if( str_detect(random_feature_method, "RFF") ){
    S <- generate_RFF( observable_feature_num_max, time, gamma, G_scaled, random_feature_method)
  } else if(random_feature_method == "ReLU"){
    # sd is 1/sqrt(raw_feature_num) for scaling. 
    W <- (rnorm( raw_feature_num * observable_feature_num_max, mean = 0, sd = 1/sqrt(raw_feature_num))) %>%
      matrix(nrow = observable_feature_num_max)
    # Construct factors by 1-layer NN with ReLU activaation
    S <- sqrt(2) * pmax(t(t(G_scaled %*% t(W))), 0) 
  } else if(random_feature_method == "erf"){
    # sd is 1/sqrt(raw_feature_num) for scaling. If not, most features become 1 or -1.
    W <- (rnorm( raw_feature_num * observable_feature_num_max, mean = 0, sd = 1/sqrt(raw_feature_num))) %>%
      matrix(nrow = observable_feature_num_max)
    S <- 2 * pnorm( t(t(G_scaled %*% t(W))) * sqrt(2), mean = 0, sd = 1) - 1
  } else{
    stop("Invalid input in observable_feature_num_list")
  }
  
  # Out-of-Sample analysis
  df_stats <- data.frame()
  for(observable_feature_num in observable_feature_num_list ){
    cat(paste0("# --------------- Processing: ", observable_feature_num, " ---------------# \n"))
    df_stats <- df_stats %>% 
      bind_rows(
        RF_regression(y_scaled, S, observable_feature_num, ridge_param_list, window_size)
      )
  }
  saveRDS( df_stats, path_join( c(path_output, random_feature_method,
                                  paste0("df_stats_", i, ".rds")) ))

  return(df_stats)
}
stopCluster(cl)
end_time <- Sys.time()
print(end_time - start_time)


### aggregate across simulation
df_stats <- data.frame()
for(i in 1:num_of_simul){
  df_stats <- df_stats %>%
    bind_rows(
      readRDS(path_join( c(path_output, random_feature_method,
                           paste0("df_stats_", i, ".rds"))) ) %>%
        mutate(simul_i = i)
    )
}

threshold <- 0.01
df_stats <- df_stats %>%
  group_nest(c, observable_feature_num, ridge_param, window_size) %>%
  mutate(stats = map(data, ~{aggreagate_stats(., threshold)})) %>%
  unnest(stats) %>% 
  select(-data)
saveRDS( df_stats, path_join( c(path_output, random_feature_method,
                                paste0("df_stats.rds"))))

