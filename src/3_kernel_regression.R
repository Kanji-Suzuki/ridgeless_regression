# ---------------------------#
# Kernel Regression
# ---------------------------#
rm(list=ls())
gc()
gc()
library(pacman)
pacman::p_load(tidyverse, stringr, stringi, csv, ggplot2, mvtnorm, readxl, ExtDist,
               zoo, xts, parallel, beepr, fs, corpcor, doParallel, doRNG, progress, doSNOW,  matrixStats)
options(max.print=100000) 
options(digits=6) 
#----------------------------#
path_input <- "./data/intermediate"
path_output<- "./data/output"
source('R/function_kernel_regression.R')
#----------------------------#

### Read data. 
G_scaled <- readRDS( path_join(c(path_input, "G_scaled.rds")) )
y_scaled <- readRDS( path_join(c(path_input, "y_scaled.rds")) )
dates <- readRDS( path_join(c(path_input, "dates.rds")) )
time <- readRDS( path_join(c(path_input, "time.rds")) )

### set parameters
observable_feature_num_max <- 600
gamma <- 2
window_size <- 12 * 1
raw_feature_num <- ncol(G_scaled)


### Kernel Regression
# Out-of-Sample analysis
for( random_feature_method in c("RFFNormal", "RFFLaplace", "RFFCaucy", "ReLU", "erf") ){
  print( paste0("Processing: ", random_feature_method) )
  df_stats <- kernel_regression(y_scaled, G_scaled, window_size, random_feature_method) 
  saveRDS( df_stats, path_join( c(path_output, random_feature_method,
                                  paste0("df_kernel_stats.rds")) ))
}
