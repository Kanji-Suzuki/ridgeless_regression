# ---------------------------#
# Data Cleaning
# ---------------------------#
rm(list=ls())
gc()
gc()
library(pacman)
pacman::p_load(tidyverse, stringr, stringi, csv, ggplot2, mvtnorm, readxl, ExtDist,
               zoo, xts, parallel, beepr, fs, corpcor, doParallel, doRNG, progress, doSNOW)
options(max.print=100000) 
options(digits=6) 
#----------------------------#
path_input <- "./data/raw"
path_output<- "./data/intermediate"
#----------------------------#

### Read data. Data is retrieved from KELLY et al.
df <- read.csv(path_join (c(path_input, "GYdata.csv")),
               header = FALSE)
dates <- df[, 1]
y <- df[, ncol(df)]
G <- df[, 2:(ncol(df)-1)] %>% 
  as.data.frame() %>% 
  mutate(V16 = lag(y))
time <- length(dates)

### Vol Filtering
vol_period_G <- 36
vol_period_y <- 12
G_scaled <- data.frame()
y_scaled <- rep(0, (time-vol_period_y) ) 
for(i in (vol_period_G+1):time){
  G_scaled <- G_scaled %>% 
    bind_rows( as.vector( 
      G[i, ] / apply(X=G[1:i, ], MARGIN = 2, FUN = function(x) as.vector( sqrt(var(x, na.rm = TRUE))) ) )
    )
}
G_scaled <- G_scaled %>% 
  as.matrix()
for(i in (vol_period_y+1):time){
  y_scaled[i-vol_period_y] <- y[i] / sqrt(mean((y[(i-vol_period_y):(i-1)])**2, na.rm=TRUE))
}

# cut off first vol_period_G recrods as the characteristics are only available after vol_period_G
y_scaled <- y_scaled[(vol_period_G-vol_period_y+1):length(y_scaled)]
dates <- dates[(vol_period_G+1):length(dates)]
time <- length(dates)

saveRDS(G_scaled,  path_join(c(path_output, "G_scaled.rds")) )
saveRDS(y_scaled,  path_join(c(path_output, "y_scaled.rds")) )
saveRDS(dates,  path_join(c(path_output, "dates.rds")) )
saveRDS(time,  path_join(c(path_output, "time.rds")) )

