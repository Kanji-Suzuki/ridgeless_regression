# ---------------------------#
# function for plot
# ---------------------------#

plot_stats <- function(df_stats, df_kernel_stats){
  #' This function plots the stats from the prediction result
  #' 
  #' @param df_stats Prediction data
  #' @param df_kernel_stats Ridge Parameter
  #'
  #' @return ggplot obj
  
  plot1 <- df_stats %>% 
    ggplot(aes(x=c, y=R2, color=as.factor(ridge_param))) +
    geom_vline(xintercept = 1,
               color = "grey", linetype="dashed", alpha = 0.7)+
    geom_hline(yintercept = df_kernel_stats$R2,
               color = "black", linetype="longdash", alpha = 0.7)+
    geom_line() +
    labs(x="Complexity", y="R2", color="Ridge") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size = 14, "Helvetica"),
          axis.title = element_text(size = 18, "Helvetica", face = "bold"), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
    )
  
  if(max(abs(pull(df_stats, R2))) > 10000){
    plot1 <- plot1 +
      coord_cartesian(ylim = c(-20000, 10))
  }
  
  

  plot2 <- df_stats %>% 
    ggplot(aes(x=c, y=average_beta_hat_norm, color=as.factor(ridge_param))) +
    geom_vline(xintercept = 1,
               color = "grey", linetype="dashed", alpha = 0.7)+
    geom_line()  +
    labs(x="Complexity", y="Norm of the Coeffs", color="Ridge") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size = 14, "Helvetica"),
          axis.title = element_text(size = 18, "Helvetica", face = "bold"), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
    )
  

  plot3 <- df_stats %>% 
    ggplot(aes(x=c, y=timing_str_rt_average, color=as.factor(ridge_param))) +
    geom_vline(xintercept = 1,
               color = "grey", linetype="dashed", alpha = 0.7)+
    geom_hline(yintercept = df_kernel_stats$timing_str_rt_average,
               color = "black", linetype="longdash", alpha = 0.7)+
    geom_line() +
    labs(x="Complexity", y="Expected Return", color="Ridge") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size = 14, "Helvetica"),
          axis.title = element_text(size = 18, "Helvetica", face = "bold"), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
    ) 
  
  if(max(abs(pull(df_stats, timing_str_rt_average))) > 1){
    plot3 <- plot3 +
      coord_cartesian(ylim = c(-1, 1))
  }
  
  plot4 <- df_stats %>% 
    ggplot(aes(x=c, y=timing_str_rt_vol, color=as.factor(ridge_param))) +
    geom_vline(xintercept = 1,
               color = "grey", linetype="dashed", alpha = 0.7)+
    geom_hline(yintercept = df_kernel_stats$timing_str_rt_vol,
               color = "black", linetype="longdash", alpha = 0.7)+
    geom_line() +
    labs(x="Complexity", y="Volatility", color="Ridge") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size = 14, "Helvetica"),
          axis.title = element_text(size = 18, "Helvetica", face = "bold"), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
    )
  
  if(max(pull(df_stats, timing_str_rt_vol)) > 100){
    plot4 <- plot4 +
      coord_cartesian(ylim = c(0, 100))
  }
  
  plot5 <- df_stats %>% 
    ggplot(aes(x=c, y=sqrt(12)*timing_str_rt_SR, color=as.factor(ridge_param))) +
    geom_vline(xintercept = 1,
               color = "grey", linetype="dashed", alpha = 0.7)+
    geom_hline(yintercept = sqrt(12)*df_kernel_stats$timing_str_rt_SR,
               color = "black", linetype="longdash", alpha = 0.7)+
    geom_line() +
    labs(x="Complexity", y="Sharp Ratio", color="Ridge") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size = 14, "Helvetica"),
          axis.title = element_text(size = 18, "Helvetica", face = "bold"), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
    )
  
  if(max(abs(sqrt(12)*pull(df_stats, timing_str_rt_SR))) > 1){
    plot5 <- plot5 +
      coord_cartesian(ylim = c(-1, 1))
  }
  
  g <- grid.arrange(plot1, plot2, plot3, plot4, plot5, 
               ncol=2)
  
  return(g)
}


