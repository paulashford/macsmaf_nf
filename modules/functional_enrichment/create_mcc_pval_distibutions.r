# experiments/e019/create_mcc_pval_distibutions.r
# 20 06 2023
# How to choose cut off for MCC and/or p-values... plot distibutions for ideas!
# install.packages("ggpubr")
# install.packages("ggpmisc")

library(data.table)
library(tidyverse)

# https://rpkgs.datanovia.com/ggpubr/
library(ggpubr)

theme_set(
  theme_minimal() +
    theme(legend.position = "top")
  )


# file root
root_dir <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019'
# plot out...
plot_dir <- file.path( root_dir, 'plots' )  
# parsing functions for gprofiler related 
source( file.path( root_dir,'gprofiler_functions.R' ) )
# stats-specific functions
source( file.path( root_dir,'stats_functions.R' ) )

# dict with all enrichment file paths and helper fns
source( file.path( root_dir, 'files.R' ) )

# monet_type <- 'M1'; net <- 'cpdb' 
# what to summarise?
# M1
net_source	<- 'cpdb'; monet_type <- 'M1'
# net_source	<- 'humanbase'; monet_type <- 'M1'
# net_source	<- 'string'; monet_type <- 'M1'
# R1
# net_source	<- 'cpdb'; monet_type <- 'R1'
# net_source	<- 'humanbase'; monet_type <- 'R1'
# net_source	<- 'string'; monet_type <- 'R1'
df_final <- load_final( final_dir = final_out_dir, monet_type = monet_type, net_source = net_source )

# Histogram distributions for:

# *MCC*
# plot Matthews correlation coefficients categorised by source
# https://vivdas.medium.com/visualizing-continous-data-with-ggplot2-in-r-2e4b7f433f67
hist_mcc <- ggplot( df_final, aes( x = mcc, fill = factor( source ) ) ) + 
	geom_histogram( binwidth=.01, alpha=0.5, position = 'identity' ) +
   	labs( title = "Histogram for MCC", fill = "Group" )
hist_mcc
plot_name <- paste0( 'hist_mcc_', monet_type, '_', net_source, '.jpg' )
save_plot( plot_dir, plot_name )

# *p-values*
# plot p-vals categorised by source
pvalplot <- ggplot( df_final, aes( x = p_value, fill = factor( source ) )) + 
	geom_histogram( binwidth = 0.0005, alpha=0.5, position = 'identity') +
   	labs(title  = "Histogram enrichment p_vals", fill="Group")
pvalplot
plot_name <- paste0( 'hist_pvals_', monet_type, '_', net_source, '.jpg' )
save_plot( plot_dir, plot_name )


View(df_final)
# http://www.sthda.com/english/articles/32-r-graphics-essentials/131-plot-two-continuous-variables-scatter-graph-and-alternatives/

# b <- ggplot(df_final, aes(x = mcc, y = p_value))
b <- ggplot(df_final, aes(x = log2(mcc), y = log2(p_value)))

# Scatter plot with regression line
b + geom_point()+
  geom_smooth(method = "lm") 
     
# Add a loess smoothed fit curve
b + geom_point()+
  geom_smooth(method = "loess") 

# b + geom_point(shape = 18)
# Add regression line and confidence interval
# Add correlation coefficient: stat_cor()
ggscatter(df_final, x = "mcc", y = "p_value",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson", 
           label.x = 3, label.y = 30) 

# Change color and shape by groups (source)
b + geom_point(aes(color = source, shape = source))+
  geom_smooth(aes(color = source, fill = source), method = "lm") +
  geom_rug(aes(color =source)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

# Remove confidence region (se = FALSE)
# Extend the regression lines: fullrange = TRUE
b + geom_point(aes(color = source, shape = source)) +
  geom_rug(aes(color =source)) +
  geom_smooth(aes(color = source), method = lm, 
              se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  ggpubr::stat_cor(aes(color = source), label.x = 3)

facet_mcc_p <- b + geom_point(aes(color = source, shape = source))+
  geom_smooth(aes(color = source, fill = source), 
              method = "lm", fullrange = TRUE) +
  facet_wrap(~source) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_bw()

facet_mcc_p

plot_dir	<- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/plots/'
# plot_file	<- 'mcc_dist_cat_source.jpg'
plot_file	<- 'facet_log2_mcc_log2_p.jpg'
# plot 		<- mccplot
plot 		<- facet_mcc_p	
filename	<- file.path(plot_dir, plot_file)
