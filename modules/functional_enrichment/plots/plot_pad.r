# experiments/e019/plots/plot_pad.r
# refer ggplot2 cheat sheet (ipad/useful/code)
# jul2023
# https://www.notion.so/woof7/Combining-variables-p-val-and-MCC-df64d93849cd4f688f15e677c4976beb?pvs=4
# see experiments/e019/create_final_output_tables_rank_aggregation.R

# RankProd
# https://bioconductor.org/install/
if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
# BiocManager::install(version = "3.14")
BiocManager::install(c("RankProd"))

# topklists - see:
# https://academic-oup-com.libproxy.ucl.ac.uk/bib/article/20/1/178/4091291?login=true#supplementary-data
# and Notion!
# https://www.notion.so/woof7/Combining-variables-p-val-and-MCC-df64d93849cd4f688f15e677c4976beb?pvs=4
install.packages("TopKLists")

library( RankProd )
library( TopKLists )
library( data.table )
library( tidyverse )

# https://rpkgs.datanovia.com/ggpubr/
library( ggpubr )

theme_set(
  theme_minimal() +
    theme(legend.position = "top")
  )

# -----------------------------------------------------
# GET DATA
# -----------------------------------------------------

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

# -----------------------------------------------------
# PLOT STUFF...
# -----------------------------------------------------
# data
p  	<- 		ggplot( df_final )
# p	<-		ggplot( filter( df_final, rank == 1 ) )

# hist mcc
# aesthetics mapping
pae <- p + 	aes( x = mcc, fill = factor( source ) ) 
# geom
p <- pae + 	geom_histogram( binwidth=.01, alpha=0.5, position = 'identity' ) 
# label
p <- p + 	labs( title = "Histogram for MCC", fill = "Group" )

# p vs mcc 
pae <- p + 	aes( y = mcc, x = -log10(p_value), color = factor( source ) )
p <- pae + 	geom_point() + geom_smooth( method = "lm" ) 
# p 	<- pae+ geom_point() + geom_smooth(method = "loess") 

# p vs mcc correlation  
pae <- p + 	aes( y = mcc, x = -log10(p_value), color = factor( source ) )
p <- pae + 	ggscatter( add = "reg.line", conf.int = TRUE, add.params = list(fill = "lightgray"), ggtheme = theme_minimal() )
p <- p + 	stat_cor( method = "pearson", label.x = 3, label.y = 30 ) 

# -----------------------------------------------------
# PLOT STUFF WITH PUBR EXTENSION
# https://rpkgs.datanovia.com/ggpubr/
# -----------------------------------------------------

# Groupped scatter plot
# https://rpkgs.datanovia.com/ggpubr/reference/ggscatter.html#details-1
#::::::::::::::::::::::::::::::::::::::::::::::::::::
df_filt <-  filter( df_final, rank == 1 )
p <- ggscatter(
		df_filt, y = "mcc", x = "p_value", color = "source" ,
  		palette = "jco",
  		add = "reg.line",
		cor.method = "pearson",
		font.label = c(22, "plain"),
  	) +
  		facet_wrap(~source) +
  		stat_cor(method = "pearson",  label.y = 1) +
  		stat_regline_equation(label.y = 1.03)
p <- ggpar(p, xscale = "log10")
p <- p + theme_pubr( base_size = 14 )
p
plot_name <- paste0( 'grouped_scatter_mcc_log10(pval)_', monet_type, '_', net_source, '.jpg' )
save_plot( plot_dir, plot_name )

# Simple scatter plot with correlation coefficient and
# regression line
ggscatter( df_final, y = mcc, x = -log10(p_value), add = "reg.line") +
  stat_cor( label.x = 3, label.y = 34 ) +
  stat_regline_equation(label.x = 3, label.y = 32)


# Scatter plot with correlation coefficient
#:::::::::::::::::::::::::::::::::::::::::::::::::# Load data
data("mtcars")
df <- mtcars
df$cyl <- as.factor(df$cyl)
sp <- ggscatter(df, x = "wt", y = "mpg",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 30)
sp + theme_pubr( base_size = 26 )


# try PCA
df_r1 <- filter(df_final, rank==1)
df_r1s <- df_r1 %>%
	mutate( label=map_chr(module, function(x) paste0('e019_M1_', x) ) ) %>%
	select( label, mcc, p_value )
df_r1s <- column_to_rownames( df_r1s, "label" )
pc1	<- prcomp( df_r1s, scale = TRUE )

#reverse the signs of the scores
# pc1$x <- -1*pc1$x

#display the first six scores
head(pc1$x)

biplot(pc1, scale=0)

pc1$sdev^2 / sum( pc1$sdev^2 )

# ----------------------------------------------------------------------------------------------------------------
# ranking / rank product
# see experiments/e019/create_final_output_tables_rank_aggregation.R
# ----------------------------------------------------------------------------------------------------------------
# # https://www.notion.so/woof7/Combining-variables-p-val-and-MCC-df64d93849cd4f688f15e677c4976beb?pvs=4

drp <- df_final %>% 
	group_by( module, source ) %>%
	mutate( row_rank_p = row_number( p_value ) ) %>%
	mutate( row_rank_mcc = row_number( desc( mcc ) ) ) %>%
	select( module, source, monet_type, term_id, rank, perc_rank, p_value, row_rank_p, mcc, row_rank_mcc, term_name )

# pivot wider
drpw <- pivot_wider( drp, names_from = source, values_from = c(row_rank_p, row_rank_mcc ) )

# test set
stype <- "GO:BP"
testdf <- filter(drp, module == 1 & source == stype ) %>%
		select( module, source, term_id, row_rank_p, row_rank_mcc )

# testdfw <- pivot_wider( testdf, names_from = source, values_from = c(row_rank_p, row_rank_mcc )  )

# -----------------------------------------
# TOPKLISTS - implemented: experiments/e019/create_final_output_tables_rank_aggregation.R
# see experiments/e019/create_final_output_tables_rank_aggregation.R
# -----------------------------------------
# /Users/ash/Dropbox/_iPad/UCL/useful_ipad/useful-stats/TopKLists.pdf
# https://www.notion.so/woof7/Combining-variables-p-val-and-MCC-df64d93849cd4f688f15e677c4976beb?pvs=4

# term-spaces
stype <- "GO:BP"
# stype <- "KEGG"
# stype <- "REAC"
mod <- 4

df_filt <- filter( drp, source == stype ) %>%
			select( "term_id" )
spaceall <- unique( df_filt$term_id )

df_filt <- filter( drp, source == stype & module == mod )

input_p		<- arrange( df_filt, row_rank_p )
input_mcc 	<- arrange( df_filt, row_rank_mcc )

input		<- list( input_p$term_id, input_mcc$term_id )
space		<- list( spaceall, spaceall )

agg.mc		<- MC( input, space )
agg.cemc	<- CEMC( input , space )

agg=list( MC1=agg.mc$MC1.TopK, MC2=agg.mc$MC2.TopK, MC3=agg.mc$MC3.TopK, CEMC=agg.cemc$TopK)
do.call(cbind, agg)

# --
# Above code for TopKLists MC1/2/3 etc now in function:aggregate_by_topk_lists(..) [experiments/e019/stats_functions.R]
source( file.path( root_dir,'stats_functions.R' ) )
df_test<- aggregate_by_topk_lists(df_final )

# remove old rank cols
# df_test <- select( df_test, -c(rank, perc_rank) )
# just keep top term for GO/KEGG/REAC (if any/all)
df_out <- filter(df_test, (term_id == top_go) | (term_id == top_kegg) | (term_id == top_reac) )
df_out <- select( df_out, -c( top_go, top_kegg, top_reac, p_value, mcc, row_rank_mcc, row_rank_p ) )

# widen just with essential cols top terms
df_outw <- pivot_wider( df_out, names_from = source, values_from = c(term_id, term_name) )
df_outw <- replace_na( df_outw, list( "term_id_REAC" = '-', "term_id_KEGG" = '-', "term_id_GO:BP" = '-' ) )
df_outw <- replace_na( df_outw, list( "term_name_REAC" = '-', "term_name_KEGG" = '-', "term_name_GO:BP" = '-' ) )

df_final_top <- select( df_final, -c(rank, perc_rank, parents) )

df_final_top_test <- inner_join( df_final_top, 
							df_outw,
							by = c( "module" = "module", "monet_type" = "monet_type" )
		)
