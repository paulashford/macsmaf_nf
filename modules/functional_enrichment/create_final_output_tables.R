# experiments/e019/create_final_output_tables.R
# 15 06 2023
# NOTE: RankAggregation (for Top 1 lists) is here: experiments/e019/create_final_output_tables_rank_aggregation.R
# UPDATED 14 08 2023: Fix bug where not genes in R1 tables

# Following process_GOslims.R we can get the enrichments and the slims together for post-process and output
# some of this stuff was in process_GOslims.R before, when it really is a new analysis component.
# slims: 
# 	Dropbox/bioinf/MACSMAF/datasets/d014/goa_map_slim.sh
#	Dropbox/bioinf/MACSMAF/datasets/d014/v3

library( data.table )
library( tidyverse )

# file root
root_dir <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019'
# parsing functions for gprofiler related 
source( file.path( root_dir,'gprofiler_functions.R' ) )

# dict with all enrichment file paths and helper fns
source( file.path( root_dir,'files.R' ) )

# what to summarise?
# M1
net_source	<- 'cpdb'; monet_type <- 'M1'
# net_source	<- 'humanbase'; monet_type <- 'M1'
# net_source	<- 'string'; monet_type <- 'M1'
# R1
# net_source	<- 'cpdb'; monet_type <- 'R1'
# net_source	<- 'humanbase'; monet_type <- 'R1'
# net_source	<- 'string'; monet_type <- 'R1'

# Note calls below are for running the full data combination/cleaning/analysis steps
# Alternatively, can just load existing final tables here (for plots, analysis etc):
my_final_table <- load_final( final_dir =  file.path( root_dir, 'functional_enrichment_tables' ), net_source = net_source, monet_type = monet_type )

# ... or re-run all steps:
tib_slims	<-	load_goslims( slim_dir=go_slim_dir, net_source=net_source, monet_type=monet_type )
tib_enrich	<-  load_gp_enrichments( rda_path=gp_run_dir, net_source=net_source, monet_type=monet_type )
tib_modannot <- proc_modules( modules_dir = modules_dir, net_source = net_source, monet_type = monet_type ) 
# final enrichment file names / dir
out_dir 	<- file.path( root_dir, 'functional_enrichment_tables' )
out_file	<- paste0('func_enrichments_', monet_type, '_', net_source, '.tsv')

# combined gp enrichment with the subsequent GP:BP slims analysis (using PIR database)
comb_slim 	<- left_join( 	tib_enrich, 
							tib_slims, 
							by = c( 'func_module_number' = 'func_module_number' ), 
							na_matches = 'never', 
							keep = FALSE 
						) %>%
    rename( module = func_module_number )

# filter out very general terms
comb_filt	<- comb_slim %>%
	group_by( source ) %>%
	filter( term_size < ( 0.05 * effective_domain_size ) )

# stats
comb_filt <- comb_filt %>%
    group_by( module ) %>%
    mutate( perc_rank = percent_rank( p_value ) ) %>%
    mutate( cumelative_dist = cume_dist( p_value ) )

# Tidy cols and filter perc_rank 0.25 
comb_filt <- comb_filt %>%
    filter( perc_rank <= 0.25 )

# combine module genes to column
df_genes_mod <- tib_modannot %>%
	group_by( module, monet_type ) %>%
	mutate( genes = paste0( unique(gene), collapse = ";" ) ) %>%
	select( -gene ) %>%
	distinct()

# long form (adds mod annot gene list back!)
comb_long	<- inner_join( df_genes_mod, comb_filt, by = c( "module" = "module" ) ) %>%
	ungroup()
# comb_long	<- inner_join( tib_modannot, comb_filt, by = c( "module" = "module" ) )
# 	ungroup() %>%
# 	mutate( monet_type = monet_type )

# Create overall rank per module
# NOTE: RankAggregation (for Top 1 lists) is here: experiments/e019/create_final_output_tables_rank_aggregation.R
# comb_summary <- comb_summary %>%
comb_summary <- comb_long %>%
	group_by( monet_type, module ) %>%
    mutate( rank = rank( p_value, ties.method = "average" ) ) %>%
    arrange( module, rank )

# add tp, fp etc cols
comb_summary2 <- add_tp_tn_fp_fn_MC( comb_summary, inc_mcc = FALSE )

# MCC Matthew's Correlation Coefficient calcs
# tried various many purrr and dplyr and data.table operations but can't seem to get any success
# updating 1 column based on these nice new simple columns (fp, tp etc) using my straightforward function calc_mcc
# reverting to data.table /for loop
dt_comb_summary2 <- data.table( comb_summary2 )
for( i in 1:nrow( dt_comb_summary2 ) ){
	dbl_mcc <- calc_mcc(
						dt_comb_summary2[i]$tp,
						dt_comb_summary2[i]$fn,
						dt_comb_summary2[i]$fp,
						dt_comb_summary2[i]$tn
					)
	dt_comb_summary2[ i, mcc := dbl_mcc ]

}

# finalise
comb_summary_final <- as_tibble( dt_comb_summary2 )
comb_summary_final <- relocate( comb_summary_final, genes, .after = cumelative_dist )

# write full table
write_tibble( tib = comb_summary_final, out_dir = out_dir, file_name = out_file, sep = '\t' )
# rank1 only
# NOTE: RankAggregation (for Top 1 lists) is here: experiments/e019/create_final_output_tables_rank_aggregation.R
# write_tibble( 	tib = filter( comb_summary_final, rank == 1 ), 
# 				out_dir = out_dir, 
# 				file_name = paste0('rank1_',out_file), 
# 				sep = '\t' )


# ------------------------------------------------------------------------------------
# BUG: col order issue for temp R1s with genes='-' - check same number (ok 29)
library(data.table)
library(tidyverse)
# file root
root_dir <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019'
bug_dir 	<- file.path( root_dir, 'functional_enrichment_tables', 'col_order_bug' )
out_dir 	<- file.path( root_dir, 'functional_enrichment_tables' )

# parsing functions for gprofiler related 
source( file.path( root_dir,'gprofiler_functions.R' ) )

# dict with all enrichment file paths and helper fns
source( file.path( root_dir,'files.R' ) )

# sort r1 order in these temp versions where genes='-'
# monet_type 	<- 'R1'; net_source <- 'cpdb'
# monet_type 	<- 'R1'; net_source <- 'humanbase'
monet_type 	<- 'R1'; net_source <- 'string'

# read file and correct
out_file		<- paste0('func_enrichments_', monet_type, '_', net_source, '.tsv')
comb_summary_final <- get_tibble( root_dir = bug_dir, file_name = out_file )
rel_comb_summary_final <- relocate( comb_summary_final, module, monet_type, experiment_info, .before = experiment_info )
# write full table
write_tibble( tib = rel_comb_summary_final, out_dir = out_dir, file_name = out_file, sep = '\t' )

# read file and correct rank1
rank1_file 	<- paste0('rank1_',out_file)
comb_summary_final <- get_tibble( root_dir = bug_dir, file_name = rank1_file )
rel_comb_summary_final <- relocate( comb_summary_final, module, monet_type, experiment_info, .before = experiment_info )
# write full table
write_tibble( tib = rel_comb_summary_final, out_dir = out_dir, file_name = rank1_file, sep = '\t' )


# ------------------------------------------------------------------------------------

# looking for something else?
# experiments/e019/create_final_output_tables_archived.R