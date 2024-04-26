# experiments/e019/create_final_output_tables_rank_aggregation.R
# 27 07 2023
# UPDATE 14 08 2023
# 1) Fix missing R1 gene lists bug
# 	archived vers: e019/functional_enrichment_tables/archived_r1_enrich_missing_genes
# 2) HUGO update

# This uses output from: experiments/e019/create_final_output_tables.R
# Assumes "final" functional enrichment files are available
# This will run RankAggregation, favouring MCC rank if tied rank 1/2, using MC3 method from
# TopKLists, and report one function per module having a sig. functional enrichment 
# as combo of GO:BP KEGG and Reac (1, 2 or 3 of these)

# TopKLists - see:
# https://academic-oup-com.libproxy.ucl.ac.uk/bib/article/20/1/178/4091291?login=true#supplementary-data
# /Users/ash/Dropbox/_iPad/UCL/useful_ipad/useful-stats/TopKLists.pdf
# https://www.notion.so/woof7/Combining-variables-p-val-and-MCC-df64d93849cd4f688f15e677c4976beb?pvs=4

# install.packages("TopKLists")
library( TopKLists )
library( data.table )
library( tidyverse )

# file root
root_dir <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019'
# parsing functions for gprofiler related 
source( file.path( root_dir,'gprofiler_functions.R' ) )
# dict with all enrichment file paths and helper fns
source( file.path( root_dir,'files.R' ) )
# stats fns
source( file.path( root_dir, 'stats_functions.R' ) )
# out files
out_dir 	<- file.path( root_dir, 'functional_enrichment_tables' )

# what to summarise?

# cpdb
net_source	<- 'cpdb'; monet_type <- 'M1'
net_source	<- 'cpdb'; monet_type <- 'R1'

# humanbase
net_source	<- 'humanbase'; monet_type <- 'M1'
net_source	<- 'humanbase'; monet_type <- 'R1'

# string
net_source	<- 'string'; monet_type <- 'M1'
net_source	<- 'string'; monet_type <- 'R1'

# has rows for KEGG, REAC and GO with all p-vals etc (1-3 rows per module)
detail_out_file	<- paste0( 'func_rank1_agg_detail_', monet_type, '_', net_source, '.tsv' )

# 1-row-per module summary
summary_out_file	<- paste0( 'func_rank1_agg_summary_', monet_type, '_', net_source, '.tsv' )
# ... with hugos and hugo descriptions
summary_out_file_hugo	<- paste0( 'func_rank1_agg_summary_hugo_', monet_type, '_', net_source, '.tsv' )

######################################
# LOAD EXISTING PREVIOUS VERSIONS
# Just load existing final tables here (for plots, analysis etc)
######################################
# RANK1 AGG SUMMARY
df_final 	<- load_final_rank1( final_dir =  file.path( root_dir, 'functional_enrichment_tables' ), net_source = net_source, monet_type = monet_type, summary = TRUE )

# RANK1 AGG DETAIL
df_finaldet	<- load_final_rank1( final_dir =  file.path( root_dir, 'functional_enrichment_tables' ), net_source = net_source, monet_type = monet_type, summary = FALSE )
# Unranked 'pre-final'
df_prefinal <- load_final( final_dir =  file.path( root_dir, 'functional_enrichment_tables' ), net_source = net_source, monet_type = monet_type )

######################################
# GENE NAMES TO HUGO UPDATE
# note: needed to fix missing gene list in r1 bug too!
# archived vers: e019/functional_enrichment_tables/archived_r1_enrich_missing_genes
# 14 08 2023
######################################
# (a) output unique genes for conversion wiht g:convert
# https://biit.cs.ut.ee/gprofiler/convert
	# # separate gene names to rows (just do for summary)
	# df_final_longgenes <- separate_rows( df_final, genes, sep=";" )

	# # Suggest doing humanbase, string etc sep as they use diff gene ids!
	# # 1st net
	# # all_genes <- df_final_longgenes$genes
	# # subsequent...
	# all_genes <- c( all_genes, df_final_longgenes$genes )

	# # finalise and save for conversion
	# unq_genes <- unique( all_genes )
	# hugo_dir <- file.path( root_dir, 'functional_enrichment_tables', 'hugo_gene_convert' )
	# hugo_file <- paste0('orig_gene_id_', net_source , '.csv')
	# # write file
	# write.csv( unq_genes, file.path( hugo_dir, hugo_file ), row.names=FALSE, quote=FALSE )
	# g:convert (online)
	# ...
# (b) read/join conversions!
	hugo_dir <- file.path( root_dir, 'functional_enrichment_tables', 'hugo_gene_convert' )
	converted_file <- paste0('gProfiler_hsapiens_', net_source , '.csv')
	df_gconv <- get_tibble( root_dir = hugo_dir, file_name = converted_file, sep = "," )
	# long-form of unconverted enrich file
	df_final_longgenes <- separate_rows( df_final, genes, sep=";" )
	# oin conv hugo
	df_gconv$initial_alias <- as.character(df_gconv$initial_alias)
	df_final_long_conv	<- left_join( df_final_longgenes, df_gconv, by = c( "genes" = "initial_alias" ) ) %>%
		select( -c(namespace, name))
	
	# collapse gene names cols
	df_final_hugo <- df_final_long_conv %>%
		group_by( module, monet_type ) %>%
		mutate( genes_orig = paste0( unique(genes), collapse = ";" ) ) %>%
		mutate( genes_hugo = paste0( unique(converted_alias), collapse = ";" ) ) %>%
		mutate( genes_descriptions = paste0( unique(description), collapse = " ~ " ) ) %>%
		select( -c(genes,converted_alias, description ) )%>%
		distinct()

	# Save summary file with hugos
	write_tibble( tib =df_final_hugo, out_dir = out_dir, file_name = summary_out_file_hugo, sep = '\t' )


######################################
# RE-RUN THE RANK AGGREGATION STEPS
# for new agg calcs and regeneration
######################################
# Run rank agg for p-val, mcc
df_rankagg <- aggregate_by_topk_lists( df_prefinal )

# Join rank agg to df_prefinal
df_prefinal_rankagg 	<- inner_join( 	df_prefinal,
										select(df_rankagg, -c( p_value, mcc, term_name ) ),
										by = c( "experiment_info" = "experiment_info",
												"monet_type" = "monet_type",
												"module" = "module", 
												"source" = "source",
												"term_id" = "term_id"
											)
							)

# Remove rows if not a rank 1 term...
df_final_rankagg 	<- filter( df_prefinal_rankagg, 
							( term_id == top_go ) | ( term_id == top_kegg ) | ( term_id == top_reac ) )
# tidy
df_final_rankagg 	<- relocate( df_final_rankagg, goslim_pir, genes, .after = top_reac )
# not nec here - wide form gets 'GO:BP' in col name which makes processing awks!
# df_final_rankagg 	<- rename_with( df_final_rankagg, ~ gsub( ":", "_", .x, fixed = TRUE ) )

# Save detail file
write_tibble( tib = df_final_rankagg, out_dir = out_dir, file_name = detail_out_file, sep = '\t' )

# This will have 1,2 or 3 rows per module dep whether GO:BP, KEGG, Reac available...
# This form removes all enrichment values and pivots terms to wide form
df_final_rankagg_simple <- select( 	df_final_rankagg, 
									c( "module", "monet_type", "experiment_info",
										"source", "term_id", "term_name", 
										"goslim_pir" , "genes" 
									)
								)

# pivot wide
df_final_rankagg_wide <- pivot_wider( df_final_rankagg_simple, 
										names_from = source, 
										values_from = c( term_id, term_name ) 
									)
# tidy
# wide form gets 'GO:BP' in col name which makes processing awks!
df_final_rankagg_wide <- rename_with( df_final_rankagg_wide, ~ gsub( ":", "_", .x, fixed = TRUE ) )
df_final_rankagg_wide <- relocate( df_final_rankagg_wide, goslim_pir, genes, .after = term_name_KEGG )
# function summary of all non-NA func names:
df_final_rankagg_wide <- unite( df_final_rankagg_wide, term_name_GO_BP, term_name_KEGG, term_name_REAC, 
									col = "function_summary", 
									sep = " / ", 
									na.rm = TRUE )
df_final_rankagg_wide <- relocate( df_final_rankagg_wide, function_summary, .after = experiment_info )

# Save summary file
write_tibble( tib = df_final_rankagg_wide, out_dir = out_dir, file_name = summary_out_file, sep = '\t' )


# not nec - will use unite(...na.rm=TRUE)
# df_final_rankagg_wide <- replace_na( df_final_rankagg_wide, 
# 										list( "term_id_REAC" = '-', "term_id_KEGG" = '-', "term_id_GO:BP" = '-' ) )
# df_final_rankagg_wide <- replace_na( df_final_rankagg_wide, 
# 										list( "term_name_REAC" = '-', "term_name_KEGG" = '-', "term_name_GO:BP" = '-' ) )

# # all cols avail from colnames(df_final_rankagg)
# "module", "monet_type", "experiment_info",
# "significant", "p_value", "term_size",
# "query_size","intersection_size", "precision",
# "recall","term_id", "source", 
# "term_name", "effective_domain_size" "source_order", 
# "parents", "goslim_pir","perc_rank",
# "cumelative_dist", "genes", "rank", 
# "P", "N", "tp", 
# "fp","fn","tn", 
# "tot_pop_valid", "mcc", "row_rank_p", 
# "row_rank_mcc", "top_go", 
# "top_kegg","top_reac"
