# files.R
# 15 06 2023
# enrichement files, slims files and helper functions

library(collections)   # has dict function

# g:profiler RDA file dictionary for each Monet type (for simplicity/laziness!)
gp_run_dir 		<- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/gprofiler_enrichments'

# GOslims PIR
go_slim_dir 	<- '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v3'

# modules dir
modules_dir		<- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/modules-files'

# final output tables dir
final_out_dir 	<- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/functional_enrichment_tables'

# Module files
# M1
dmod_m1	<- dict( items= c(	'modules-M1-cpdb-coval-0.0.dat',
							'modules-M1-humanbase-coval-0.4.dat',
							'modules-M1-string-coval-0.2.dat'	),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)
# R1
dmod_r1	<-	dict( items= c(	'modules-R1-cpdb-coval-0.4.dat',
							'modules-R1-humanbase-coval-0.5.dat',
							'modules-R1-string-coval-0.4.dat'	),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)


# RDAs from eg Dropbox/bioinf/MACSMAF/experiments/e019/run_gprofiler.R
# M1
dgp_m1	<- dict( items= c(	'gprofiler_enrichments_e019_M1_mod_modules-M1-cpdb-coval-0.0.dat_GOBP-KEGG-REAC.rda',
							'gprofiler_enrichments_e019_M1_mod_modules-M1-humanbase-coval-0.4.dat_GOBP-KEGG-REAC.rda',
							'gprofiler_enrichments_e019_M1_mod_modules-M1-string-coval-0.2.dat_GOBP-KEGG-REAC.rda' 
							),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)
# R1
dgp_r1	<-	dict( items= c(	'gprofiler_enrichments_e019_R1_mod_modules-R1-cpdb-coval-0.4.dat_GOBP-KEGG-REAC.rda',
							'gprofiler_enrichments_e019_R1_mod_modules-R1-humanbase-coval-0.5.dat_GOBP-KEGG-REAC.rda',
							'gprofiler_enrichments_e019_R1_mod_modules-R1-string-coval-0.4.dat_GOBP-KEGG-REAC.rda'
							),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)

# GO slims from go_slim_dir
# M1
dgoslim_m1	<- dict( items= c(	'go_bp_mod_rank_0.25_M1_cpdb_goslim_pir.txt',
								'go_bp_mod_rank_0.25_M1_humanbase_goslim_pir.txt',
								'go_bp_mod_rank_0.25_M1_string_goslim_pir.txt'
							),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)
# R1
dgoslim_r1	<-	dict( items= c(	'go_bp_mod_rank_0.25_R1_cpdb_goslim_pir.txt',
								'go_bp_mod_rank_0.25_R1_humanbase_goslim_pir.txt',
								'go_bp_mod_rank_0.25_R1_string_goslim_pir.txt'
							),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)

# Final output tables
dfinal_m1	<- dict( items= c(	'func_enrichments_M1_cpdb.tsv',
								'func_enrichments_M1_humanbase.tsv',
								'func_enrichments_M1_string.tsv'
							),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)
# R1
dfinal_r1	<-	dict( items= c(	'func_enrichments_R1_cpdb.tsv',
								'func_enrichments_R1_humanbase.tsv',
								'func_enrichments_R1_string.tsv'
							),
				 keys = c(	'cpdb',
							'humanbase',
							'string'
							)
				)

# ** Final output tables - RANK1 aggregate SUMMARY**
dfinal_m1_r1agg	<- dict( items= c(	'func_rank1_agg_summary_M1_cpdb.tsv',
									'func_rank1_agg_summary_M1_humanbase.tsv',
									'func_rank1_agg_summary_M1_string.tsv'
							),
				 		keys = c(	'cpdb',
									'humanbase',
									'string'
							)
						)
# R1
dfinal_r1_r1agg	<-	dict( items= c(	'func_rank1_agg_summary_R1_cpdb.tsv',
									'func_rank1_agg_summary_R1_humanbase.tsv',
									'func_rank1_agg_summary_R1_string.tsv'
							),
						keys = c(	'cpdb',
									'humanbase',
									'string'
							)
				)

# ** Final output tables - RANK1 aggregate DETAIL**
dfinal_m1_r1agg_detail	<- dict( items= c(	'func_rank1_agg_detail_M1_cpdb.tsv',
									'func_rank1_agg_detail_M1_humanbase.tsv',
									'func_rank1_agg_detail_M1_string.tsv'
							),
				 		keys = c(	'cpdb',
									'humanbase',
									'string'
							)
						)
# R1
dfinal_r1_r1agg_detail	<-	dict( items= c(	'func_rank1_agg_detail_R1_cpdb.tsv',
									'func_rank1_agg_detail_R1_humanbase.tsv',
									'func_rank1_agg_detail_R1_string.tsv'
							),
						keys = c(	'cpdb',
									'humanbase',
									'string'
							)
				)




# Load relevant enrichment R1/M1 cpdb etc data
load_gp_enrichments	<- function( rda_path = '/', net_source = 'network', monet_type = 'X' ){
	require(tidyverse)
	if ( monet_type == 'M1' ){
		dgp <- dgp_m1
	}else if ( monet_type == 'R1' ) {
		dgp <- dgp_r1
	}else if ( monet_type == 'K1' ) {
		dgp <- dgp_k1
	} else{
		return( "Monet types M1, R1, (K1)" )
	}
	enrich_file <- dgp$get( net_source )
	load( file.path( rda_path, enrich_file ) )
 
	# split query field into exp / module number
 	gp_enrich	<- gp_enrich$result %>%
    	separate( query, 
                sep = "_mod_", 
                into = c( "experiment_info", "func_module_number" ), 
                remove =  TRUE, 
                convert = TRUE, 
                extra = "warn", 
                fill = "warn" )


	return( as_tibble(gp_enrich ) )
}

# Load relevant GO SLIMS R1/M1 cpdb etc data
load_goslims	<- function( slim_dir = '/', net_source = 'network', monet_type = 'X' ){
	require(tidyverse)
	if ( monet_type == 'M1' ){
		dgoslim <- dgoslim_m1
	}else if ( monet_type == 'R1' ) {
		dgoslim <- dgoslim_r1
	}else if ( monet_type == 'K1' ) {
		dgoslim <- dgoslim_k1
	} else{
		return( "Monet types M1, R1, (K1)" )
	}
	# Read slims
	slim_file <- dgoslim$get( net_source )
	df_goslim <- read_delim( file.path( slim_dir, slim_file ), delim = '\t', skip = 4 ) 
	names(df_goslim) <- c("func_module_number", "goslim_pir")
	return( df_goslim)
	# return(slim_file)
}

# Load final RANK1 AGGREGATE results tables for doing post-analysis such as plots
load_final_rank1	<- function( final_dir = '/', net_source = 'network', monet_type = 'X', summary = TRUE ){
	require(tidyverse)
	
	if ( monet_type == 'M1' ){
		if (summary == TRUE) 	{ dfinal <- dfinal_m1_r1agg }
		if (summary == FALSE) 	{ dfinal <- dfinal_m1_r1agg_detail }
	}else if ( monet_type == 'R1' ) {
		if (summary == TRUE) 	{ dfinal <- dfinal_r1_r1agg }
		if (summary == FALSE) 	{ dfinal <- dfinal_r1_r1agg_detail }
	}else if ( monet_type == 'K1' ) {
		if (summary == TRUE) 	{ dfinal <- dfinal_k1_r1agg }
		if (summary == FALSE) 	{ dfinal <- dfinal_k1_r1agg_detail }	
		return( "Monet types M1, R1, (K1)" )
	}
	# Read slims
	final_file <- dfinal$get( net_source )
	df_final <- get_tibble( root_dir = final_dir, file_name = final_file , sep = '\t' ) 
	return( df_final )
}

# Load final results tables for doing post-analysis such as plots
load_final	<- function( final_dir = '/', net_source = 'network', monet_type = 'X' ){
	require(tidyverse)
	if ( monet_type == 'M1' ){
		dfinal <- dfinal_m1
	}else if ( monet_type == 'R1' ) {
		dfinal <- dfinal_r1
	}else if ( monet_type == 'K1' ) {
		dfinal <- dfinal_k1
	} else{
		return( "Monet types M1, R1, (K1)" )
	}
	# Read slims
	final_file <- dfinal$get( net_source )
	df_final <- get_tibble( root_dir = final_dir, file_name = final_file , sep = '\t' ) 
	return( df_final )
}

# long form modules for annotation or putting a gene names field back into enrichments file
# from /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e017/process_modules.R
proc_modules <- function( modules_dir ="/", net_source = 'network', monet_type = 'X' ){
	require(tidyverse)
	if ( monet_type == 'M1' ){
		dmod <- dmod_m1
	}else if ( monet_type == 'R1' ) {
		dmod <- dmod_r1
	}else if ( monet_type == 'K1' ) {
		dmod <- dmod_k1
	} else{
		return( "Monet types M1, R1, (K1)" )
	}

	# module file 
	module_file <- dmod$get( net_source )
	
    # Note: this no longer applies here -> [modue file needs 'dummy' header row C1, C2, C3 with at least (or more) cols as max genes per module (workaround!)]
   	# Process modules to display format - previous read.table problematic rows (would newline on odd cols and mess up table!)
	# raw_mod <- read.table( file.path( modules_dir, module_file ),  fill = TRUE, header = FALSE, quote = "" )
    # df_mod = as_tibble( raw_mod )
	# remove dummy 2nd col which is always 1
    # df_mod = mutate( df_mod, V2 = NULL ) 
	# New code adapted from: https://stackoverflow.com/questions/74942222/how-do-you-make-readr-read-all-the-columns-in-a-csv-file-when-the-file-has-a-var
	df_mod <- dplyr::tibble( rawdata = readr::read_lines( file.path( modules_dir, module_file ) ) )  %>%
    	# insert line identifier
    	dplyr::mutate( fileline = dplyr::row_number(),  
                 # remove " which in .csv are text delimiters and are processed correctly when improting as .csv but not as text/lines 
                  rawdata = stringr::str_remove_all( rawdata, pattern = '\"' ) ) %>%
		# split acording to the delimiter (a comma in this case)
		tidyr::separate_rows( rawdata, sep = '\t' ) %>% 
		# build groupings according to line identifier
		dplyr::group_by( fileline ) %>% 
		# insert a line identifier per group (imported line)
		dplyr::mutate( colnum = dplyr::row_number() ) %>% 
		# release grouping to avoid unwanted behaviour down stream
		dplyr::ungroup() %>% 
		# make the table wide
		tidyr::pivot_wider( names_from = colnum, values_from = rawdata ) %>%
		# you might want to rename the columns at this stage and/or drop the "fileline" column
		rename_with( ~ paste0("col_", .x, recycle0 = TRUE) )%>%
		rename( module = `col_1`) %>%
		select( -c( col_fileline, col_2 ) )
		
    # pivot long form with module prefixes
    df_mod_long <- df_mod %>%
		pivot_longer( !c( "module" ), names_to = "dummy", values_to = "gene", names_prefix = "C" , values_drop_na = TRUE ) %>%
		# pivot_longer( !c( "V1" ), names_to = "dummy", values_to = "gene", names_prefix = "C" , values_drop_na = TRUE ) %>%
        # pivot_longer( !c( "module" ), names_to = "dummy", values_to = "gene", names_prefix = "C" , values_drop_na = TRUE ) %>%
        filter( !gene == "" ) %>%
		select( -dummy ) %>%
        mutate( monet_type =  monet_type ) 
		# mutate( V1 = paste0( monet_type, "_", V1 ) )
	# sometimes base R better/easier	
    # names( df_mod_long ) <- c( 'module', 'gene', 'monet_type' )
	df_mod_long$module <- as.numeric(df_mod_long$module)
    
	return(df_mod_long)
}