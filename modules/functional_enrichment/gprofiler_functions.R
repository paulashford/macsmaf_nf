# experiments/e019/gprofiler-multi/gprofiler_functions.R
# Programmatic use of g:Profiler as the full module queries are causing time-outs on web-page
# parsing / enrichment functions
# 05-24 05 2023

# function to concat path and read files to tibbles
get_tibble <- function( root_dir, file_name, sep = "\t" ){
    require( readr )
    return( read_delim( file.path( root_dir, file_name ), delim = sep ) )
}
write_tibble <- function( tib, out_dir, file_name, sep = "\t" ){
    require( readr )
    write_delim( tib, file = file.path( out_dir, file_name ), quote = "none", delim = sep )
}

# Read module file and convert to named list format appropriate for using directly in gost
parse_module_file <- function( base_dir, filename, module_prefix="module_" ){
	require(readr)
	require(tidyr)
	require(dplyr)
	require(purrr)

	# read Monet module file into characterlines as character vector (don"t want df/dt/tibble yet as unequal num genes per module - ie not relational/"tidy")
	modules <- read_lines( file.path(base_dir, filename) )	

	# conv vec into tibble...
	modules_tib <- tibble( modules )
	
	# Use tidyr command for splitting cells
	# As only 3 cols specified, the gene list is left alone in single col with  extra="merge"
	modules_tib <- modules_tib %>%
		separate( col="modules", into=c("mod", "junk", "genes"), sep="\t", extra="merge" ) %>%
		# junk col not needed
		select( -junk )

	# Add string prefix to module number and reorder cols
	modules_tib <- modules_tib %>%
    	mutate( module=map_chr( mod, function(x) paste0(module_prefix, x) ) ) %>%
    	select( -mod ) %>%
    	relocate( genes, .after = last_col() )

	# Convert gene list to space sep for readability
	modules_tib <- modules_tib %>%
    	rowwise() %>%
    	mutate( gl = paste0( gsub( "\t", " ", genes ) ) ) %>%
    	select( -genes )

	# pivot the genes in gl col to rows (tried to do direct col action eval into list but no luck...)
	modlong <- modules_tib %>%
    	separate_rows( gl, sep=" " ) %>%	
    	group_by( module )

	# Nest genelists within module groups
	modnest <- modlong %>%
    	nest( genes = gl)

	# for each row move nested tibbles into lists
	modlist <- modnest %>%
    	rowwise() %>%
    	mutate( genelist = list( tibble::deframe( genes ) ) )

	# get named list with names of each genelist being the module and return
	return( setNames( modlist$genelist, modlist$module ) )

}

# gprofiler gost wrapper includes boilerplate params for our purposes
# Note multi_query = TRUE, doesn"t quite do what expected as it will summarise all modules into top terms [see https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html]
# 	whereas we want enriched terms grouped by modules - provided pass named list of modules into query this is what is done in this wrapper
run_gprofiler_enrichment <- function(query, 
										sig_level=0.01, 
										sources=c( "GO:MF","GO:CC","GO:BP","KEGG","REAC","HPA","CORUM","HP","WP" ),
										exclude_go_iea=TRUE,
										multq=FALSE
										){
	require(gprofiler2)
	print( paste0( "Running g:profiler for ", length(query), " modules..." ) )

	gostres <- gost( query = query,
                	organism = "hsapiens", 
                	ordered_query = FALSE, 
                	multi_query = multq, 
                	significant = TRUE, 
                	exclude_iea = exclude_go_iea, 
                	measure_underrepresentation = FALSE, 
                	evcodes = FALSE, 
                	user_threshold = sig_level, 
                	correction_method = "g_SCS", 
                	domain_scope = "annotated", 
                	custom_bg = NULL, 
                	numeric_ns = "", 
                	sources = sources,
                	as_short_link = FALSE )

	return(gostres)

}

# Some post-processing of gprofiler results data.frame (gp_enrich$result) columns & filtering
# max_term: maximum GO term size for those that are very generic as fraction of effective domain size
# Previously (see e017/proccess_gprofiler.R) using web-version multi-module-run, resultant file was very wide and needed pivoting appropriately
# Programmatic runs of gp return a long-form (sane) table which just needs module number/experiment splitting out 
# Here is format used in e017:
#  head -n 2 long_form_v2_gProfiler_hsapiens_M1_09-02-2023_15-17-59__multiquery_simple_HDR.csv_M1.csv
# "","source","term_name","term_id","term_size","effective_domain_size","experiment_info","func_module_number","adjusted_p_value","query_size","intersection_size"
# "1","GO:MF","NADH dehydrogenase (ubiquinone) activity","GO:0008137",41,61806,"e017_M1",9,2.09101860194375e-130,59,41
# 
post_process_gprofiler_enrichment <- function( gpresults, max_term = 0.05 ){
	require(tidyr)

	gpr <- gpresults %>%
		# column naming
		rename( experiment_info = query ) %>%
		# experiment/module number fields
		separate( experiment_info, sep = "_mod_", into = c("experiment_info", "func_module_number") )
	
	# Max term size grouped by source (GO:BP, REAC etc)
	gpr <- gpr %>%
		group_by(source) %>%
		filter( term_size < (max_term * effective_domain_size) ) %>%
		# sort
		arrange(func_module_number, p_value)	
	
	# add rankings
	gpr <- gpr %>%
	   group_by(source, func_module_number) %>%
    	mutate( perc_rank = percent_rank(p_value) ) %>%
    	mutate( cumelative_dist = cume_dist(p_value) )

	return(gpr)

}

#  Add TP, FP etc confusion matrix and Matthew's Correlation Coefficient (MCC) to g:profiler output
# MCC 'generally considered one of best measures (of combining confusion matrix TP, FP etc into single measure)' [https://en.wikipedia.org/wiki/Phi_coefficient #cite_note-Powers2011-10]
# Powers, David M. W. (10 October 2020). "Evaluation: from precision, recall and F-measure to ROC, informedness, markedness and correlation". arXiv:2010.16061 [cs.LG].
# Also  from : https://en.wikipedia.org/wiki/Phi_coefficient#cite_note-Powers2011-10
#       "As explained by Davide Chicco in his paper "Ten quick tips for machine learning in computational biology" [12] (BioData Mining, 2017) and "The advantages of the Matthews correlation coefficient (MCC) over F1 score and accuracy in binary classification evaluation" [35] (BMC Genomics, 2020), 
#       "the Matthews correlation coefficient is more informative than F1 score and accuracy in evaluating binary classification problems,
#       " because it takes into account the balance ratios of the four confusion matrix categories (true positives, true negatives, false positives, false negatives).[12][35]"
# 	Example use: https://www.statology.org/matthews-correlation-coefficient-in-r/
add_tp_tn_fp_fn_MC	<- function( gpresults, inc_mcc=TRUE ){
	require(tidyr)

	gpresults <- gpresults %>%
		ungroup() %>%
		# P (all positives)
		mutate( P = term_size ) %>%
		# N (all negatives) [using total population = P + N,  hence N = total population - P]
		mutate( N = effective_domain_size - P ) %>%
		# TP
		mutate( tp = intersection_size ) %>%
		# FP
		mutate( fp = query_size - intersection_size ) %>%
		# FN
		mutate( fn = term_size - intersection_size ) %>%
		# TN
		mutate( tn = N - fp ) %>%
		# validation (should equal effective domain size)
		mutate( tot_pop_valid = tp + fp + fn + tn )
		
	# if ( inc_mcc == TRUE ){
	# 	require(mltools)
	# 	require(data.table)
	# 	# confusion matrix
	
	}
	calc_mcc <- function( tp, fn, fp, tn ){ 
		require(mltools)
		cm <- matrix( c( tp, fn, fp, tn ), nrow=2 ) 
		matthews_c	<- mcc( confusionM = cm  )
		return( as.numeric(matthews_c ) )
	}

	# for nested tibble (don't ask!)
	# calc_mcc_from_nest <- function( tib ){ 
	# 	require(mltools)
	# 	# force into useable format
	# 	confl <- as.list(tib$confusion)
		
	# 	cm <- matrix( c( confl$tp, confl$fn, confl$fp, confl$tn ), nrow=2 ) 
	# 	matthews_c	<- as.numeric(mcc( confusionM = cm  ))
	# 	assign('matthews_c',matthews_c,envir=.GlobalEnv)
	# }

		
	# 	gpresults <- gpresults[ , MCC := .( mcc( confusionM =  matrix( c( tp, fn, fp, tn ), nrow = 2 ) ) ) ]

	# 	<- gpresults %>%
	# 		rowwise() %>%
	# 		mutate( MCC = ( confusionM =  matrix( c( gpresults$tp, gpresults$fn, gpresults$fp, gpresults$tn ), nrow = 2 ) ) )
	# }
 
# 	return(gpresults)
# }

# simply return module num and concatenated GO terms for top ranked
# use for input to GOAtools for slims
get_top_go_terms_by_module <- function( gpr, min_perc_rank=0.25 ){
	require(tidyr)
	require(dplyr)
	
	go_bp_mod <- gpr %>%
    	ungroup() %>%
    	filter( source == "GO:BP" ) %>%
    	filter( perc_rank <= min_perc_rank ) %>% 
    	select( func_module_number, term_id ) %>%
    	group_by( func_module_number ) %>%
    	summarise( lst = paste0( term_id, collapse = ";" ), .groups = "keep" ) 
    
	return(go_bp_mod)
}

#  Pivot the g:Profiler enrichment table to long form 
#  From e017 as func: proc_enrich() - note this isn't required in command line gprofiler runs (eg e019)
#  These are the 3 names_to to capture (from 15/02/23 update; previously only using adjusted_p_value)
#   	adjusted_p_value__e017_R1_mod_2 = col_double(),
#   	query_size__e017_R1_mod_2 = col_double(),
#   	intersection_size__e017_R1_mod_2 = col_double(),
pivot_gp_long <- function( root_dir, enrich_file, pfilt = 0.01, split_exp_module_field = FALSE ){
	require(readr)
	require(tidyr)
	require(dplyr)

    # data frame R1 or M1
    df_mult <- read_delim( file.path( root_dir, enrich_file ), delim = "," )

    # Pivot:  "adjusted_p_value", "query_size", "intersection_size" for each GO/Reac etc term's gene-set in the rows
    df_long <- df_mult %>%
                select( c( "source", "term_name", "term_id", "term_size", "effective_domain_size" ) | starts_with( "adjusted_p_value" ) | starts_with( "query_size" ) | starts_with( "intersection_size" ) ) %>%
                pivot_longer( !c( "source", "term_name", "term_id", "term_size", "effective_domain_size" ), 
                    names_to = c( ".value", "module_info" ),
                    names_sep = "__",
                    values_drop_na = FALSE ) %>%
				# filter by p-value
				filter( adjusted_p_value < pfilt )
    
    # Split exper info into 2 fields (expID, module)
	if ( split_exp_module_field ){
		df_long <- df_long %>%
        	separate( 	module_info, 	
            		    sep = "_mod_", 
                		into = c( "experiment_info", "func_module_number" ), 
                		remove =  TRUE, 
                		convert = TRUE, 
                		extra = "warn", 
                		fill = "warn" )
	}
    

    # df_long <- arrange( df_long, func_module_number, adjusted_p_value )
    # If want to print all rows...
    # print.data.frame(df_long)
    return(df_long)
}

# from /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e017/process_modules.R
proc.modules <- function( root_dir, mod_file, type = 'none' ){
	print('use process_modules.R')
}
