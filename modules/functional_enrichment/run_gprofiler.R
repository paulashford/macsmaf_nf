# experiments/e019/gprofiler-multi/run_gprofiler.R
# Programmatic use of g:Profiler as the full module queries are causing time-outs on web-page
# 05-24 05 2023

# https://www.tidyverse.org/
# install.packages("tidyverse")
# these are in gprofiler_functions.R as requires():
# library(gprofiler2)
# library(readr)
# library(tidyr)
# library(dplyr)
# library(purrr)

# Note there is also a data.table backend to tidyverse
# https://dtplyr.tidyverse.org/
# library(data.table)
# library(dtplyr)
#  as.data.table(modlist) 
# etc

# parsing functions for gprofiler
source('gprofiler_functions.R')

# where are module files?
base_dir <- "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/modules-files/"
# output dir
out_dir <- "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/gprofiler_enrichments"

# module file names
# mprefix <- "e019_M1_mod_"
# filename <- "modules-M1-cpdb-coval-0.0.dat" # run / 780 modules
# filename <- "modules-M1-humanbase-coval-0.4.dat" # run /  647 modules
# filename <- "modules-M1-string-coval-0.2.dat" # run / 943 modules

# mprefix <- "e019_R1_mod_"
# filename <- "modules-R1-cpdb-coval-0.4.dat" # run / 367 modules
# filename <- "modules-R1-humanbase-coval-0.5.dat" # run 265 modules
# filename <- "modules-R1-string-coval-0.4.dat" # run 373 modules

# mprefix <- "e019_K1_mod_"

# enrichment sources
sources = c( "GO:BP","KEGG","REAC" )

# parse the modules
modqry <- parse_module_file( base_dir, filename, module_prefix = mprefix )
# run g:profiler
gp_enrich <- run_gprofiler_enrichment( modqry, 
										sig_level = 0.01,
										sources=sources, 
										exclude_go_iea = TRUE, 
										multq = FALSE)

# if want to inspect...
View( gp_enrich$result )

# save results
outfile <- paste0( "gprofiler_enrichments_", mprefix, filename, "_", paste0( sources,collapse="-" ), ".rda" )
outfile <- gsub( ":", "", outfile )
save( gp_enrich, file=file.path( out_dir, outfile ) )

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

	return( gostres )
	

}