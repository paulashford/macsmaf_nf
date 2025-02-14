# GO slim pre-processing and other GO slim funcs 
# 13 01 2025

# full data / owltools info refer to:  macsr_nf/macsr_nf_dev/subworkflows/go_slim/go_slim_README.txt

# from macsr_nf/macsr_nf_dev/subworkflows/go_slim/go_slim_README.txt
# Uses GOSlim annotation pre-downloaded from [...] to simplify GO terms found from func. enrichment of a
# set of PPI network modules, i.e., 
#     distinct numbered sets of genes from modules (or elsewhere!) -> 
#     GO:BP enrichment per module number ->
#     GO-slim'd enrichment per module number (this subworkflow)

# python3 scripts/map_to_slim.py --association_file=${data_dir}/${assoc_file}.tsv ${go_dir}/go-basic.obo ${go_dir}/goslim_pir.obo > ${data_dir}/${out_file}

# simply return module num and concatenated GO terms for top ranked
# use for input to GOAtools for slims
get_top_go_terms_by_module <- function( gpr, min_perc_rank = 0.25, go_domain = "GO:BP" ){
	require(tidyr)
	require(dplyr)
	
	go_bp_mod <- gpr %>%
    	ungroup() %>%
    	filter( source == {{ go_domain }} ) %>%
    	filter( perc_rank <= min_perc_rank ) %>% 
    	select( func_module_number, term_id ) %>%
    	group_by( func_module_number ) %>%
    	summarise( lst = paste0( term_id, collapse = ";" ), .groups = "keep" ) 
    
	return(go_bp_mod)
}

# get_top_go_terms_by_module <- function( gpr, min_perc_rank=0.25 ){
# 	require(tidyr)
# 	require(dplyr)
	
# 	go_bp_mod <- gpr %>%
#     	ungroup() %>%
#     	filter( source == "GO:BP" ) %>%
#     	filter( perc_rank <= min_perc_rank ) %>% 
#     	select( func_module_number, term_id ) %>%
#     	group_by( func_module_number ) %>%
#     	summarise( lst = paste0( term_id, collapse = ";" ), .groups = "keep" ) 
    
# 	return(go_bp_mod)
# }
