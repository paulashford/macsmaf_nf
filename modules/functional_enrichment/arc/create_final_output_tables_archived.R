# experiments/e019/create_final_output_tables_archived.R
# This code at end of script experiments/e019/create_final_output_tables.R
# was not required but may still have useful info about process to get to the final tables
# and summary counts etc...  
# 05 07 2023

# <...>


# testm1c <- get_tibble(out_dir, 'func_enrichments_M1_cpdb.tsv')
# testm1h <- get_tibble(out_dir, 'func_enrichments_M1_humanbase.tsv')
# testm1s <- get_tibble(out_dir, 'func_enrichments_M1_string.tsv')

# testr1c <- get_tibble(out_dir, 'func_enrichments_R1_cpdb.tsv')
# testr1h <- get_tibble(out_dir, 'func_enrichments_R1_humanbase.tsv')
# testr1s <- get_tibble(out_dir, 'func_enrichments_R1_string.tsv')
# colnames(etc)

# -----------------------------------------------------------------------------
	

# test <- comb_summary2 %>%
# 	rowwise() %>%
# 	map_dbl( list(.tp, .fn, .fp, .tn), ~ calc_mcc(.tp, .fn, .fp, .tn) )
# 		# mutate( mccval = map_dbl(ac))

# test <- comb_summary2 %>%
# 	rowwise() %>%
# 	mutate( matthews_corr = calc_mcc( c_across( all_of( c('tp', 'fn', 'fp', 'tn') ) )) )
	
# # transmute(rowwise(comb_summary2), total = sum(c_across( all_of( c('fp', 'tp') ) ) ))

# dt <- data.frame(comb_summary2)
# test <- lapply( dt[ c("tp", "fn", "fp", "tn" )], cm )

	# distinct(-genes)
	# v1 test summarise( module_mod, module, source, term_id, term_name, goslim_pir_025, perc_rank, adjusted_p_value, perc_rank, genes = paste0( unique(gene), collapse = ";" ),term_size, effective_domain_size) %>%
    # summarise( module_mod, module, source, term_id, term_name, goslim_pir_025, adjusted_p_value,  perc_rank, genes = paste0( unique(gene), collapse = ";" ), query_size, term_size, intersection_size,  effective_domain_size ) %>%
    # v1 test distinct(  module_mod, module, source, term_id, term_name, goslim_pir_025, perc_rank, adjusted_p_value, perc_rank, genes ,term_size, effective_domain_size, .keep_all = TRUE )
    # distinct( module_mod, module, source, term_id, term_name, goslim_pir_025, adjusted_p_value,  perc_rank, genes, query_size, term_size, intersection_size,  effective_domain_size, .keep_all = TRUE )


##########################################################
##  NOTE: 06 06 23
## COVERAGE of term in module in addition to p-value ranking
## not sure exactly what calc this is, but this is a placeholder reminder :)
###############################################################
# see the notes and refs about MCC here:
# experiments/e019/gprofiler_functions.add_tp_tn_fp_fn_MC()

library( mltools )
mcc(preds = NULL, actuals = NULL, TP = NULL, FP = NULL, TN = NULL,
      FN = NULL, confusionM = NULL)




# g:Profiler parsed files in long-form (parse with: experiments/e019/proccess_gprofiler.R)
# v1
# gp_r1   <- 'gprofiler_multi/long_form_gProfiler_hsapiens_R1_13-12-2022_13-18-07__multiquery_simple_HDR.csv_R1.csv'
# gp_m1   <- 'gprofiler_multi/long_form_gProfiler_hsapiens_M1_09-02-2023_15-17-59__multiquery_simple_HDR.csv_M1.csv'
# NOTE was done with ANNOTATED GENES unlike ALL for m1/r1 - leave for now
# gp_k1   <- 'gprofiler_multi/long_form_gProfiler_hsapiens_K1_26-01-2023_18-48-10__multiquery_simple_HDR.csv_K1.csv'

# v2 ("v2" has additional columns query_size and intersection_size)
gp_r1 <- "gprofiler_multi/long_form_v2_gProfiler_hsapiens_R1_13-12-2022_13-18-07__multiquery_simple_HDR.csv_R1.csv"
gp_m1 <- "gprofiler_multi/long_form_v2_gProfiler_hsapiens_M1_09-02-2023_15-17-59__multiquery_simple_HDR.csv_M1.csv" 

# Module annotation files (basically just long form of module no/info and genes - process_modules.R)
mod_annot_r1 <- "annot/mod_annot_R1.csv"
mod_annot_m1 <- "annot/mod_annot_M1.csv"

# Which Monet method to process?
# monet_type <- "R1"
# gp <- gp_r1
# mod_annot <- mod_annot_r1

monet_type <- "M1"
gp <- gp_m1
mod_annot <- mod_annot_m1

# GO slims
# R1 original
slim_file <- '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/go_bp_mod_q025_goslim_pir.txt'
# R1 new (same)
slim_file <- "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2.1/????"
# M1
slim_file <- "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2.1/go_bp_mod_q025_goslim_pir_M1.txt"

# read g:Profiler datasets of long-form results parsed with experiments/e019/proccess_gprofiler.R
df_enrich <- get_tibble( root_dir, gp )
#enrich_long_r1  <- read_enrich_long( root_dir, gp_r1 )

# read module annot files
man <- get_tibble( root_dir, mod_annot )
#man_r1 <- read_delim( '', delim = "," )

# unique(enrich_long_r1$effective_domain_size)
# [1] 61806
# unique(enrich_long_m1$effective_domain_size)
# [1] 61806
# NA in v2 so far....unique(enrich_long_k1$effective_domain_size)
    #  [1] 16964 17194 16205 10461  8064  3385  7827  4756 10976 19943 14830
# Max term size
max_term = 0.05 *  unique( df_enrich$effective_domain_size )

# filter sources (for now)
# unique(enrich_long_r1$source) # note for M1 I didn't test "TF" or "MIRNA"
# "REAC"  "KEGG"  "GO:BP" "HP"    "GO:CC" "TF"    "WP"    "MIRNA" "GO:MF"  "CORUM" "HPA" 
df_enrich <- filter( df_enrich, source %in% c( 'GO:BP', 'REAC', 'KEGG' ) )

# Remove super-generic terms
df_slim <- filter( df_enrich, term_size < max_term )

# Get gene-based list of GO terms
# Use GOATOOLS for mapping eg /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/goa_map_slim_test.sh

# Join annot and long enrichment file  R1: [174,328 × 11] note is >800k without max_term and source filters!
comb <- inner_join( man, df_slim, by = c( "module" = "func_module_number" ) )
# sort by module and p-vals
comb <- arrange( comb, module, adjusted_p_value )
# for testing
# View(filter( comb,module < 11 & source == "GO:BP" )) 
View(filter( comb,module < 11 )) 
# Add Reactome-GO mapping (is a bit useless for immune mappings!  )
# SKIP FOR NOW
# df_reac <- read_delim( '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/cross_references/reactome2go_modified', delim = "~", skip = 1)
# comb_r1 <- left_join( comb_r1, df_reac, by = c( "term_id" = "pathway_id" ), na_matches = "never", keep = TRUE )

# Rank p-values WITHIN EACH MODULE
# Rank p-values by 
# R-help You use a data frame to create multiple columns so you can wrap
# this up into a function:
# my_quantile <- function(x, probs) {
#   tibble(x = quantile(x, probs), probs = probs)
# }
comb <- comb %>%
    group_by( module ) %>%
    mutate( perc_rank = percent_rank( adjusted_p_value ) ) %>%
    mutate( cumelative_dist = cume_dist( adjusted_p_value ) )

# Generate GOSLIM input file (see bottom of this script)
# add header to line 5 of file eg d014/v2.1/go_bp_mod_q025_goslim_pir_M1.txt:
# module	goslim_pir_025

# GO SLIM PIR FOR GO:BP q025 [/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/goa_map_slim_test.sh - see bottom of this script for dataset creation for the .sh script]
df_goslim <- read_delim( slim_file, delim = "\t", skip = 4 ) 
comb_slim <- left_join( comb, df_goslim, by = c( "module" = "module" ), na_matches = "never", keep = FALSE )

# Tidy cols and filter perc_rank 0.25 
comb_slim <- comb_slim %>%
    filter( perc_rank <= 0.25 ) %>%
    select( module_mod, module, gene, source, term_id, term_name, goslim_pir_025, adjusted_p_value,  perc_rank, query_size, term_size, intersection_size,  effective_domain_size )
# write_delim( comb_r1_slim, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/comb_r1_slim.tsv", quote = "none", delim = "\t")
# write_delim( comb_r1_slim, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2/comb_r1_slim.tsv", quote = "none", delim = "\t")
write_tibble( comb_slim, out_dir = out_dir, file_name = paste0( "comb_", monet_type, "_slim.tsv" ) )

# Combine genes into column
comb_summary <- comb_slim %>%
    group_by( module, source ) %>%
    # v1 test summarise( module_mod, module, source, term_id, term_name, goslim_pir_025, perc_rank, adjusted_p_value, perc_rank, genes = paste0( unique(gene), collapse = ";" ),term_size, effective_domain_size) %>%
    summarise( module_mod, module, source, term_id, term_name, goslim_pir_025, adjusted_p_value,  perc_rank, genes = paste0( unique(gene), collapse = ";" ), query_size, term_size, intersection_size,  effective_domain_size ) %>%
    # v1 test distinct(  module_mod, module, source, term_id, term_name, goslim_pir_025, perc_rank, adjusted_p_value, perc_rank, genes ,term_size, effective_domain_size, .keep_all = TRUE )
    distinct( module_mod, module, source, term_id, term_name, goslim_pir_025, adjusted_p_value,  perc_rank, genes, query_size, term_size, intersection_size,  effective_domain_size, .keep_all = TRUE )

# Create overall rank per module
comb_summary <- comb_summary %>%
    group_by( module ) %>%
    mutate( rank = rank( adjusted_p_value, ties.method = "average" ) ) %>%
    arrange( module, rank )
# write_delim( comb_r1_summary, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/comb_r1_summary.tsv", quote = "none", delim = "\t")
# write_delim( comb_r1_ummary, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2/comb_r1_summary.tsv", quote = "none", delim = "\t")
write_tibble( comb_summary, out_dir = out_dir, file_name = paste0( "comb_summary_ranked_by_module_", monet_type, ".tsv" ) )

# top ranked (function not working yet - do individually below!!)
# comb_top <- function( comb_summary, source_filter, rank_filter = 1, grouping = "module" ){
#     comb_summary_filt <- comb_summary %>%
#     filter( source == source_filter ) %>%    
#     group_by( grouping ) %>%
#     mutate( paste0( "rank_", source_filter ) = rank( adjusted_p_value, ties.method = "average" ) ) %>%
#     filter( paste0( "rank_", source_filter )  == 1 )
# }

# Remove GO terms except top ranked (which will have GOSLIM anyway!)
mucomb_summary_GOBP <- comb_summary %>%
    filter( source == "GO:BP" ) %>%
    group_by( module ) %>%
    mutate( rankGOBP = rank( adjusted_p_value, ties.method = "average" ) ) %>%
    filter( rankGOBP == 1 )
    # arrange( module, rank )
# write_delim( comb_r1_summary_GOBP, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/comb_r1_summary_GOBP.tsv", quote = "none", delim = "\t")
# write_delim( comb_r1_summary_GOBP, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2/comb_r1_summary_GOBPr1.tsv", quote = "none", delim = "\t")
write_tibble( comb_summary_GOBP, out_dir = out_dir, file_name = paste0( "comb_summary_GOBP_rank_1_", monet_type, ".tsv" ) )

# v2 Remove KEGG terms except top ranked 
comb_summary_KEGG <- comb_summary %>%
    filter( source == "KEGG" ) %>%
    group_by( module ) %>%
    mutate( rankKEGG = rank( adjusted_p_value, ties.method = "average" ) ) %>%
    filter( rankKEGG == 1 )
#write_delim( comb_r1_summary_KEGG, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2/comb_r1_summary_KEGGr1.tsv", quote = "none", delim = "\t")
write_tibble( comb_summary_KEGG, out_dir = out_dir, file_name = paste0( "comb_summary_KEGG_rank_1_", monet_type, ".tsv" ) )

# v2 Remove REAC terms except top ranked 
comb_summary_REAC <- comb_summary %>%
    filter( source == "REAC" ) %>%
    group_by( module ) %>%
    mutate( rankREAC = rank( adjusted_p_value, ties.method = "average" ) ) %>%
    filter( rankREAC == 1 )
# write_delim( comb_r1_summary_REAC, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2/comb_r1_summary_REACr1.tsv", quote = "none", delim = "\t")
write_tibble( comb_summary_REAC, out_dir = out_dir, file_name = paste0( "comb_summary_REAC_rank_1_", monet_type, ".tsv" ) )

# v2 combine top ranked for GOBP, KEGG and REAC
comb_summary_rank1 <- 
    bind_rows( comb_summary_GOBP, comb_summary_KEGG ) %>%
    bind_rows( comb_summary_REAC ) %>%
    arrange( module, source )
# write_delim( comb_r1_rank1, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2/comb_r1_rank1_summary.tsv", quote = "none", delim = "\t" )
write_tibble( comb_summary_rank1, out_dir = out_dir, file_name = paste0( "comb_summary_rank_1_", monet_type, ".tsv" ) )

# v1 Combine with the KEGG and REAC
# comb_r1_summary_condensed <- comb_r1_summary %>%
#     filter( source %in% c( "KEGG", "REAC" ) ) %>%
#     group_by( module ) %>%
#     mutate(rankpath = rank(adjusted_p_value, ties.method = "average" ) ) %>%
#     filter( rankpath == 1 ) %>%
#     bind_rows( comb_r1_summary_GOBPr1 ) %>%
#     arrange( module, rank )
# write_delim( comb_r1_summary_condensed, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/comb_r1_summary_condensed.tsv", quote = "none", delim = "\t")    
    

# ***************************************    
# Running GO slims GOATOOLS
# ***************************************
# original R1: /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/go_bp_mod_q025_goslim_pir.txt
# /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/goa_map_slim_test.sh

# GO GP - NO genes
# go_bp_mod <- comb %>%
go_bp_mod <- gpr2 %>%
    ungroup() %>%
    filter( source == "GO:BP" ) %>%
    filter( perc_rank <= 0.25 ) %>% 
    select( func_module_number, term_id ) %>%
    group_by( func_module_number ) %>%
    summarise( lst = paste0( term_id, collapse = ";" ), .groups = "keep" ) 
    # group terms into single list by gene
    # group_by( func_module_number )

go_bp_mod <- distinct ( go_bp_mod )
go_bp_mod <- go_bp_mod %>%
    group_by( module ) %>%
    summarise( lst = paste0( term_id, collapse = ";" ) )
    # transmute( mod_gene = paste0(module, gene, collapse = "_"), terms = lst)
# write_delim( go_bp_mod, file = "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/go_bp_mod_q025.tsv", quote = "none", delim = "\t")
write_tibble( go_bp_mod, out_dir, paste0( "go_bp_mod_q025_", monet_type,".tsv" ), sep = "\t" )


# # GO GP - with genes
# go_bp <- comb_r1 %>%
#     # filter( module < 11 & source == "GO:BP" ) %>%
#     filter( source == "GO:BP" ) %>%
#     filter( perc_rank <= 0.25 ) %>%
#     select( c( "module", "gene", "term_id", "term_name" ,"adjusted_p_value ") ) 
    
    
    # # group terms into single list by gene
    # group_by( module ) %>%
    # summarise(lst= paste0(term_id, collapse = ";"))
    # # transmute( mod_gene = paste0(module, gene, collapse = "_"), terms = lst)




go_bp_united <- go_bp %>%
    unite( col='module_gene', c('module', 'gene'), sep='_' )

# write_delim( test1, "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/test2.txt", delim='\t' )
# write_delim( go_bp, "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/go_bp.txt", delim='\t' )
write_delim( go_bp, "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/test/go_bp_mod_gene.txt", delim='\t' )







# # Write id files
# id_type = 'BP'
# id_name <- file.path( root_dir, 'id_sets', paste0( monet_label, '_', id_type, '.txt' ) )

# write.csv(
#     unique( filter( df_slim, source == paste0('GO:', id_type ) )$term_id ),
#     file = id_name,
#     row.names = FALSE, 
#     quote = FALSE )



# # Get GAF files for just the enriched GO terms
# cd /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011

# # Note - remove the column header "x" from id lists, or everything matches :)
# # grep -F -f id_sets/gP_M1_BP.txt /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf > id_sets/gP_M1_BP.gaf    
# grep -F -f id_sets/gP_R1_BP.txt /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf > id_sets/gP_R1_BP.gaf    

# # (2) I need to provide the relevant GO terms in the .gaf file ie:
#     # https://github.com/owlcollab/owltools/wiki/Map2Slim
#     # using an existing slim
#     # owltools go.obo --gaf annotations.gaf --map2slim --subset goslim_pombe --write-gaf annotations.mapped.gaf
#     # (a) So use ID list to subset the human general gaf file
# # grep -F -f idlist_R1_GOBP_GOMF_gProfiler_13-12-2022_13-18-07.txt goa_human.gaf > goa_human_R1_term_GOBP_GOMF_t1307.txt
# cd /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019
# # Initial M1 experiments/e019/id_sets/gP_M1_BP_v1.gaf had no -w flag (probably doesn't matter)
# grep -F -f id_sets/gP_M1_BP.txt /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf > id_sets/gP_M1_BP_v1.gaf
# grep -w -F -f id_sets/gP_M1_BP.txt /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf > id_sets/gP_M1_BP.gaf
# grep -w -F -f id_sets/gP_R1_BP.txt /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf > id_sets/gP_R1_BP.gaf    

# # Define GOSlim variables
# # Helful thread! https://github.com/owlcollab/owltools/issues/124
# go_set=/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/go.obo
# slim_set=/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goslims/goslim_pir
# # slim_set=/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goslims/goslim_generic
# # enriched_gaf=/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/id_sets/gP_M1_BP_v1.gaf
# # enriched_gaf=/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/id_sets/gP_M1_BP.gaf
# # enriched_gaf=/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/id_sets/gP_R1_BP.gaf
# enriched_gaf=/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf

# # owltools $go_set --gaf $enriched_gaf --map2slim --subset $slim_set --write-gaf annotations.mapped_gP_M1_BP_v1.gaf > annotations.mapped_gP_M1_BP_v1.console
# owltools $go_set --gaf $enriched_gaf --map2slim --subset $slim_set --write-gaf annotations.mapped_goa_human.gaf > annotations.mapped_goa_human.console





# # data.table form for my sanity!
# dt_slim <- data.table( df_slim )
# # num terms after removing 'super' terms
# # 11,776
# str(df_slim)
# #'data.frame':   37750 obs. of  9 variables:
# str(dt_slim)
# # Classes ‘data.table’ and 'data.frame':  37750 obs. of  9 variables:
# str(unique(dt_slim$term_id))
# # chr [1:11776] "REAC:R-HSA-5668541" "KEGG:04672" "REAC:R-HSA-5669034" ...
# # Just GO for now
# length( unique( dt_slim[source %in% c("GO:BP","GO:MF"), term_id] ) )
# # 4,316

# # GOSlim terms (GO:BP and GO:MF - need to map Reac and KEGG)
# write.csv(
#     unique( dt_slim[source %in% c( "GO:BP", "GO:MF" ), term_id] ),
#     file = '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/idlist_R1_GOBP_GOMF_gProfiler_13-12-2022_13-18-07.txt',
#     #file = '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/df_R1_term_GOBP_GOMF_ids_long_form_gProfiler_hsapiens_K1_26-01-2023_18-48-10.txt', 
#     row.names = FALSE, 
#     quote = FALSE )

# # GO slims using PIR subset http://geneontology.org/docs/go-subset-guide/
# # also goslims/goslim_generic.obo

# # Attempting to SLIM the network GO terms...
# cd /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011
# # (1) I got the parameters wrong on this one - the idfile, it seems, is for GO sim IDs from user defined slim (I think!)
# # ie # https://github.com/owlcollab/owltools/wiki/Map2Slim
#     # Example command lines:
#     #     using a custom slim from an id file:
#     # owltools go.obo --gaf annotations.gaf --map2slim --idfile slim.terms --write-gaf annotations.mapped.gaf
#     ## INCORRECT for ref:
# owltools 
#     goslims/goslim_pir.obo 
#     --gaf goa_human.gaf 
#     --map2slim 
#     --idfile df_R1_term_GOBP_GOMF_ids_long_form_gProfiler_hsapiens_K1_26-01-2023_18-48-10.txt
#     --write-gaf annotations.mapped_R1_GOBP_GOMF_term_ids_long_form_gProfiler_hsapiens_K1_26-01-2023_18-48-10.gaf
#  # owltools goslims/goslim_pir.obo --gaf goa_human.gaf --map2slim --idfile df_R1_term_GOBP_GOMF_ids_long_form_gProfiler_hsapiens_K1_26-01-2023_18-48-10.txt --write-gaf annotations.mapped_R1_GOBP_GOMF_term_ids_long_form_gProfiler_hsapiens_K1_26-01-2023_18-48-10.gaf

# # (2) I need to provide the relevant GO terms in the .gaf file ie:
#     # https://github.com/owlcollab/owltools/wiki/Map2Slim
#     # using an existing slim
#     # owltools go.obo --gaf annotations.gaf --map2slim --subset goslim_pombe --write-gaf annotations.mapped.gaf
#     # (a) So use ID list to subset the human general gaf file
#     # grep -F -f idlist_R1_GOBP_GOMF_gProfiler_13-12-2022_13-18-07.txt goa_human.gaf > goa_human_R1_term_GOBP_GOMF_t1307.txt
#     cd /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019

#     grep -F -f id_sets/gP_M1_BP.txt /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goa_human.gaf > id_sets/gP_M1_BP.gaf
#     owltools /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/go-basic.obo --gaf id_sets/test1.gaf --map2slim --subset /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/goslims/goslim_pir.obo --write-gaf annotations.mapped_test1.gaf
 

# df_test3 <- filter(
#     df_enrich, 
#     func_module_number == 1 & 
#     source == "GO:BP" & 
#     term_size < max_term )
# write.csv(df_test3$term_id , file = '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/df_test3_ids.txt', row.names = FALSE, quote = FALSE )
# # cd /Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011
# # owltools go.obo --gaf goa_human.gaf --map2slim --idfile df_test3_ids.txt --write-gaf annotations.mapped_df_test3.gaf
# # tib_test3 <- read_delim( '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/annotations.mapped_df_test3.gaf' , delim = "\t", skip = 5 )
# tib_test3 <- read_delim( '/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d011/annotations.mapped_testgen.gaf' , delim = "\t", skip = 5 )
# # Try using a GOSlim ontology this time... :)
# #owltools goslims/goslim_generic.obo --gaf goa_human.gaf --map2slim --idfile df_test3_ids.txt --write-gaf annotations.mapped_testgen.gaf





# # From 
# # some dply tests
# df_test2 <- enrich_long_r1 %>%
#     filter(source == "GO:BP") %>%
#     arrange(source, adjusted_p_value)

# # For GOslim testing  x-scrivener-item:///Users/ash/Dropbox/bioinf/MACSMAF/scriv/macsmaf_research.scriv?id=60A8403E-C3BB-4203-9282-14712EFCD513
# df_test <- filter(df_enrich, func_module_number == 1 & source == "GO:BP")
# df_test$term_id
# # lots of GOslim mappings as has generic terms... cut the generics (say >5%)
#  unique(df_test$effective_domain_size)

