# experiments/e019/proccess_modules.R
# FROM e017
# 15 06 23
# Process modules from DREAM/Monet ready for g:Profiler analysis
# NOTE: in e015 this function was (confusingly) in process_gprofiler.R 

library(tidyverse)
library(vegan)
library(data.table)

# bash process to find max cols, add header and simplify filenames
# How many columns (num. genes + 3) in each module and which is largest?
# FROM e017
# 15 06 23/ash/Dropbox/bioinf/MACSMAF/experiments/e019/monet_results_outputs
# cd /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/modules-files
# [92] awk -F'\t' '{print NF}' modules-M1-cpdb-coval-0.0.dat | sort -n
# [99] awk -F'\t' '{print NF}' modules-M1-humanbase-coval-0.4.dat | sort -n
# [99] awk -F'\t' '{print NF}' modules-M1-string-coval-0.2.dat | sort -n
# [102] awk -F'\t' '{print NF}' modules-R1-cpdb-coval-0.4.dat | sort -n
# [87] awk -F'\t' '{print NF}' modules-R1-humanbase-coval-0.5.dat | sort -n
# [94] awk -F'\t' '{print NF}' modules-R1-string-coval-0.4.dat | sort -n

# NOTE: need a 'dummy' header with at least max num of cols as genes eg
# cd /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/modules-files
# NO LONGER USED - see e019/files.R.procmodules()
        # cat ../monet_dummy_header_102_wide.txt modules-M1-cpdb-coval-0.0.dat > modules-M1-cpdb-coval-0.0_dhdr.tsv
        # cat ../monet_dummy_header_102_wide.txt modules-M1-humanbase-coval-0.4.dat > modules-M1-humanbase-coval-0.4_dhdr.tsv
        # cat ../monet_dummy_header_102_wide.txt modules-M1-string-coval-0.2.dat > modules-M1-string-coval-0.2_dhdr.tsv
        # cat ../monet_dummy_header_102_wide.txt modules-R1-cpdb-coval-0.4.dat > modules-R1-cpdb-coval-0.4_dhdr.tsv
        # cat ../monet_dummy_header_102_wide.txt modules-R1-humanbase-coval-0.5.dat > modules-R1-humanbase-coval-0.5_dhdr.tsv
        # cat ../monet_dummy_header_102_wide.txt modules-R1-string-coval-0.4.dat > modules-R1-string-coval-0.4_dhdr.tsv


# FROM e017
# 15 06 23r root of e019/
# mkdir monet_results_outputs
# mv 2022-* monet_results_outputs
# mkdir annot/ 

# files
# FROM e017
# 15 06 23 '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019'

# FROM e017
# 15 06 23et modules (see e019.sh)
mod_r1 <- 'mod_R1_csr-edge-data_simple_dHDR.tsv'
mod_m1 <- 'mod_M1_csr-edge-data_simple_dHDR.tsv'
mod_k1 <- 'mod_K1_csr-edge-data_simple_dHDR.tsv'
mod_k1_n150 <- 'mod_K1_n150_csr-edge-data_simple_dHDR.tsv'
mod_k1_n200 <- 'mod_K1_n200_csr-edge-data_simple_dHDR.tsv'
mod_k1_n250 <- 'mod_K1_n250_csr-edge-data_simple_dHDR.tsv'

# which one to use?
# mod_type <- "R1"
# mod_process <- mod_r1
# mod_type <- "M1"
# mod_process <- mod_m1
# mod_type <- "K1"
# mod_process <- mod_k1
# mod_type <- "K1_n150"
# mod_process <- mod_k1_n150
# mod_type <- "K1_n200"
# mod_process <- mod_k1_n200
mod_type <- "K1_n250"
mod_process <- mod_k1_n250


# Process this module to long form for g:Profiler
proc_mod <- proc.modules( root_dir, mod_process , type = mod_type )
write.csv( proc_mod, file.path( root_dir, "annot", paste0( "mod_annot_", mod_type, ".csv") ) )

# 15 06 23
# see /Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/files.R proc_modules()
# proc.modules <- function( root_dir, mod_file, type = 'none' ){
#     # modue file needs 'dummy' header row C1, C2, C3 with at least (or more) cols as max genes per module (workaround!)
#    # Process modules to display format
#     raw_mod <- read.table( file.path( root_dir, mod_file ),  fill = TRUE, header = TRUE, quote = "" )
#     #df_mod <- read_tsv(file.path(root_dir, mod_file), col_names = FALSE, guess_max = Inf)
#     df_mod = as_tibble( raw_mod )
#     # df_mod = mutate( df_mod, C2 = NULL ) 
#     df_mod = mutate( df_mod, dummy = NULL ) 
#     # print.data.frame( df_mod )
#     df_mod
#     # pivot long form with module prefixes
#     df_mod_long <- df_mod %>%
#         pivot_longer( !c( "module" ), names_to = "dummy", values_to = "gene", names_prefix = "C" , values_drop_na = TRUE ) %>%
#         filter( !gene == "" ) %>%
#         mutate( dummy = NULL ) %>%
#         mutate( module_mod = paste0( type, "_", module ))
#     # print.data.frame(df_mod_long)
#     return(df_mod_long)
# }
