# experiments/e019/process_GOslims.R
#  30/05/2023
# Use long form of gp enrichments for FULL module set to create input to owltools 
# Parse results from owltools map2slim

# ** See d014 and d011 **
# Note that for full modules here (e019) gProfiler runs were scripted, not web-based mult-queries
# Data already long in saved RDA files, but minor column name renaming should suffice

library(data.table)
# library(readr)
# library(dplyr)
library(tidyverse)

# file root
root_dir <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019'
# out directory
# out_dir <- "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v2.1/"
out_dir <- "/Users/ash/Dropbox/bioinf/MACSMAF/datasets/d014/v3/"

# parsing functions for gprofiler
source( file.path( root_dir, 'gprofiler_functions.R' ) )

# dict with all enrichment file paths and helper fns
source( file.path( root_dir, 'files.R' ) )

# grprofiler runs
# monet_type <- 'M1'; net_source <- 'cpdb'
# monet_type <- 'M1'; net_source <- 'humanbase'
# monet_type <- 'M1'; net_source <- 'string'
# monet_type <- 'R1'; net_source <- 'cpdb'
# monet_type <- 'R1'; net_source <- 'humanbase'
# monet_type <- 'R1'; net_source <- 'string'
monet_type <- 'K1'; net_source <- 'cpdb'

# load RDA of previously run gprofiler enrichment (ie Dropbox/bioinf/MACSMAF/experiments/e019/run_gprofiler.R)
# see files.R
gpr     <- load_gp_enrichments( rda_path = gp_run_dir, net_source = net_source, monet_type = monet_type )

# post-process (colnames/split/filter/sort/rank)
gpr2 <-  post_process_gprofiler_enrichment( gpr, max_term = 0.05 )

# Get the concatenated top GO terms per module for slimming
min_perc_rank = 0.25
go_for_slims <- get_top_go_terms_by_module( gpr2, min_perc_rank )
write_tibble( go_for_slims, out_dir, paste0( "go_bp_mod_rank_", min_perc_rank, "_", monet_type, "_", net_source, ".tsv" ), sep = "\t" )

