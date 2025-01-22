#!/usr/bin/env Rscript
# add_metrics.r
# Add classification metrics and MCC to g:Profiler enrichment results
# 17 01 2025
# Called from main.nf; implemented in add_tp_tn_fp_fn_mcc

# Load required libraries
library(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Required arguments: input_file output_rds output_tsv")
}

input_file <- args[1]
output_rds <- args[2]
output_tsv <- args[3]

# Source functions file
source(file.path(Sys.getenv("NXF_SCRIPT_DIR"), "subworkflows", "mod_func_enrich", "gprofiler_enrichment_functions.r"))

# Read the input data
enrichment_data <- readRDS(input_file)

# Add the metrics
enrichment_with_metrics <- add_tp_tn_fp_fn_mcc(enrichment_data)

# Save the results in both formats
saveRDS(enrichment_with_metrics, file = output_rds)
write_tsv(enrichment_with_metrics, file = output_tsv) 