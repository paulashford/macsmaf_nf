#!/usr/bin/env Rscript
# add_metrics.r
# Add classification metrics and MCC to g:Profiler enrichment results
# 17 01 2025
# Called from main.nf; implemented in add_tp_tn_fp_fn_mcc

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Required arguments: input_file output_file")
}

input_file <- args[1]
output_file <- args[2]

# Source the functions
source(file.path(Sys.getenv("NXF_SCRIPT_DIR"), "gprofiler_enrichment_functions.r"))

# Read the input data
enrichment_data <- readRDS(input_file)

# Add the metrics
enrichment_with_metrics <- add_tp_tn_fp_fn_mcc(enrichment_data)

# Save the results
saveRDS(enrichment_with_metrics, file = output_file) 