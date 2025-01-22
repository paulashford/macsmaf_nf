#!/usr/bin/env Rscript
# post_process_gp.r
# Post-process g:Profiler enrichment results

suppressPackageStartupMessages({
    library(tidyverse)
})

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop(
        "Five arguments must be supplied:\n",
        "1. gprofiler_results: g:Profiler results file (RDS)\n",
        "2. method: network method (e.g. K1)\n",
        "3. db: network database (e.g. cpdb)\n",
        "4. cutoff: network cutoff value\n",
        "5. out_file: output filename\n",
        call. = FALSE
    )
}

# Assign arguments to named variables
gprofiler_results_file <- args[1]
method <- args[2]
db <- args[3]
cutoff <- args[4]
out_file <- args[5]

# Get script directory using Nextflow's environment variable
script_dir <- Sys.getenv("NXF_SCRIPT_DIR", unset = NA)
if (is.na(script_dir)) {
    stop("NXF_SCRIPT_DIR environment variable not set")
}

# Source functions file
source(file.path(script_dir, 'subworkflows', 'mod_func_enrich', 'gprofiler_enrichment_functions.r'))

# Load results
if (file.exists(gprofiler_results_file)) {
    cat("Loading RDS file:", gprofiler_results_file, "\n")
    
    tryCatch({
        gp_results <- readRDS(gprofiler_results_file)
        cat("Successfully loaded g:Profiler results\n")
    }, error = function(e) {
        stop("Failed to read RDS file: ", e$message)
    })
} else {
    stop("Input RDS file does not exist: ", gprofiler_results_file)
}

# Post-process results
processed_results <- post_process_gprofiler_enrichment(gp_results$result)

# Add method, db and cutoff info
processed_results <- processed_results %>%
    mutate(
        network_method = method,
        network_db = db,
        network_cutoff = as.numeric(cutoff)
    ) %>%
    relocate(network_method, network_db, network_cutoff, .before = everything())

# Save results
saveRDS(processed_results, file = out_file) 