#!/usr/bin/env Rscript
# mod_parse.r
# Module functional enrichment and annotation
# Parse a HUGO gene combined cut-off module file, e.g. parse K1_string, or K1_cpdb etc
# 14 01 2025
# Refer: 
#   macsr_nf/macsr_nf_dev/subworkflows/mod_func_enrich/mfe_README.txt
#   macsr_nf/macsr_nf_dev/subworkflows/mod_func_enrich/modferan.r

# Load required packages
suppressPackageStartupMessages({
    library(tidyverse)
})

# Get script directory using Nextflow's environment variable
script_dir <- Sys.getenv("NXF_SCRIPT_DIR", unset = NA)
if (is.na(script_dir)) {
    stop("NXF_SCRIPT_DIR environment variable not set")
}

# Source functions file
source(file.path(script_dir, 'gprofiler_enrichment_functions.r'))

# Command line Args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop(
        "Two arguments must be supplied:\n",
        "1. mod_file: a HUGO gene combined cut-off module file (e.g. for K1_string)\n",
        "2. module_prefix: prefix module naming string (e.g. 'network_modules_')",
        call. = FALSE
    )
}

mod_file <- args[1]
module_prefix <- args[2]

# Validate input file exists
if (!file.exists(mod_file)) {
    stop("Input file does not exist: ", mod_file)
}

# Parse module file
tib_parsed <- parse_module_file(
    filename = mod_file, 
    split_cut_off_values = TRUE, 
    module_prefix = module_prefix
)

# Save parsed modules
saveRDS(tib_parsed, file = 'parsed_modules.rds', version = 2)

# Add debug information about the saved object
if (Sys.getenv("DEBUG") == "true") {
    write_delim(
        tib_parsed,
        file = "parsed_modules.tsv",
        delim = "\t"
    )
    # Print debug information about the saved object
    cat("Saved RDS file info:\n")
    cat("File size:", file.size("parsed_modules.rds"), "bytes\n")
    cat("Object structure:\n")
    str(tib_parsed)
}
