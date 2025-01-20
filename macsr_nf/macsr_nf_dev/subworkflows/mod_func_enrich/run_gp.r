#!/usr/bin/env Rscript
# run_gp.r
# Module functional enrichment and annotation using g:Profiler
# Performs enrichment analysis on pre-processed module data
# 14 01 2025
# refer: 
#	macsr_nf/macsr_nf_dev/subworkflows/mod_func_enrich/mfe_README.txt
# 	macsr_nf/macsr_nf_dev/subworkflows/mod_func_enrich/modferan.r

suppressPackageStartupMessages({
    library(tidyverse)
})

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    stop(
        "Six arguments must be supplied:\n",
        "1. pre_proc_mod_file: pre-processed module file (RDS output from mod_parse.r)\n", 
        "2. filter_cut_off: network cut-off value to filter on (string)\n",
        "3. sources: gProfiler databases to use for enrichment (comma-separated list)\n",
        "4. signif_level: gProfiler significance level (numeric)\n",
        "5. exclude_iea: exclude GO IEA terms (TRUE/FALSE)\n",
        "6. out_file: output filename\n",
        call. = FALSE
    )
}

# Assign arguments to named variables
pre_proc_mod_file <- args[1]
filter_cut_off <- args[2]
sources <- args[3]
signif_level <- as.numeric(args[4])
exclude_iea <- args[5]
out_file <- args[6]

# Validate inputs
if (!file.exists(pre_proc_mod_file)) {
    stop("Input file does not exist: ", pre_proc_mod_file)
}

if (is.na(signif_level) || signif_level <= 0 || signif_level > 1) {
    stop("Significance level must be between 0 and 1")
}

# Get script directory using Nextflow's environment variable
script_dir <- Sys.getenv("NXF_SCRIPT_DIR", unset = NA)
if (is.na(script_dir)) {
    stop("NXF_SCRIPT_DIR environment variable not set")
}

# Source functions file
source(file.path(script_dir, 'gprofiler_enrichment_functions.r'))

# Load and process module data
if (file.exists(pre_proc_mod_file)) {
    cat("Loading RDS file:", pre_proc_mod_file, "\n")
    cat("File size:", file.size(pre_proc_mod_file), "bytes\n")
    
    tryCatch({
        tib_parsed <- readRDS(pre_proc_mod_file)
        cat("Successfully loaded RDS file\n")
    }, error = function(e) {
        stop("Failed to read RDS file: ", e$message, "\n",
             "File path: ", pre_proc_mod_file, "\n",
             "File exists: ", file.exists(pre_proc_mod_file))
    })
} else {
    stop("Input RDS file does not exist: ", pre_proc_mod_file)
}

# Add debug output before gprofiler query
print("DEBUG: Converting module tibble to gprofiler query")
print(paste("DEBUG: Filter cutoff value:", filter_cut_off))
print(paste("DEBUG: Input tibble dimensions:", nrow(tib_parsed), "x", ncol(tib_parsed)))

gp_qry <- convert_module_tibble_to_gp_query(tib_parsed, filter_cut_off = filter_cut_off)

# Check if query is empty and exit gracefully
if (length(gp_qry) == 0) {
    cat("WARNING: No valid modules to process - skipping g:Profiler analysis\n")
    # Create empty output file to signal skip
    file.create(out_file)
    quit(status = 0)
}

# Add debug output after conversion
print("DEBUG: Conversion complete")
print(paste("DEBUG: Number of queries generated:", length(gp_qry)))

# Run g:Profiler enrichment analysis
exclude_iea_bool <- as.logical(tolower(exclude_iea))  # Convert string to logical more safely
gp_enrich <- run_gprofiler_enrichment(
    gp_qry,
    sig_level = signif_level,
    sources = strsplit(sources, ",")[[1]], 
    exclude_go_iea = exclude_iea_bool,
    multq = FALSE 
)

# Save results
saveRDS(gp_enrich, file = out_file)

# Debug output if enabled
if (Sys.getenv("DEBUG") == "true") {
    write_delim(
        gp_enrich$result,  # Access the result data frame
        file = paste0(out_file, ".tsv"),
        delim = "\t"
    )
}
