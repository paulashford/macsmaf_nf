#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/workflows
// export NF_CONFIG=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/conf/base.config
// nextflow run macs_nf_dev.nf -c "${NF_CONFIG}"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include {PASCAL_GWAS} from '../macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/main.nf'
include { PASCAL_GWAS } from '../subworkflows/pascal_gwas'
include { CONVERT_IDS } from '../subworkflows/id_conversion'
include { GO_SLIM } from '../subworkflows/go_slim'
include { MOD_FUNC_ENRICH } from '../subworkflows/mod_func_enrich'

// Add validation functions
def validateModFuncEnrichInputs(net_methods, net_dbs, net_cutoffs, modules_dir, file_prefix) {
    // Validate methods
    if (!net_methods) {
        error "Network methods cannot be null"
    }
    // Validate dbs
    if (!net_dbs) {
        error "Network dbs cannot be null"
    }
    // Validate cutoffs
    if (!net_cutoffs) {
        error "Network cutoffs cannot be null"
    }
    // Verify cutoffs strings can be parsed as valid numerics
    if (!net_cutoffs.every { it.toString().isNumber() }) {
        error "Network cutoffs must be valid numeric values"
    }

    // Validate directory paths
    if (!modules_dir) {
        error "Network modules directory path cannot be null"
    }
    if (!file(modules_dir).exists()) {
        error "Network modules directory does not exist: ${modules_dir}"
    }

    // Validate file prefix
    if (!file_prefix) {
        error "Network file prefix cannot be null"
    }
}

def validateGprofilerParams(sources, sig_threshold, exclude_iea) {
    // Validate sources
    if (!sources) {
        error "gprofiler_sources parameter is required"
    }

    // Validate significance threshold
    if (sig_threshold == null) {
        error "gprofiler_sig parameter is required"
    }
    if (!(sig_threshold instanceof Number)) {
        error "gprofiler_sig must be a number"
    }
    if (sig_threshold < 0 || sig_threshold > 1) {
        error "gprofiler_sig must be between 0 and 1"
    }

    // Validate IEA exclusion parameter
    if (exclude_iea == null) {
        error "gprofiler_exclude_iea parameter is required"
    }
    if (!(exclude_iea instanceof Boolean)) {
        error "gprofiler_exclude_iea must be a boolean"
    }
}

def validateOutputPaths(output_dir) {
    // Validate output directory
    if (!output_dir) {
        error "Output directory path cannot be null"
    }
    
    def output_path = file(output_dir)
    if (!output_path.exists()) {
        log.info "Creating output directory: ${output_dir}"
        if (!output_path.mkdirs()) {
            error "Failed to create output directory: ${output_dir}"
        }
    }
    if (!output_path.isDirectory()) {
        error "Output path exists but is not a directory: ${output_dir}"
    }
    if (!output_path.canWrite()) {
        error "Output directory is not writable: ${output_dir}"
    }
}

workflow {
    log.info """
    Parameters loaded:
    net_dbs: ${params.net_dbs}
    net_methods: ${params.net_methods}
    net_cutoffs: ${params.net_cutoffs}
    """

    // Validate inputs
    validateModFuncEnrichInputs(
        params.net_methods,
        params.net_dbs,
        params.net_cutoffs,
        params.nf_network_modules_dir,
        params.net_file_prefix
    )

    // Validate gprofiler parameters
    validateGprofilerParams(
        params.gprofiler_sources,
        params.gprofiler_sig,
        params.gprofiler_exclude_iea
    )

    // Validate output paths
    validateOutputPaths(params.nf_enrichment_dir)   

    // Create channels from the parameters
    ch_network_databases = Channel.fromList(params.net_dbs ?: [])
    ch_analysis_methods = Channel.fromList(params.net_methods ?: [])
    ch_network_cutoffs = Channel.fromList(params.net_cutoffs ?: [])

    // Check if channels are empty
    ch_network_databases.ifEmpty { error "No network types specified in params.net_dbs" }
    ch_analysis_methods.ifEmpty { error "No network methods specified in params.net_methods" }
    ch_network_cutoffs.ifEmpty { error "No cutoff values specified in params.net_cutoffs" }

    // Create databases and methods pairs
    ch_net_dbs_methods = ch_network_databases
        .combine(ch_analysis_methods)
        .map { db, method -> [ method, db ] }

    // Execute process with validated inputs
    MOD_FUNC_ENRICH(
        ch_net_dbs_methods,
        ch_network_cutoffs,
        params.nf_network_modules_dir,
        params.net_file_prefix,
        params.gprofiler_sources,
        params.gprofiler_sig,
        params.gprofiler_exclude_iea
    )

    // Add GO_SLIM workflow call using MOD_FUNC_ENRICH output
    GO_SLIM(
        MOD_FUNC_ENRICH.out.enrichment_results,
        params.go_slim_min_perc_rank ?: 0.25  // Default to 0.25 if not specified
    )

    // Output for checking and validation
    MOD_FUNC_ENRICH.out.modules
        .view { file -> 
            "DEBUG: Parsed modules output:\n" +
            "  File: ${file}"
        }

    MOD_FUNC_ENRICH.out.enrichment_results
        .view { file -> 
            "DEBUG: Enrichment results:\n" +
            "  File: ${file}"
        }
}

