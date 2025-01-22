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
include { RANK_ANNOT } from '../subworkflows/rank_annot'

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

    // Create channels from the parameters
    ch_network_databases = Channel.fromList(params.net_dbs ?: [])
                                .ifEmpty { error "No network types specified in params.net_dbs" }
    ch_analysis_methods = Channel.fromList(params.net_methods ?: [])
                                .ifEmpty { error "No network methods specified in params.net_methods" }
    ch_network_cutoffs = Channel.fromList(params.net_cutoffs ?: [])
                                .ifEmpty { error "No cutoff values specified in params.net_cutoffs" }

    // Create databases and methods pairs
    ch_net_dbs_methods = ch_network_databases
        .combine(ch_analysis_methods)
        .map { db, method -> [ method, db ] }

    // Validate gprofiler parameters and run functional enrichment
    validateGprofilerParams(
        params.gprofiler_sources,
        params.gprofiler_sig,
        params.gprofiler_exclude_iea
    )

    MOD_FUNC_ENRICH(
        ch_net_dbs_methods,
        ch_network_cutoffs,
        params.nf_network_modules_dir,
        params.net_file_prefix,
        params.gprofiler_sources,
        params.gprofiler_sig,
        params.gprofiler_exclude_iea
    )

    // Add GO_SLIM workflow using MOD_FUNC_ENRICH output
    GO_SLIM(
        MOD_FUNC_ENRICH.out.processed_enrichment,  // tuple(method, db, cutoff, file) with metrics
        params.go_slim_min_perc_rank ?: 0.25
    )

    // Add RANK_ANNOT workflow
    RANK_ANNOT(
        MOD_FUNC_ENRICH.out.processed_enrichment,  // tuple(method, db, cutoff, file) with metrics
        GO_SLIM.out.mapped_slim_results,           // tuple(method, db, cutoff, file) after mapping to slim
        MOD_FUNC_ENRICH.out.modules,               // tuple(method, db, cutoff, file) containing gene lists
        params.max_term_size ?: 0.05
    )

    // Debug output
    MOD_FUNC_ENRICH.out.processed_enrichment
        .view { method, db, cutoff, file -> 
            "DEBUG: Processed enrichment results:\n" +
            "  Method: ${method}\n" +
            "  Database: ${db}\n" +
            "  Cutoff: ${cutoff}\n" +
            "  File: ${file}"
        }

    GO_SLIM.out.go_slim_results
        .view { method, db, cutoff, file -> 
            "DEBUG: GO Slim results:\n" +
            "  Method: ${method}\n" +
            "  Database: ${db}\n" +
            "  Cutoff: ${cutoff}\n" +
            "  File: ${file}"
        }
    
    RANK_ANNOT.out.ranked_results
        .view { method, db, cutoff, file -> 
            "DEBUG: Ranked and annotated results:\n" +
            "  Method: ${method}\n" +
            "  Database: ${db}\n" +
            "  Cutoff: ${cutoff}\n" +
            "  File: ${file}"
        }

    // Output for checking and validation
    MOD_FUNC_ENRICH.out.modules
        .view { file -> 
            "DEBUG: Parsed modules output:\n" +
            "  File: ${file}"
        }

    MOD_FUNC_ENRICH.out.processed_enrichment
        .view { file -> 
            "DEBUG: Processed enrichment results:\n" +
            "  File: ${file}"
        }
}

