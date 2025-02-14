#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GO_SLIM } from '../subworkflows/go_slim'

// Test workflow that only runs GO_SLIM
workflow {
    // Create a channel from the input enrichment results file
    ch_enrichment_results = Channel.fromPath(params.test_enrichment_file)
        .ifEmpty { error "No enrichment results file found at: ${params.test_enrichment_file}" }

    // Run GO_SLIM workflow
    GO_SLIM(
        ch_enrichment_results,
        params.go_slim_min_perc_rank ?: 0.25
    )

    // Debug output
    GO_SLIM.out.go_slim_results
        .view { file -> 
            "DEBUG: GO Slim results:\n" +
            "  File: ${file}"
        }
    
    GO_SLIM.out.mapped_slim_results
        .view { file -> 
            "DEBUG: Mapped GO Slim results:\n" +
            "  File: ${file}"
        }
}