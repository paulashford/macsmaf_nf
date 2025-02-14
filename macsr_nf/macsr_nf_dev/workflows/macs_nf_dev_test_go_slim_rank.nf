#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GO_SLIM } from '../subworkflows/go_slim'
include { RANK_ANNOT } from '../subworkflows/rank_annot'

// Test workflow that runs GO_SLIM and RANK_ANNOT together
workflow {
    // Create a channel from the input enrichment results file (using RDS)
    ch_enrichment_results = Channel.fromPath(params.test_enrichment_file.replaceAll('\\.tsv$', '.rds'))
        .ifEmpty { error "No enrichment results RDS file found at: ${params.test_enrichment_file.replaceAll('\\.tsv$', '.rds')}" }
        // Parse the filename to extract method, db, and cutoff
        .map { file -> 
            def matcher = file.name =~ /([A-Z][0-9])_([a-z]+)_([0-9.]+)_enrichment_with_metrics\.rds$/
            if (matcher.matches()) {
                def (method, db, cutoff) = [matcher[0][1], matcher[0][2], matcher[0][3]]
                return tuple(method, db, cutoff, file)
            } else {
                error "File name doesn't match expected pattern: ${file.name}\nExpected pattern: {method}_{db}_{cutoff}_enrichment_with_metrics.rds\nExample: K1_cpdb_0.9_enrichment_with_metrics.rds"
            }
        }

    // Run GO_SLIM workflow with the full tuple
    GO_SLIM(
        ch_enrichment_results,  // Pass the full tuple (method, db, cutoff, file)
        params.go_slim_min_perc_rank ?: 0.25
    )

    // Run RANK_ANNOT workflow
    RANK_ANNOT(
        ch_enrichment_results,
        GO_SLIM.out.go_slim_results,
        params.max_term_size ?: 0.05
    )

    // Debug output
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
} 