#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { CONVERT_IDS } from '../id_conversion/main.nf'

workflow {
    // Test parameters
    // params.test_network_file = "${params.nf_out_dir}/pre_processed_networks/TEST/HEAD3_network_modules_K1_humanbase.dat"
    params.test_network_file = "${params.nf_out_dir}/pre_processed_networks/network_modules_K1_humanbase.dat"
    params.test_network_db = 'humanbase'
    params.test_network_file_type = 'original'
    
    // Create input channels
    ch_network_file = Channel.fromPath(params.test_network_file)

    // Create mapping channel with correct mapping type
    ch_mapping = Channel.value([
        file("${params.map_hgnc}"),
        'hgnc',
        'ensembl_gene_id'
    ])

    // Run CONVERT_IDS workflow
    CONVERT_IDS(
        ch_network_file,
        params.test_network_db,
        params.test_network_file_type,
        ch_mapping
    )

    // // Collect and merge all converted files
    // CONVERT_IDS.out
    //     .collectFile(
    //         name: 'merged_converted_network.txt',
    //         storeDir: "${params.nf_out_dir}/converted_networks",
    //         skip: 1,  // Skip header lines
    //         keepHeader: true  // Keep the header from the first file
    //     )
    //     .view { merged_file ->
    //         """
    //         DEBUG: Merged converted network file created: 
    //         ${merged_file}
    //         """
    //     }
} 

workflow.onComplete {
    log.info """
    Workflow subworkflows/id_conversion/test_id_conversion.nf completed
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    """
}

workflow.onError {
    log.error """
    Test failed
    Error: ${workflow.errorMessage}
    """
} 