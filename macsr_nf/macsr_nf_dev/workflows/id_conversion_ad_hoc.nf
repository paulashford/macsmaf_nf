#!/usr/bin/env nextflow
// Simple workflow to process 1 network at a time to convert gene IDs
// Not part of the main workflow - can be run ad-hoc; examples with comments included below.
// Feb 2025 
// p.ashford@ucl.ac.uk
// export proj_root=~/git/macsmaf/macsr_nf/macsr_nf_dev
// export NF_CONFIG="${proj_root}/conf/base.config"
// source "${proj_root}venv/bin/activate"
// cd ${proj_root}/workflows
// nextflow run id_conversion_ad_hoc.nf -c "${NF_CONFIG}"

nextflow.enable.dsl=2
include { CONVERT_IDS } from '../subworkflows/id_conversion/main.nf'

// ad-hoc parameters
	// what is geneID type in network to convert? (i.e. is it the original net (i.e. source IDs defined by dict), or an already processed net, etc)
    // params.input_network_file_type = 'original' // (alternatively: 'net_id_type_orig')
	
	// params.input_network_db = 'humanbase'
	// params.input_network_method = 'K1'
	// params.input_network_file = "${params.nf_out_dir}/pre_processed_networks/network_modules_K1_humanbase.dat"
    // what map file to use for translating IDs
	// params.map_file = "${params.map_hgnc}"
	// params.map_file_type = 'hgnc' // alternatively 'biomart' (also 'uniprotkb' if need uniprot acc) - defined in dict MAPPING_FILE_TYPES
	// what is the type of geneID to convert to?
	// params.id_map_to = 'ensembl_gene_id'
	// params.id_output_dir = "${params.id_conversion_out_dir}/${params.map_file_type}"

	// params.input_network_db = 'string'
	// params.input_network_method = 'K1'
	// params.input_network_file = "${params.nf_out_dir}/pre_processed_networks/network_modules_K1_string.dat"
	// // what map file to use for translating IDs
	// params.map_file = "${params.map_biomart}"
	// params.map_file_type = 'biomart' 
	// params.id_map_to = 'Gene_stable_ID' // ensembl_gene_id
	// params.id_output_dir = "${params.id_conversion_out_dir}/${params.map_file_type}"

	// params.input_network_db = 'cpdb' (not required for orignal geneID to ENSG as CPDB uses ENSG IDs in its network)
	
	// Run for a channel of files in input dir
	params.input_network_dir = "${params.nf_out_dir}/compiled_module_datasets/v01/pre_processed_networks"
	// params.input_network_dir = "${params.nf_out_dir}/pre_processed_networks"
	params.input_network_file_type = 'net_id_type_pre_proc'
	params.map_file = "${params.map_hgnc}"
	params.map_file_type = 'hgnc'
	params.id_map_to = 'symbol'
	params.id_output_dir = "${params.id_conversion_out_dir}/${params.map_file_type}"

workflow {
    // Add debug parameter explicitly
    params.debug = params.debug ?: false

    // Ensure output directory exists
    file(params.id_output_dir).mkdir()

    // Create input channel for all .dat files in pre_processed_networks directory
    ch_network_files = Channel
        .fromPath("${params.input_network_dir}/*.dat")
        .map { file -> 
            // Extract method and db from filename (assuming format: network_modules_K1_humanbase.dat)
            def parts = file.name.toString().split('_')
            def method = parts[2]  // K1
            def db = parts[3].replace('.dat', '')  // humanbase
            return tuple(file, method, db)
        }

    // Create mapping channel with correct mapping type
    ch_mapping = Channel.value(
        [
            file(params.map_file),
            params.map_file_type,
            params.id_map_to
        ] 
    )

    // Run CONVERT_IDS workflow for each file
    CONVERT_IDS(
        Channel.value(params.id_output_dir),
        Channel.value(params.id_map_to),
        ch_network_files.map{ it[1] },  // method
        ch_network_files.map{ it[0] },  // file
        ch_network_files.map{ it[2] },  // db
        Channel.value(params.input_network_file_type),
        ch_mapping
    )

	// onError {
	// 	log.error """
	// 		ID conversion failed for ${params.input_network_file} using ${params.map_file_type} map file to convert to ${params.id_map_to}.
	// 		Error: ${workflow.errorMessage}
	// 	"""
	// } 
	// onComplete {
	// 	log.info """
	// 		Workflow subworkflows/id_conversion/input_id_conversion.nf completed using ${params.map_file_type} map file to convert to ${params.id_map_to}.
	// 		Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}	
	// 		Duration: ${workflow.duration}
	// 	"""
	// }

	
}
// ----------------------------------------------------------------------------------------------------
// PARAMETER VALUES - EXAMPLE COPIES BELOW FOR REFERENCE ONLY - USE the .config file or constants.py
// ----------------------------------------------------------------------------------------------------
// Mapping files should be in .config file  - copies below for reference:
//	params {
//		// ID mapping biomart for EnsemblIDs and gene/entrez IDs etc
//		map_biomart					=	"${params.nf_datasets_dir}/biomart/biomart01/mart_export.txt"
//		// ID map using HGNC table (overlaps biomart for some IDs eg entrez- use biomart primarily unless HGNC specifc)
//		map_hgnc					= 	"${params.nf_datasets_dir}/hgnc/hgnc_complete_set.txt"
//		// ID mapping tables for genes IDs and UniProt
//		// map_uniprot 				= 	"${params.nf_datasets_dir}/uniprot/HUMAN_9606_idmapping_cut_cols_with_header.tab"
//	}

// Valid input parameters defined here: macsr_nf_dev/macsr_nf_dev/constants.py
// copies below for reference:
// MAPPING_FILE_TYPES = ['hgnc', 'biomart', 'uniprotkb'] (choose 1 that contains both the geneID to convert from and geneID to convert to)
// HGNC_VALID_TYPES = ['entrez_id', 'ensembl_gene_id', 'hgnc_id', 'symbol', 'refseq_accession' ]
// BIOMART_VALID_TYPES = ['Gene_stable_ID', 'Gene_stable_ID_version', 'Transcript_stable_ID', 'Transcript_stable_ID_version', 'Protein_stable_ID', 'Protein_stable_ID_version', 'NCBI_gene_ID', 'Gene_description', 'HGNC_ID', 'UniProtKB_ID']
// UNIPROTKB_VALID_TYPES = ['UniProtKB-AC', 'UniProtKB-ID', 'EntrezGeneID', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional_PubMed']
// dict for lookup of type of gene IDs for each network db
// NET_DICT_TYPES = {
//     'net_id_type_orig': NET_DB_ORIGINAL_GENE_ID_TYPE_DICT,
//     'net_id_type_pre_proc': NET_DB_PRE_PROCESSED_GENE_ID_TYPE_DICT,
//     'id_map_type': MAPPING_TYPE_ESSENTIAL_COLS_DICT
// }

