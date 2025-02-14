17 01 2025 g:Profiler enrichment for K1 modules
Successful run following debugging macsr_nf_dev
NOTE: 
	api_base_url='https://biit.cs.ut.ee/gprofiler_beta'
	in: subworkflows/mod_func_enrich/gprofiler_enrichment_functions.r
	
	Now amended so can pass as NXF sys variable (or in run_gprofiler_enrichment function call)
	# g:Profiler API base URL
	nxf_api_base_url <- Sys.getenv("NXF_GPROFILER_API_URL", unset = NA)
	if (!is.na(nxf_api_base_url)) {
		api_base_url <- nxf_api_base_url
	}

cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/workflows
export NF_CONFIG=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/conf/base.config
nextflow run macs_nf_dev.nf -c "${NF_CONFIG}"

Completed at: 17-Jan-2025 16:44:21
Duration    : 24m 1s
CPU hours   : 6.0
Succeeded   : 115


--------------
params {
	config_profile_name        = 	'MACSMAF macsr_nf_dev'
    config_profile_description = 	'Parameters for network module functional enrichment, ranking and annotation'

	// user dir
	user_root_dir				= 	'/Users/ash'
	
	// Path to general datasets (eg for ID tables)
	nf_datasets_dir				= 	"${params.user_root_dir}/data/funvar_pipeline/datasets"

	// git proj dir
    root_proj_dir				= 	"${params.user_root_dir}/git/macsmaf/macsr_nf/macsr_nf_dev"
	script_dir					= 	"${params.root_proj_dir}/script"
	nf_out_dir					=	"${params.root_proj_dir}/output"
	// outdir  					= "${projectDir}/output"

	// script refs
	sed_simplify_script			= 	"${params.script_dir}/simplify_label.sed"
	sed_split_label_script		= 	"${params.script_dir}/split_label_cutoff.sed"
	sed_pp_modules_lines_script	=	"${params.script_dir}/preproc_module_lines.sed"
	
	// Networks
	net_dbs						= 	[ 'cpdb', 'humanbase', 'string' ]
	net_methods					=	[ 'K1' ]
	// net_methods					=	[ 'M1', 'R1', 'K1' ]
	net_cutoffs					=	[ '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9' ]
	// Path to the DREAM modules, original IDs (nextflow datasets directory, which may link to others)
	nf_network_modules_dir		= 	"${params.user_root_dir}/Dropbox/bioinf/MACSMAF/datasets/d024/id_conversion"
	net_file_prefix 			= 	'network_modules_'


	//g:Profiler params
	gprofiler_sources			=	"GO:BP,KEGG,REAC"
	gprofiler_sig				=	0.01
	gprofiler_exclude_iea		=	true


	// ID mapping biomart for EnsemblIDs and gene/entrez IDs etc
	map_biomart					=	"${params.nf_datasets_dir}/biomart/biomart01/mart_export.txt"
	// ID map using HGNC table (overlaps biomart for some IDs eg entrez- use biomart primarily unless HGNC specifc)
	map_hgnc					= 	"${params.nf_datasets_dir}/hgnc/hgnc_complete_set.txt"
	// ID mapping tables for genes IDs and UniProt
	map_uniprot 				= 	"${params.nf_datasets_dir}/uniprot/HUMAN_9606_idmapping_cut_cols_with_header.tab"

	
	// exper_id
	// nf_exper_id					= 'e023'
	// string: the net type {string, cpdb or humanbase}
    // nf_network_type 			= 'string'
    // string: the DREAM net module detect algorithm {M1, R1, K1}
	// nf_module_detect 			= 'K1'
    //  string: The output directory where the results will be saved
	
	// Gene Ontology (GO) and GOslims
	go_dir						= 	"${params.nf_datasets_dir}/go"
	go_obo						=	"${params.go_dir}/go-basic.obo"
	go_slim						=	"${params.go_dir}/goslim_pir.obo"
	gaf_human					=	"${params.go_dir}/gaf_human.gaf"
	// owltools exec
	// owl_tools					= 	"/Users/ash/git/bioinf/owltools/owltools"
	
	// ID MAP
	// nf_id_mapping_file			= 

	// nf_module_prefix			= "${params.nf_exper_id}_${params.nf_module_detect}_${nf_network_type}_mod_"
	

	
    // input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    
}

	