// subworkflow for g:Profiler functional enrichment
nextflow.enable.dsl=2
workflow MOD_FUNC_ENRICH {
	take:
		ch_net_dbs_methods
		ch_net_cutoffs
		net_modules_dir
		net_file_prefix
		gprofiler_sources
		gprofiler_sig
		gprofiler_exclude_iea
	
	main:
		// Get network files and validate they exist
		network_files = get_network_filenames(
			ch_net_dbs_methods,
			net_file_prefix,
			net_modules_dir
		)
		
		// Convert output to a channel of files
		network_file = network_files
			.splitText()   
			.map { it.trim() } 
			.filter { it.length() > 0 } 
			.map { file(it) }  
			.ifEmpty { error "No module network files found matching pattern" }
		
		// Parse modules and combine with method/db info
		parsed_modules = parse_modules(network_file, net_file_prefix)
			.combine(ch_net_dbs_methods)
			.combine(ch_net_cutoffs)
			.map { modules, method, db, cutoff -> 
				tuple(method, db, cutoff, modules)
			}
		
		// Run gprofiler with combined inputs
		gprofiler_results = run_gprofiler(
			parsed_modules,          // tuple(method, db, cutoff, modules)
			gprofiler_sources,       // val(sources)
			gprofiler_sig,          // val(signif_level)
			gprofiler_exclude_iea    // val(exclude_iea)
		)
	
		// Post process enrichment results
		processed_enrichment = post_process_gprofiler_enrichment(
			gprofiler_results
		)

		// Add classification metrics and MCC
		processed_enrichment = add_tp_tn_fp_fn_mcc(
			processed_enrichment
		)


	emit:
		modules = parsed_modules
		enrichment_results = gprofiler_results
		processed_enrichment = processed_enrichment.processed_results
}

// get_network_filenames: for given netdb type (string, cpdb, ... ) and method (K1, ...) get matching networks 
// available in net_modules_dir. Use of glob and files() returns a list, so if multiple types of network for 
// a given netdb / method (e.g. each using different gene identifiers), then each will be processed through channel
process get_network_filenames {
	debug true
	input:
		val(netdb_method)
		val(net_file_prefix)
		val(net_modules_dir)

	output:
		path 'net_file_names'

	script:
		"""
		for file in ${net_modules_dir}/${net_file_prefix}${netdb_method[0]}_${netdb_method[1]}*.dat; do
			if [ -f "\$file" ]; then
				echo "\$file" >> net_file_names
			fi
		done
		if [ ! -s net_file_names ]; then
			echo "No matching files found for pattern: ${net_modules_dir}/${net_file_prefix}${netdb_method[0]}_${netdb_method[1]}*.dat" >&2
			exit 1
		fi
		"""
}

process parse_modules {
	debug true    // this will print the stdout from the script section on Terminal
	publishDir "${params.nf_out_dir}", mode: 'copy'

	input:
		path network_file
		val net_file_prefix

	output:
		path "parsed_modules.rds"

	script:
	"""
	export NXF_SCRIPT_DIR="${params.root_proj_dir}/subworkflows/mod_func_enrich"
	Rscript "${params.root_proj_dir}/subworkflows/mod_func_enrich/mod_parse.r" \\
		"${network_file}" \\
		"${net_file_prefix}"
	"""
}

process run_gprofiler {
	debug true
	publishDir "${params.nf_out_dir}/enrichment", mode: 'copy'

	input:
		tuple val(method), val(db), val(cutoff), path(parsed_modules)
		val(sources)
		val(signif_level)
		val(exclude_iea)

	output:
		tuple val(method), val(db), val(cutoff), path("${out_prefix}_gp_enrich.rds"), emit: enrichment_results

	script:
	out_prefix = "${method}_${db}_${cutoff}"
	"""
	# First check if input file is empty
	if [ ! -s ${parsed_modules} ]; then
		echo "WARNING: Empty parsed modules file - skipping g:Profiler analysis"
		touch ${out_prefix}_gp_enrich.rds
		exit 0
	fi

	Rscript ${params.root_proj_dir}/subworkflows/mod_func_enrich/run_gp.r \\
		${parsed_modules} \\
		${cutoff} \\
		${sources} \\
		${signif_level} \\
		${exclude_iea.toString().toLowerCase()} \\
		${out_prefix}_gp_enrich.rds
	"""
}

process post_process_gprofiler_enrichment {
	debug true
	publishDir "${params.nf_out_dir}/enrichment", mode: 'copy'

	input:
		tuple val(method), val(db), val(cutoff), path(gprofiler_results)

	output:
		tuple val(method), val(db), val(cutoff), path("${out_prefix}_processed_enrichment.rds"), emit: processed_results

	script:
	out_prefix = "${method}_${db}_${cutoff}"
	"""
	# First check if input file is empty
	if [ ! -s ${gprofiler_results} ]; then
		echo "WARNING: Empty g:Profiler results file - skipping post-processing"
		touch ${out_prefix}_processed_enrichment.rds
		exit 0
	fi

	Rscript ${params.root_proj_dir}/subworkflows/mod_func_enrich/post_process_gp.r \\
		${gprofiler_results} \\
		${method} \\
		${db} \\
		${cutoff} \\
		${out_prefix}_processed_enrichment.rds
	"""
}

process add_tp_tn_fp_fn_mcc {
	debug true
	publishDir "${params.nf_enrichment_dir}", mode: 'copy'

	input:
		tuple val(method), val(db), val(cutoff), path(processed_enrichment)

	output:
		tuple val(method), val(db), val(cutoff), path("${out_prefix}_enrichment_with_metrics.rds"), emit: processed_results
		tuple val(method), val(db), val(cutoff), path("${out_prefix}_enrichment_with_metrics.tsv"), emit: enrichment_results

	script:
	out_prefix = "${method}_${db}_${cutoff}"
	"""
	# First check if input file is empty
	if [ ! -s ${processed_enrichment} ]; then
		echo "WARNING: Empty processed enrichment file - skipping metrics calculation"
		touch ${out_prefix}_enrichment_with_metrics.rds
		touch ${out_prefix}_enrichment_with_metrics.tsv
		exit 0
	fi
	
	Rscript ${params.root_proj_dir}/subworkflows/mod_func_enrich/add_metrics.r \
		${processed_enrichment} \
		"${out_prefix}_enrichment_with_metrics.rds" \
		"${out_prefix}_enrichment_with_metrics.tsv"
	"""
}