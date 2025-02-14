// id_conversion.nf
// Subworkflow for converting gene IDs eg used in different net databases

process convert_gene_ids {
    //debug { params.debug }
    tag "${network_db}:${network_file_type}"
    publishDir "${params.id_conversion_out_dir}", mode: 'copy'

    // Activate virtual environment 
    beforeScript "source ${params.root_proj_dir}/venv/bin/activate"

    input:
    tuple val(module), val(mapping)
    val network_db
    val network_file_type
    
    output:
    path 'converted_ids', emit: converted

    script:
    def (module_id, ids_to_convert) = module
    def (mapping_file, map_file_type, mod_id_type_out, id_type) = mapping
    """
    export DEBUG=${params.debug}
    macsr_nf_dev convert-ids-with-mapfile \\
        --map_file $mapping_file \\
        --map_file_type $map_file_type \\
        --ids_to_convert '$ids_to_convert' \\
        --id_type '$id_type' \\
        --approved_only \\
        --col_filter \\
        --outfile 'converted_id_map_table.tsv' \\
        --id_type_out $mod_id_type_out > converted_ids1
    echo "$module_id""\t"\$( <converted_ids1 ) | tr ' ' '\\t' | sed 's/nan//g' | tr -s '[:blank:]' > converted_ids
    """
}

process preproc_module_file {
    // This just makes the 1st tab ~ so that splitCsv will give 2 columns: (module_id, tab-sep-gene-list)
    // network-modules-K1-cpdb-0.0.dat:4    ENSG00000169962 ENSG00000173662 ENSG00000179002 --->
    // network-modules-K1-cpdb-0.0.dat:4~ENSG00000169962    ENSG00000173662 ENSG00000179002
    input:
    path module_file
    val sed_pp_modules_lines_script

    output:
    path 'module_file_preproc.dat'

    script:
    """
    # Split on first tab, preserve rest of line
    awk -F'\\t' '{
        module_id=\$1
        \$1=""  # Remove first field
        sub(/^[ \\t]+/, "", \$0)  # Remove leading whitespace
        print module_id"~"\$0
    }' $module_file > module_file_preproc.dat
    """
}

process determine_id_type {
    // Activate virtual environment 
    beforeScript "source ${params.root_proj_dir}/venv/bin/activate"

    input:
    val network_db
    val network_file_type

    output:
    stdout

    script:
    def dict_type = network_file_type == 'original' ? 'net_id_type_orig' : 'net_id_type_pre_proc'
    """
    macsr_nf_dev dict-lookup --dict_type $dict_type --dict_key $network_db
    """
}

workflow PREPROC_MODULE_FILE {
    take:
    ch_modfiles
    sed_pp_modules_lines_script
    
    main:
    ch_ppmod    =   preproc_module_file(ch_modfiles, sed_pp_modules_lines_script )

    //ch_labelled     = label_network_module_files( ch_nets ) 
    
    emit:
    ch_ppmod
        
}

workflow CONVERT_IDS {
    take:
        out_dir         // Channel
        id_map_to       // Channel
        input_network_method  // Channel
        network_file    // Channel
        network_db      // Channel
        network_file_type  // Channel
        ch_mapping     // Channel

    main:
        if (params.debug) {
            log.info "DEBUG: Starting ID conversion workflow"
            log.info "DEBUG: Network DB: $network_db"
            log.info "DEBUG: Network file type: $network_file_type"
        }

        // Get ID type from Python constants
        nf_mod_id_type = determine_id_type(network_db, network_file_type)
            .map { it.trim() }

        // Preprocess module files
        ch_pp = PREPROC_MODULE_FILE(
            network_file,
            params.sed_pp_modules_lines_script
        )

        // Process module lines and create a channel for each module
        ch_modules = ch_pp
            .splitCsv(header: false, sep: '~')
            .map { row -> 
                def (module_id, gene_list) = row
                if (params.debug) {
                    log.info "DEBUG: Processing module: ${module_id}"
                }
                [module_id, gene_list.trim()]
            }

        // Prepare mapping parameters
        ch_mapping_with_id = ch_mapping
            .combine(nf_mod_id_type)
            .map { map_file, map_type, target_type, id_type -> 
                [map_file, map_type, target_type, id_type]
            }

        // Convert IDs for each module by combining with mapping info
        ch_converted = convert_gene_ids(
            ch_modules.combine(ch_mapping_with_id)
                .map { module_id, genes, map_file, map_type, target_type, id_type ->
                    [ [module_id, genes], [map_file, map_type, target_type, id_type] ]
                },
            network_db,
            network_file_type
        )

        ch_converted
            .collectFile(   
                name: "network_modules_" + input_network_method + "_" + network_db + "_converted_" + network_file_type + "_to_" + id_map_to + ".dat",
                storeDir: file(out_dir)
            )
            // .tap { 
            //     if (params.debug) {
            //         log.info "DEBUG: collectFile output: ${it}"
            //     }
            // }
        
    emit:
        ch_converted
}
