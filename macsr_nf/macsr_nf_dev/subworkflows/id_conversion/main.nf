//
// include { PREPROC_NET_MODULES } from '../pascal_gwas/main.nf'
// Subworkflow for converting gene IDs eg used in different net databases

// NextFlow submodule testing:
// cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/id_conversion
// export NF_CONFIG=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/conf/base.config
// nextflow run main.nf -c "${NF_CONFIG}"

// subworkflow test
// cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/id_convert
// export NF_CONFIG=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/conf/base.config
// nextflow run main.nf -c "${NF_CONFIG}"

// Hi Ash, I looked at the file I got from you before and it was the same format for the GWAS data. 
// The only thing was that it contained all the modules from all methods. I'm attaching it below for your reference. Could you please generate the same file but with all three methods and networks (and the same naming system)?

// Name    chi2Pvalue  empPvalue
// ENSG-R1-humanbase-coval-0.3.dat_158 3.34001336E-4   3.42E-4
// ENSG-R1-humanbase-coval-0.5.dat_122 3.34001336E-4   3.36E-4
// ENSG-M1-cpdb-coval-0.8.dat_202  5.1915018E-4    1.88E-3
// ENSG-R1-humanbase-coval-0.4.dat_239 5.73738664E-4   1.72E-4
// ENSG-R1-cpdb-coval-0.5.dat_49   8.42244723E-4   2.63E-3
// ENSG-M1-cpdb-coval-0.8.dat_89   8.78298669E-4   1.158E-2
// ENSG-M1-cpdb-coval-0.1.dat_98   1.20799046E-3   8.31E-3
// ENSG-R1-humanbase-coval-0.3.dat_67  1.56715245E-3   1.12E-3
// ENSG-M1-cpdb-coval-0.7.dat_128  1.81988014E-3   4.72E-3

// combine netwokr module files in dir (eg for individual module files with different edge threshold cut-offs)
// create a 'labelled' net module file with filename in 1st col, so that file now contains metadata info from path name
// eg: grep -R "${dataset_dir}" --include="*.dat" -e '' > "${exper_dir}/k1_modules_grep.txt" 2>/dev/null

process convert_ids {
    input:
    // module lines
    tuple val(module_id), val(ids_to_convert)
    // mapping params
    // Channel.value( [ file(nf_map_file), nf_map_type, nf_mod_id_type, nf_mod_id_out_type ] )
    tuple path(mapping_file), val(map_file_type), val(mod_id_type), val(mod_id_type_out)
    
    output:
    // path 'converted_id_map_table.tsv'
    path 'converted_ids'

    script:
    """
    macsr_nf_dev convert-ids-with-mapfile \
        --map_file $mapping_file \
        --map_file_type $map_file_type \
        --ids_to_convert $ids_to_convert \
        --id_type $mod_id_type \
        --approved_only \
        --col_filter \
        --outfile 'converted_id_map_table.tsv' \
        --id_type_out $mod_id_type_out > converted_ids1
    echo $module_id"\t"\$( <converted_ids1 ) | tr ' ' '\t' | sed 's/nan//g' | tr -s '[:blank:]' > converted_ids
    """
    // printf "%s %s" $module_id \$tmp > converted_ids
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
    //sed -E 's/(.*\\_[0-9]*)(\\t)/\\1~/' $module_file > module_file_preproc.dat
    """
    $sed_pp_modules_lines_script $module_file > module_file_preproc1.dat
    tr '\t' ',' < module_file_preproc1.dat > module_file_preproc.dat
    """
}

workflow CONVERT_IDS {
    
    main:
    // mapping file/type 
    // nf_map_file = "${params.map_biomart}"
    // nf_map_type = 'biomart'
    nf_map_file = "${params.map_hgnc}"
    nf_map_type = 'hgnc'

    // The module file(s) to process and the ID type used for module genes
    sed_pp_modules_lines_script = "${params.sed_pp_modules_lines_script}"
    // nf_mod_files            =   "/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/output/d024_network_modules_K1_cpdb.dat"
    // nf_mod_files            =   "/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/output/d024_network_modules_K1_humanbase.dat"
    // nf_mod_files            =   "/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/output/d024_network_modules_K1_string.dat"
    // For HUGO conv...
    // nf_mod_files    =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/cpdb/network_modules_K1_cpdb_ENSG.dat"
    // nf_mod_files    =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/humanbase/network_modules_K1_humanbase_ENSG.dat"
    // nf_mod_files        =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/humanbase/network_modules_K1_humanbase_ENSG_del_empty_modules.dat"
    nf_mod_files        =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/string/network_modules_K1_string_ENSG_del_empty_modules.dat"

    mod_file                =   file(nf_mod_files)
    // what are gene IDs in module file?
    // nf_mod_id_type          =   'entrez_id'
    nf_mod_id_type          =   'ensembl_gene_id'
    // nf_mod_id_type          =   'Protein_stable_ID'
    // what gene IDs do we convert to?
    // nf_mod_id_out_type      =   'ensembl_gene_id'
    nf_mod_id_out_type      =   'symbol'
    // nf_mod_id_out_type      =   'Gene_stable_ID'
    // composite_mod_filename  =   mod_file.baseName + "_pp" + ".dat"
    composite_mod_filename  =   mod_file.baseName + "_HUGO" + ".dat"
    // The module file(s)' Channel
    ch_modfiles     = Channel.fromPath( nf_mod_files )    
    ch_pp = PREPROC_MODULE_FILE( ch_modfiles, sed_pp_modules_lines_script )    
    
    // mapping Channel
    ch_mapping  = Channel.value( [ file(nf_map_file), nf_map_type, nf_mod_id_type, nf_mod_id_out_type ] )
    // The module lines in a module file
    ch_modules  =  ch_pp
                        .splitCsv( header:false, sep:'~' )
                        .map { row -> [ row[0], row[1] ] }

    ch_conv = CONVERT_GENE_IDS( ch_modules, ch_mapping )

    ch_conv
        .collectFile(   name: composite_mod_filename, 
                        storeDir: "${workflow.launchDir}/output"
        )
        .subscribe  { str -> println " PROCESS_MODULES: ${str}" }
    
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

workflow CONVERT_GENE_IDS {
    take:
    ch_modules
    ch_mapping

    main:
    ch_conv = convert_ids( ch_modules, ch_mapping )
                // .combine()

    emit:
    ch_conv
}
