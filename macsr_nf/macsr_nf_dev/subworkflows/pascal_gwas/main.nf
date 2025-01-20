// Subworkflow for running the Pascal / GWAS tests for sets of modules
nextflow.enable.dsl=2
// TODO needs the PASCAL runner process - this is not implemented here
// cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas
// export NF_CONFIG=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/conf/base.config
// nextflow run main.nf -c "${NF_CONFIG}"

// Hi Ash, I looked at the file I got from you before and it was the same format for the GWAS data. 
// The only thing was that it contained all the modules from all methods. 
// I'm attaching it below for your reference. 
// Could you please generate the same file but with all three methods and networks (and the same naming system)?

// Name    chi2Pvalue  empPvalue
// ENSG-R1-humanbase-coval-0.3.dat_158 3.34001336E-4   3.42E-4
// ENSG-R1-humanbase-coval-0.5.dat_122 3.34001336E-4   3.36E-4

// combine netwokr module files in dir (eg for individual module files with different edge threshold cut-offs)
// create a 'labelled' net module file with filename in 1st col, so that file now contains metadata info from path name
// eg: grep -R "${dataset_dir}" --include="*.dat" -e '' > "${exper_dir}/k1_modules_grep.txt" 2>/dev/null

process label_network_module_files {
    input:
    path(network_module_file)

    output:
    path 'labelled_network_modules_file.dat'

    script:
    """
    grep -H ${network_module_file} -e '' > 'labelled_network_modules_file.dat'
    """
}

// DREAM/Monet outputs format -> simpler format in the combined network module file (ie the 1st column) eg
// /path/to/datasets/humanbase/2023-04-20-095041__K1__result-modules__network-data-humanbase-coval-0.6.dat:1	1.0	4071	857	11030	1277	7168	8572	8614
// --->
// network-modules-K1-humanbase-0.6.dat:1	1.0	4071	857	11030	1277	7168	8572	8614
// sed -E 's/^.*\_\_([M|R|K][1])\_\_.*(cpdb|humanbase|string).*([0-9].[0-9]).dat/network-modules-\1-\2-\3.dat/g' labelled_network_modules_file > labelled_network_module_file_simplified.dat 
// sed -E 's/^.*\_\_([M|R|K][1])\_\_.*(cpdb|humanbase|string).*([0-9].[0-9]).dat/network-modules-\1-\2-\3.dat/g' labelled_network_modules_file > labelled_network_module_file_simplified.dat
process simplify_network_name_column {
    input:
    path labelled_network_modules_file
    val sed_simplify_script

    output:
    path 'labelled_network_module_file_simplified.dat'

    script:
    """
    $sed_simplify_script $labelled_network_modules_file > labelled_network_module_file_simplified_1.dat
    awk '!(\$2 = "")' labelled_network_module_file_simplified_1.dat > labelled_network_module_file_simplified.dat
    """
}

// enforce single tab delim and output with informative file name inc network and method
// change : to _ for compatibility with downstrean scripts for cut off determination
//  [network-modules-K1-cpdb-0.0.dat:118 -> network-modules-K1-cpdb-0.0.dat_118]
process finalise_preproc {
    input:
    path labelled_network_module_file_simplified
    val composite_modules_filename

    output:
    path "*.dat"
    // path 'module_file_finalised_preproc.dat'

    script:
    """
    sed 's/:/\\_/' $labelled_network_module_file_simplified > labelled_network_module_file_simplified1
    tr -s '[:blank:]' < labelled_network_module_file_simplified1 > labelled_network_module_file_simplified_tr
    tr ' ' '\t' < labelled_network_module_file_simplified_tr > $composite_modules_filename
    """
}

workflow PASCAL_GWAS {

    main:
    nf_network_modules_dir      =   "${params.nf_network_modules_dir}"
    nf_module_detect            =   "${params.nf_module_detect}"
    nf_network_type             =   "${params.nf_network_type }"
    sed_simplify_script         =   "${params.sed_simplify_script}"
    nf_network_modules_file     =   file(nf_network_modules_dir)
    composite_modules_filename  =   nf_network_modules_file.baseName + "_network_modules_" + nf_module_detect + "_" + nf_network_type + ".dat"
    
    ch_nets = Channel.fromPath( "$nf_network_modules_dir/$nf_module_detect/$nf_network_type/*.dat") 
    
    PREPROC_NET_MODULES( ch_nets, sed_simplify_script, composite_modules_filename )
        .collectFile(   name: composite_modules_filename, 
                        storeDir: "${workflow.launchDir}/output"
        )
        .subscribe  { str -> println "PREPROC_NET_MODULES: ${str}" }
    
}

workflow PREPROC_NET_MODULES {
    take:
    ch_nets
    sed_simplify_script
    composite_modules_filename

    main:
    ch_labelled     = label_network_module_files( ch_nets ) 
    ch_simplified   = simplify_network_name_column( ch_labelled, sed_simplify_script )
    ch_preproc      = finalise_preproc( ch_simplified, composite_modules_filename )
        
    emit:
    ch_preproc
        
}

// workflow CONVERT_NETWORK_TO_ENSG {
    
//     take:
//     ch_preproc
//     nf_network_type
    
//     main:
//     def nid_map = [ "cpdb": "ensembl_gene", "string": "ensembl_protein", "humanbase": "entrez_geneid"]
        
//     emit:
//     ch_net_ensg
        
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
// def validateInputSamplesheet(input) {
//     def (metas, fastqs) = input[1..2]

//     // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
//     def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
//     if (!endedness_ok) {
//         error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
//     }

//     return [ metas[0], fastqs ]
// }
// //
// // Generate methods description for MultiQC
// //
// def toolCitationText() {
//     // TODO nf-core: Optionally add in-text citation tools to this list.
//     // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
//     // Uncomment function in methodsDescriptionText to render in MultiQC report
//     def citation_text = [
//             "Tools used in the workflow included:",
//             "FastQC (Andrews 2010),",
            
//             "."
//         ].join(' ').trim()

//     return citation_text
// }

// def toolBibliographyText() {
//     // TODO nf-core: Optionally add bibliographic entries to this list.
//     // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
//     // Uncomment function in methodsDescriptionText to render in MultiQC report
//     def reference_text = [
//             "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            
//         ].join(' ').trim()

//     return reference_text
// }

// def methodsDescriptionText(mqc_methods_yaml) {
//     // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
//     def meta = [:]
//     meta.workflow = workflow.toMap()
//     meta["manifest_map"] = workflow.manifest.toMap()

//     // Pipeline DOI
//     if (meta.manifest_map.doi) {
//         // Using a loop to handle multiple DOIs
//         // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
//         // Removing ` ` since the manifest.doi is a string and not a proper list
//         def temp_doi_ref = ""
//         def manifest_doi = meta.manifest_map.doi.tokenize(",")
//         manifest_doi.each { doi_ref ->
//             temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
//         }
//         meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
//     } else meta["doi_text"] = ""
//     meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

//     // Tool references
//     meta["tool_citations"] = ""
//     meta["tool_bibliography"] = ""

//     // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
//     // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
//     // meta["tool_bibliography"] = toolBibliographyText()


//     def methods_text = mqc_methods_yaml.text

//     def engine =  new groovy.text.SimpleTemplateEngine()
//     def description_html = engine.createTemplate(methods_text).make(meta)

//     return description_html.toString()
// }

