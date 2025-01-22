// Subworkflow for GOSlim analysis of top 25% ranked GO terms for each modules gene set
nextflow.enable.dsl=2

process RUN_GO_SLIM {
    publishDir "${params.nf_enrichment_dir}/go_slim", mode: 'copy'

    input:
    tuple val(method), val(db), val(cutoff), path(enrichment_results_rds)
    val min_perc_rank

    output:
    tuple val(method), val(db), val(cutoff), path("go_slim_results_${method}_${db}_${cutoff}.txt"), emit: go_slim_results

    script:
    """
    Rscript ${params.root_proj_dir}/subworkflows/go_slim/go_slim.r \\
        --input '${enrichment_results_rds}' \\
        --min_perc_rank '${min_perc_rank}' \\
        --output 'go_slim_results_${method}_${db}_${cutoff}.txt'
    """
}

process MAP_TO_SLIM {
    publishDir "${params.nf_enrichment_dir}/go_slim", mode: 'copy'

    input:
    tuple val(method), val(db), val(cutoff), path(association_file)
    path go_obo
    path go_slim_obo

    output:
    tuple val(method), val(db), val(cutoff), path("*_goslim.txt"), emit: mapped_slim_results

    script:
    """
    python3 ${params.goatools_script_path}/map_to_slim.py \\
        --association_file=${association_file} \\
        ${go_obo} \\
        ${go_slim_obo} > ${association_file.baseName}_goslim.txt
    """
}

process clean_slim_results {
    publishDir "${params.nf_enrichment_dir}/go_slim", mode: 'copy'

    input:
    tuple val(method), val(db), val(cutoff), path(goslim_file)

    output:
    tuple val(method), val(db), val(cutoff), path("${method}_${db}_${cutoff}_goslim_clean.txt"), emit: cleaned_results

    script:
    """
    # Remove first 5 lines and add header
    tail -n +6 ${goslim_file} > temp.txt
    echo -e "func_module_number\tgo_slim_pir" > "${method}_${db}_${cutoff}_goslim_clean.txt"
    cat temp.txt >> "${method}_${db}_${cutoff}_goslim_clean.txt"
    """
}

workflow GO_SLIM {
    take:
    enrichment_results  // tuple(method, db, cutoff, file)
    min_perc_rank

    main:
    RUN_GO_SLIM(enrichment_results, min_perc_rank)
    
    MAP_TO_SLIM(
        RUN_GO_SLIM.out.go_slim_results,
        params.go_obo_path,
        params.go_slim_obo_path
    )

    clean_slim_results(MAP_TO_SLIM.out.mapped_slim_results)

    emit:
    go_slim_results = RUN_GO_SLIM.out.go_slim_results
    mapped_slim_results = clean_slim_results.out.cleaned_results
}