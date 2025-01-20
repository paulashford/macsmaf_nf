// Subworkflow for GOSlim analysis of top 25% ranked GO terms for each modules gene set
nextflow.enable.dsl=2

process RUN_GO_SLIM {
    publishDir "${params.nf_enrichment_dir}/go_slim", mode: 'copy'

    input:
    path enrichment_results_rds
    val min_perc_rank

    output:
    path "go_slim_results_*.txt", emit: go_slim_results

    script:
    """
    Rscript ${params.root_proj_dir}/subworkflows/go_slim/go_slim.r \\
        --input '${enrichment_results_rds}' \\
        --min_perc_rank '${min_perc_rank}' \\
        --output 'go_slim_results_${enrichment_results_rds.simpleName}.txt'
    """
}

process MAP_TO_SLIM {
    publishDir "${params.nf_enrichment_dir}/go_slim", mode: 'copy'

    input:
    path association_file
    path go_obo
    path go_slim_obo

    output:
    path "*_goslim.txt", emit: mapped_slim_results

    script:
    """
    python3 ${params.goatools_script_path}/map_to_slim.py \\
        --association_file=${association_file} \\
        ${go_obo} \\
        ${go_slim_obo} > ${association_file.baseName}_goslim.txt
    """
}

workflow GO_SLIM {
    take:
    enrichment_results_rds
    min_perc_rank

    main:
    RUN_GO_SLIM(enrichment_results_rds, min_perc_rank)
    
    MAP_TO_SLIM(
        RUN_GO_SLIM.out.go_slim_results,
        params.go_obo_path,
        params.go_slim_obo_path
    )

    emit:
    go_slim_results = RUN_GO_SLIM.out.go_slim_results
    mapped_slim_results = MAP_TO_SLIM.out.mapped_slim_results
}