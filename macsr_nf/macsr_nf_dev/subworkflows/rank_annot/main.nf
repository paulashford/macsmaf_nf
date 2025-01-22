#!/usr/bin/env nextflow

process RANK_ANNOTATE {
    publishDir "${params.nf_enrichment_dir}/rank_annotated", mode: 'copy'
    
    input:
    tuple val(method), val(db), val(cutoff), path(enrichment_file)
    tuple val(method2), val(db2), val(cutoff2), path(goslim_file)
    tuple val(method3), val(db3), val(cutoff3), path(modules_file)
    val(max_term_size)

    output:
    tuple val(method), val(db), val(cutoff), path("*_ranked_annotated.tsv"), emit: ranked_results
    tuple val(method), val(db), val(cutoff), path("*_topk_aggregated.tsv"), emit: aggregated_results
    tuple val(method), val(db), val(cutoff), path("*_rank1_aggregated.tsv"), emit: rank1_results
    tuple val(method), val(db), val(cutoff), path("*_rank_agg_final.tsv"), emit: final_results
    
    script:
    // Verify tuples match
    if (method != method2 || db != db2 || cutoff != cutoff2 ||
        method != method3 || db != db3 || cutoff != cutoff3) {
        error "Mismatched inputs: ${method}_${db}_${cutoff} vs ${method2}_${db2}_${cutoff2} vs ${method3}_${db3}_${cutoff3}"
    }
    """
    Rscript ${params.root_proj_dir}/subworkflows/rank_annot/rank_annot.r \
        --enrichment_file ${enrichment_file} \
        --goslim_file ${goslim_file} \
        --modules_file ${modules_file} \
        --output_file "${method}_${db}_${cutoff}_ranked_annotated.tsv" \
        --max_term_size ${max_term_size} \
        --perc_rank_cutoff 0.25
    """
}

workflow RANK_ANNOT {
    take:
    enrichment_results  // tuple(method, db, cutoff, file) with metrics
    goslim_results      // tuple(method, db, cutoff, file) after mapping to slim
    modules_results     // tuple(method, db, cutoff, file) with gene lists
    max_term_size      // value

    main:
    RANK_ANNOTATE(
        enrichment_results,
        goslim_results,
        modules_results,
        max_term_size
    )

    emit:
    ranked_results = RANK_ANNOTATE.out.ranked_results
    aggregated_results = RANK_ANNOTATE.out.aggregated_results
    rank1_results = RANK_ANNOTATE.out.rank1_results
    final_results = RANK_ANNOTATE.out.final_results
}
