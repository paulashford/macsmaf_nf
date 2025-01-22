#!/usr/bin/env bash
# run_workflows.sh

export NF_WORKFLOW_DIR=/Users/ash/git/macsmaf/macsr_nf
export NF_CONFIG=${NF_WORKFLOW_DIR}/macsr_nf_dev/conf/base.config
# note this is the macsr_nf_dev directory and subworkflow etc scripts should use relative paths
export NXF_SCRIPT_DIR=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev

# if need rds to tsv 
# example: 
# export input_rds=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/output/enrichment/K1_cpdb_0.9_enrichment_with_metrics.rds
# export output_tsv=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/output/enrichment/K1_cpdb_0.9_enrichment_with_metrics.tsv
# Rscript ${NF_WORKFLOW_DIR}/macsr_nf_dev/script/rds_to_tsv.r $input_rds $output_tsv

# Run specific workflows /test workflows

nextflow run ${NF_WORKFLOW_DIR}/macsr_nf_dev/workflows/macs_nf_dev.nf -c "${NF_CONFIG}"
# nextflow run ${NF_WORKFLOW_DIR}/macsr_nf_dev/workflows/macs_nf_dev.nf -c "${NF_CONFIG}" -resume
#nextflow run ${NF_WORKFLOW_DIR}/macsr_nf_dev/workflows/macs_nf_dev_test_go_slim.nf -c "${NF_CONFIG}"
# nextflow run ${NF_WORKFLOW_DIR}/macsr_nf_dev/workflows/macs_nf_dev_test_go_slim_rank.nf -c "${NF_CONFIG}"

