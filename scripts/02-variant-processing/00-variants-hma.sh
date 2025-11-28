#!/usr/bin/env bash

## Location of where the pipeline is installed
PIPE="/home/${USER}/al/tools_AL/nf-pipelines"
DIR="/home/${USER}/al/analyses/popgen"

nextflow run "${PIPE}/main.nf" \
   --pipeline 'variant' \
    -profile 'conda,afw' \
    --partition 'afw' \
    -work-dir "${DIR}/work-variants-hma" \
    --outdir "${DIR}/results" \
    --out_prefix "variants-hma" \
    --sheet "${DIR}/scripts/samplesheet-variants-hma.csv" \
    --method 'joint' \
    --vcftype 'allsites' \
    --regions "${DIR}/data/regions" \
    --mapq 20 \
    --baseq 20 \
    --genoq 20 \
    --caller new \
    --filter_opt '--SnpGap 10' \
    --view_opt '--exclude-types indels' \
    --call_opt '--group-samples /home/${USER}/al/analyses/popgen/data/groups.tsv' \
    -resume

