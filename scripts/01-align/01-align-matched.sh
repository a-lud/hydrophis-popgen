#!/usr/bin/env bash

NF="/home/${USER}/al/tools_AL/nf-pipelines"
DIR="/home/${USER}/al/analyses/popgen"
SEQ="${DIR}/data/reads"

nextflow run "${NF}/main.nf" \
    --pipeline 'alignment' \
    -profile 'conda,afw' \
    --partition 'afw' \
    -work-dir "${DIR}/work-matched-reference" \
    -N 'alastair.ludington@adelaide.edu.au' \
    -resume \
    --outdir "${DIR}/results" \
    --out_prefix "alignments-matched" \
    --seqdir "${SEQ}" \
    --sheet "matched-genome.csv" \
    --platform 'illumina' \
    --aligner 'bwa2' \
    --mapq 20

