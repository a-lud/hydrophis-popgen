#!/usr/bin/env bash

DIR="/home/${USER}/hpcfs/analysis/msmc"
NF="${DIR}/nf-pipelines"

BAM="${DIR}/bams/HEL"
COV="${DIR}/coverage/HEL-coverage.txt"
ASM="${DIR}/genome/Hydrophis_elegans.fa"

nextflow run "${NF}/main.nf" \
    --pipeline 'msmc2' \
    -profile 'conda,phoenix' \
    --partition 'icelake' \
    -work-dir "${DIR}/work-msmc2-hel" \
    -N 'alastair.ludington@adelaide.edu.au' \
    -resume \
    --outdir "${DIR}/results" \
    --out_prefix "hydrophis-elegans" \
    --genome "${ASM}" \
    --bams "${BAM}" \
    --coverage "${COV}" \
    --K 49 \
    --E 1 \
    --mapq 30 \
    --baseq 30 \
    --clock '2*2+25*1+2*3 1*2+25*1+1*2+1*3' \
    --msmctools "/hpcfs/users/${USER}/software/msmc-tools"

