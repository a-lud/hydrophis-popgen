#!/usr/bin/env bash

# Merge per-sample all-site VCFs into a single VCF file for a given list of
# samples/species. Then fill tags and filter by F_MISSING <= 0.2. F_MISSING
# works here as we only have ONE group. It fails with two or more.

DIR="/home/${USER}/al/analyses/popgen"
DATA="${DIR}/data"
OUT="${DIR}/results/variants-hma/vcf-final/per-sample-filtered/analysis-ready"

mkdir -p "${OUT}"

cd "$(dirname "${OUT}")" || exit 1

# Per-species VCF files
for S in HCU HEL HMA; do
    echo "Processing: ${S}"
    GRP_FILE="${DATA}/group-${S}.tsv"

    bcftools merge --threads 4 ${S}*vcf.gz |
    bcftools annotate -x ^INFO/MQ,INFO/DP |
    bcftools plugin fill-tags -- \
        --samples-file "${GRP_FILE}" \
        -t 'DP:1=int(sum(FORMAT/DP)),F_MISSING,TYPE' |
    bcftools view -i 'F_MISSING<=0.2' -O z -o ${OUT}/${S}.vcf.gz --write-index --threads 2
done