#!/usr/bin/env bash

# Merge per-sample all-site VCFs for 2 species, re-calculate DP and filter
# on species-group F_MISSING <= 0.2. Had to write a custom script to do this
# as BCFtools +fill-tags is not working correctly, as it assigns the cohort
# F_MISSING to species groups.

DIR="/home/${USER}/al/analyses/popgen"
OUT="${DIR}/results/variants-hma/vcf-final/per-sample-filtered/analysis-ready"

mkdir -p "${OUT}"

cd "$(dirname "${OUT}")" || exit 1

# ---------------------------------------------------------------------------- #
# Pairwise VCF files

# HCU - HEL
echo "Merging HCU and HEL..."
bcftools merge --threads 4 HCU-*vcf.gz HEL-*.vcf.gz |
bcftools annotate -x ^INFO/DP,INFO/MQ |
bcftools plugin fill-tags -- -t 'DP:1=int(sum(FORMAT/DP))' |
./f_missing-cyvcf2.py -s1 HCU -s2 HEL -o "${OUT}/HCU-HEL.vcf.gz" - &


# HCU - HMA
echo "Merging HCU and HMA..."
bcftools merge --threads 4 HCU-*vcf.gz HMA-*.vcf.gz |
bcftools annotate -x ^INFO/DP,INFO/MQ |
bcftools plugin fill-tags -- -t 'DP:1=int(sum(FORMAT/DP))' |
./f_missing-cyvcf2.py -s1 HCU -s2 HMA -o "${OUT}/HCU-HMA.vcf.gz" - &


# HEL - HMA
echo "Merging HEL and HMA..."
bcftools merge --threads 4 HEL-*vcf.gz HMA-*.vcf.gz |
bcftools annotate -x ^INFO/DP,INFO/MQ |
bcftools plugin fill-tags -- -t 'DP:1=int(sum(FORMAT/DP))' |
./f_missing-cyvcf2.py -s1 HEL -s2 HMA -o "${OUT}/HEL-HMA.vcf.gz" - &

wait

echo "Done!"