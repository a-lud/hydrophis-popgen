#!/usr/bin/env bash

# Three species VCF file.
# Merge multi-sample per-species VCFs then filter on F_MISSING_??? <= 0.2
# in each species. Adjust depth statistic after.

bcftools merge --threads 4 H*gz MGR*vcf.gz |
bcftools view -i 'F_MISSING_HCU <= 0.2 && F_MISSING_HEL <= 0.2 && F_MISSING_HMA <= 0.2 && F_MISSING_MGR' |
bcftools plugin fill-tags -O z -o HCU-HEL-HMA.vcf.gz - -- -t 'DP:1=int(sum(FORMAT/DP))'