#!/usr/bin/env bash

# By sample
#   - Filter by min 10x coverage max 2 * avg depth (sample specific)
#   - Filter by allelic balance (AB) 0.2 <= AB <= 0.8

DIR="/home/${USER}/al/analyses/popgen"
DATA="${DIR}/data"
COVDIR="${DIR}/results/variants-hma/coverage-stats"
OUT="${DIR}/results/variants-hma/vcf-final/per-sample-filtered"

mkdir -p "${OUT}"

cd "${OUT}" || exit 1

for SAMPLE in $(cat ${DATA}/*samples | sort | uniq); do
    # Get coverage values
    COV=$(tail -n+2 "${COVDIR}/${SAMPLE}-coverage.csv" | cut -d ',' -f 3)
    COVR=$(printf "%.0f" "${COV}")
    COV_2X=$(echo "${COVR} * 2" | bc)
    echo "${SAMPLE} AVG:${COVR} MAX:${COV_2X}"

    # Subset the join-called VCF for the sample then apply per-sample depth filters
    # as well as allelic-balance (AB) filters. AB filters ensure that HET calls are
    # more accurate. It ensures that sites with e.g. 1 ALT-read don't get called hets...
    bcftools view -O u --samples "${SAMPLE}" '../Hydrophis_major.sort.vcf.gz' |
    bcftools plugin setGT -- \
        -t q \
        -n './.' \
        -i "FMT/DP<10 || FMT/DP>${COV_2X}" |
    bcftools plugin setGT -O z -o "${SAMPLE}-filtered.vcf.gz" -- \
        -t q \
        -n './.' \
        -i "GT=\"het\" && (AD[0:1]/(AD[0:0]+AD[0:1]) < 0.2 || AD[0:1]/(AD[0:0]+AD[0:1]) > 0.8)"
done