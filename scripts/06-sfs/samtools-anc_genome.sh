#!/usr/bin/env bash
#SBATCH --job-name=anc_genome
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -a 1-3
#SBATCH --time=04:00:00
#SBATCH --mem=20GB
#SBATCH -o %x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/al-biohub/popgen"
BAMDIR="${DIR}/alignments-matched-updated/alignments"
OUT="${DIR}/anc-genome"

mkdir -p "${OUT}"
cd "${OUT}" || exit 1

BAM=$(find "${BAMDIR}" -name "MGR*.bam" | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${BAM%.*}")

samtools consensus \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -a \
    --min-MQ 20 \
    --min-BQ 20 \
    --show-ins no \
    --show-del yes \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -o "${BN}-anc.fa" \
    "${BAM}"