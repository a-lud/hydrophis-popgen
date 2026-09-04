#!/usr/bin/env bash
#SBATCH --job-name=ips5-Pfam-validation
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -a 1-6
#SBATCH --time=04:00:00
#SBATCH --mem=40GB
#SBATCH -o ./joblogs/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/hpcfs/analysis/popgen"
DIRQUERY="${DIR}/results/v2r/tblastn-7tm-genome-search"
OUTDIR="${DIR}/results/v2r/curated-7tm"

mkdir -p "${OUTDIR}"

QUERY=$(find "${DIRQUERY}" -type f -name '*-7tm.candidates.fa' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")

IPS="/home/${USER}/al-biohub/database/interproscan-5.74-105.0"
${IPS}/interproscan.sh \
    -i "${QUERY}" \
    -t n \
    --appl Pfam \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    -b "${OUTDIR}/$(basename ${QUERY%-*})" \
    -dp