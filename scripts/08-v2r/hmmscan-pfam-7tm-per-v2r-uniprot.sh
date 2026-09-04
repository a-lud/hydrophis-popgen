#!/usr/bin/env bash
#SBATCH --job-name=7tm-per-v2r-uniprot-genes
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -a 1-349
#SBATCH --time=04:00:00
#SBATCH --mem=20GB
#SBATCH -o ./joblogs/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/$USER/hpcfs/analysis/popgen/results/v2r/domains-per-v2r"
OUTDIR="${DIR}/pfam-hits"
DB="/home/${USER}/al-biohub/database/Pfam-A.hmm"

mkdir -p "${OUTDIR}"

FAA=$(find "${DIR}/v2r-by-species" -name '*.faa' | tr '\n' ' ' | cut -d' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${FAA}" -v2r.faa)

echo "${BN}"
hmmscan \
    -o "${OUTDIR}/${BN}.out" \
    --tblout "${OUTDIR}/${BN}.tblout" \
    --domtblout "${OUTDIR}/${BN}.domtblout" \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    "${DB}" "${FAA}"