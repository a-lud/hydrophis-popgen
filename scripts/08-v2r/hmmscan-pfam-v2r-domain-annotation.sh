#!/usr/bin/env bash
#SBATCH --job-name=hmmscan-v2r-domain-annotation
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -a 1-6
#SBATCH --time=02:00:00
#SBATCH --mem=20GB
#SBATCH -o ./joblogs/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/hpcfs/analysis/popgen/results/v2r/evm-v2r"
PEP_DIR="/home/${USER}/hpcfs/analysis/popgen/results/v2r/curtating-miniprot-alignments"
OUTDIR="/home/${USER}/hpcfs/analysis/popgen/results/v2r/pfam-curated-v2r-proteins"
DB="/home/${USER}/al-biohub/database/Pfam-A.hmm"

FAA=$(find "${PEP_DIR}" -name '*merged.pep' | tr '\n' ' ' | cut -d' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${FAA}" -miniprot-curated-merged.pep)

mkdir -p "${OUTDIR}"

echo "${BN}"
hmmscan \
    -o "${OUTDIR}/${BN}-merged.out" \
    --tblout "${OUTDIR}/${BN}-merged.tblout" \
    --domtblout "${OUTDIR}/${BN}-merged.domtblout" \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    "${DB}" "${FAA}"
