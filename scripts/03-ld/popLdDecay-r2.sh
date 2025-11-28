#!/usr/bin/env bash
#SBATCH --job-name=popLdDecay-r2
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -a 1-3
#SBATCH --time=02:00:00
#SBATCH --mem=5GB
#SBATCH -o %x_%a_%A_%j.log
#SBATCH --mail-type=FAIL

DIR="/home/${USER}/al-biohub/popgen"
VCFDIR="${DIR}/variants/variants-hma"
OUT="${DIR}/linkage"
POPLD="/home/${USER}/hpcfs/software/PopLDdecay/bin"

mkdir -p "${OUT}"; cd "${OUT}" || exit 1

# Per species population file
VCF=$(find "${VCFDIR}" -name '*vcf.gz' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${VCF%%.*}")

"${POPLD}/PopLDdecay" -InVCF "${VCF}" -OutStat "${BN}.linkage.stat.gz"
"${POPLD}/Plot_OnePop.pl" -inFile "${BN}.linkage.stat.gz" -o "${BN}"
