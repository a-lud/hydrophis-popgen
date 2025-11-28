#!/usr/bin/env bash
#SBATCH --job-name=dsuite-dtrios-no_tree
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH -o %x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/al-biohub/popgen"
OUT="${DIR}/dsuite-dtrios"
VCF="${DIR}/variants/variants-hma/HCU-HEL-HMA-MGR.vcf.gz"

mkdir -p "${OUT}"
cd "${OUT}" || exit 1

/hpcfs/users/${USER}/software/Dsuite-master/Build/Dsuite Dtrios \
    -o 'dsuite_' \
    -n 'no_tree' \
    --no-combine \
    --no-f4-ratio \
    --JKnum 100 \
    "${VCF}" \
    "${DIR}/data/dsuite-sets.txt"