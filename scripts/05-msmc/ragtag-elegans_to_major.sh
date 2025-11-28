#!/usr/bin/env bash
#SBATCH --job-name=ragtag-elegans_to_major
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --ntasks-per-core=2
#SBATCH --time=08:00:00
#SBATCH --mem=80GB
#SBATCH -o /home/a1645424/hpcfs/analysis/popgen/scripts/joblogs/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/hpcfs/analysis/popgen"
QRY="${DIR}/data/genomes/hydrophis_elegans.fa"
REF="${DIR}/data/genomes/hydrophis_major.fa"
OUT="${DIR}/results/ragtag-scaffold"

mkdir -p "${OUT}"

source "/hpcfs/users/${USER}/micromamba/etc/profile.d/micromamba.sh"
micromamba activate ragtag

# Scaffold H. elegans against H. major
ragtag.py scaffold \
    --remove-small \
    -o "${OUT}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    "${REF}" "${QRY}"

micromamba deactivate
