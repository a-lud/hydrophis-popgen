#!/usr/bin/env bash
#SBATCH --job-name=v2r-representative-miniprot
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -a 1-6
#SBATCH --time=02:00:00
#SBATCH --mem=40GB
#SBATCH -o ./joblogs/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/hpcfs/analysis/popgen"
PROT="${DIR}/data/v2r-proteins/v2r-clean-clustered_rep_seq.fasta"
OUT="${DIR}/results/v2r/protein-alignments-representative"
GENOMES="${DIR}/data/genomes-miniprot-v2r"

mkdir -p "${OUT}"

REF=$(find "${GENOMES}" -type f -name '*.fa' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${REF%.*}")

# Align UniProt V2R sequences to each genome
miniprot \
    -G 45000 \
    -t "${SLURM_CPUS_PER_TASK}" \
    --gff \
    -P "${BN}_" \
    --outc=0.5 \
    "${REF}" \
    "${PROT}" > "${OUT}/${BN}-miniprot.gff3"
