#!/usr/bin/env bash
#SBATCH --job-name=tblastn-7tm-snake-genomes
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -a 1-6
#SBATCH --time=24:00:00
#SBATCH --mem=40GB
#SBATCH -o ./joblogs/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/hpcfs/analysis/popgen"
OUTDIR="${DIR}/results/v2r/tblastn-7tm-genome-search"
QUERY="${OUTDIR}/HMA.7tm.cdhit.fa"
GENOMES="${DIR}/data/genomes-miniprot-v2r"

mkdir -p "${OUTDIR}"

cd "${OUTDIR}" || exit 1

REF=$(find "${GENOMES}" -type f -name '*.fa' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
if [[ ! -f "${REF}.ndb" ]]; then
    makeblastdb -input_type fasta -in "${REF}" -dbtype nucl -out "${REF}" -parse_seqids
fi

tblastn \
    -query "${QUERY}" \
    -db "${REF}" \
    -out "$(basename ${REF%.fa}).7tm.tblastn.outfmt6" \
    -outfmt 6 \
    -num_threads "${SLURM_CPUS_PER_TASK}"