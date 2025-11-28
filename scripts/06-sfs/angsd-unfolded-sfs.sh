#!/usr/bin/env bash
#SBATCH --job-name=sfs-ref-matched
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -a 1-3
#SBATCH --time=24:00:00
#SBATCH --mem=100GB
#SBATCH -o %x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# Estimating SFS from reference-matched BAM files. However, we're specifying
# the use of 'sites', as these are positions that we know have adequate coverage
# PER-SAMPLE, and have been properly filtered for genotype likelihood (AB). I know
# ANGSD works off genotype likelihoods (and is designed for low-cov data), but
# the goal here is to use a consistent set of sites for all analyses.
# These sites have been used in:
# - SFS estimation (here)
# - Per-species Tajima's D estimation
#
# A similar set of SNPs (without AB filtering) were used in MSMC2 analyses.

DIR="/home/${USER}/al-biohub/popgen"
BAMDIR="${DIR}/alignments-matched-updated/alignments"
VCFDIR="${DIR}/variants/matched"
OUT="${DIR}/sfs"

mkdir -p "${OUT}"
cd "${OUT}" || exit 1

BL=$(find "${BAMDIR}" -name "*bamlist" | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${BL%.*}")
SITES="${VCFDIR}/variants-${BN}/${BN}.sites"

case "${BN}" in
    "HCU")
        ANC="${DIR}/anc-genome/MGR-835-hcu_ref-anc.fa"
        ;;
    "HEL")
        ANC="${DIR}/anc-genome/MGR-835-hel_ref-anc.fa"
        ;;
    "HMA")
        ANC="${DIR}/anc-genome/MGR-835-hma_ref-anc.fa"
        ;;
esac

angsd \
    -sites "${SITES}" \
    -bam "${BL}" \
    -out "${BN}" \
    -anc "${ANC}" \
    -GL 1 \
    -doSaf 1 \
    -doCounts 1 \
    -remove_bads 1 \
    -uniqueOnly 1 \
    -only_proper_pairs 1 \
    -nThreads "${SLURM_CPUS_PER_TASK}"

# Create 1D-SFS for each species
realSFS -cores "${SLURM_CPUS_PER_TASK}" -bootstrap 100 ${BN}.saf.idx  > ${BN}.sfs