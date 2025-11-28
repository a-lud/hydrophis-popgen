#!/usr/bin/env bash
#SBATCH --job-name=pixy-dxy_fst
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --ntasks-per-core=1
#SBATCH -a 1-3
#SBATCH --time=08:00:00
#SBATCH --mem=40GB
#SBATCH -o %x_%a_%A_%j.log
#SBATCH --mail-type=FAIL

# Within species diversity statistics.
# Using two window sizes: 100kb and 50kb. The latter should be fine, as LD
# curves appear to be flattening at this point, but trying 100kb to compare.

DIR="/home/${USER}/al-biohub/popgen/variants/variants-hma"
OUT="/home/${USER}/al-biohub/popgen/div-stats"

# VCF=$(find "${DIR}" -name '*.vcf.gz' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
VCF=$(echo "${DIR}/HCU-HMA.vcf.gz ${DIR}/HEL-HMA.vcf.gz ${DIR}/HCU-HEL.vcf.gz" | cut -d' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${VCF%%.*}")

GROUP="/home/${USER}/al-biohub/popgen/data/group-${BN}.tsv"

source "/home/${USER}/hpcfs/micromamba/etc/profile.d/micromamba.sh"
micromamba activate pixy_update

pixy \
    --stats dxy fst \
    --vcf "${VCF}" \
    --populations "${GROUP}" \
    --chromosomes 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15' \
    --window_size 100000 \
    --n_cores 15 \
    --output_folder "${OUT}" \
    --output_prefix "${BN}-100kb"

pixy \
    --stats dxy fst\
    --vcf "${VCF}" \
    --populations "${GROUP}" \
    --chromosomes 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15' \
    --window_size 50000 \
    --n_cores 15 \
    --output_folder "${OUT}" \
    --output_prefix "${BN}-50kb"

micromamba deactivate