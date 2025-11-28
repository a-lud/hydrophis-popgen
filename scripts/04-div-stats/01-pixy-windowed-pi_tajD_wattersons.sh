#!/usr/bin/env bash
#SBATCH --job-name=pixy-pi
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 18
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

VCF=$(find "${DIR}" -name '*.vcf.gz' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${VCF%%.*}")

GROUP="/home/${USER}/hpcfs/analysis/popgen/data/group-${BN}.tsv"

source "/home/${USER}/hpcfs/micromamba/etc/profile.d/micromamba.sh"
micromamba activate pixy_update

parallel -j 5 --joblog ${OUT}/pixy-parallel-${BN}-100kb.log pixy \
    --stats pi watterson_theta tajima_d \
    --vcf "${VCF}" \
    --populations "${GROUP}" \
    --chromosomes {} \
    --window_size 100000 \
    --n_cores 1 \
    --output_folder "${OUT}" \
    --output_prefix "${BN}-{}-100kb" ::: chr{1..15}

parallel -j 5 --joblog ${OUT}/pixy-parallel-${BN}-50kb.log pixy \
    --stats pi watterson_theta tajima_d \
    --vcf "${VCF}" \
    --populations "${GROUP}" \
    --chromosomes {} \
    --window_size 50000 \
    --n_cores 1 \
    --output_folder "${OUT}" \
    --output_prefix "${BN}-{}-50kb" ::: chr{1..15}

micromamba deactivate