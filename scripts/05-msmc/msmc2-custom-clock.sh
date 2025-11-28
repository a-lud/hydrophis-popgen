#!/usr/bin/env bash
#SBATCH --job-name=msmc
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH -a 1-15
#SBATCH --time=24:00:00
#SBATCH --mem=30GB
#SBATCH -o %x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/${USER}/al-biohub/popgen"
OUT="${DIR}/msmc"

mkdir -p "${OUT}"
cd "${OUT}" || exit 1

SAMPLE=$(head -n 15 "${DIR}/data/groups.tsv" | cut -f 1 | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
SPECIES=${SAMPLE%-*}

CLOCK_1="4*1+25*2+1*4+1*6"

case $SPECIES in
    HMA)
        MULTHEPSEP_DIR="/home/${USER}/hpcfs/analysis/msmc/results/hydrophis-major/multHepSep"
        msmc2 \
            -t "${SLURM_CPUS_PER_TASK}" \
            -I 0,1 \
            -p "${CLOCK_1}" \
            -i 100 \
            -o "${OUT}/${SAMPLE}-clock1-msmc.txt" ${MULTHEPSEP_DIR}/${SAMPLE}.chr{1..15}.txt
        ;;
    HCU)
        MULTHEPSEP_DIR="/home/${USER}/hpcfs/analysis/msmc/results/hydrophis-curtus/multHepSep"
        msmc2 \
            -t "${SLURM_CPUS_PER_TASK}" \
            -I 0,1 \
            -p "${CLOCK_1}" \
            -i 100 \
            -o "${OUT}/${SAMPLE}-clock1-msmc.txt" ${MULTHEPSEP_DIR}/${SAMPLE}.chr{1..16}.txt
        ;;
    HEL)
        MULTHEPSEP_DIR="/home/${USER}/hpcfs/analysis/msmc/results/hydrophis-elegans/multHepSep"
        msmc2 \
            -t "${SLURM_CPUS_PER_TASK}" \
            -I 0,1 \
            -p "${CLOCK_1}" \
            -i 100 \
            -o "${OUT}/${SAMPLE}-clock1-msmc.txt" ${MULTHEPSEP_DIR}/${SAMPLE}.chr{1..15}.txt
        ;;
esac

