#!/usr/bin/env bash
#SBATCH --job-name=genomics_general-fd
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=04:00:00
#SBATCH --mem=16GB
#SBATCH -o %x_%j.log
#SBATCH --mail-type=FAIL

# The speciation order between the three species is not clear. fd measures
# introgression from P3 into P2. Therefore, we try all permutations of the
# populations to ensure we're not missing anything. Output
# files are named using the following format statistics_P1-P2-P3-O.csv.

source "/home/${USER}/hpcfs/micromamba/etc/profile.d/micromamba.sh"
micromamba activate gg

GG="/home/${USER}/hpcfs/software/genomics_general"

DIR="/home/${USER}/al-biohub/popgen"
VCF="${DIR}/variants/variants-hma/HCU-HEL-HMA-MGR-biallelic.vcf.gz"
OUT="${DIR}/fd"

cd "${OUT}" || exit 1

# Convert VCF to '.geno' object
python "${GG}/VCF_processing/parseVCF.py" -i "${VCF}" | bgzip > 'HCU-HEL-HMA-MGR-biallelic.geno.gz'

# fd statistics
python "${GG}/ABBABABAwindows.py" \
     -w 100000 \
     -m 1 \
     -P1 'HCU' -P2 'HEL' -P3 'HMA' -O 'MGR' \
     --popsFile "${DIR}/data/groups.tsv" \
     -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
     -o "fd_HCU-HEL-HMA-MGR.csv" \
     -f 'phased' \
     -T "${SLURM_CPUS_PER_TASK}" \
     --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
     -w 50000 \
     -m 1 \
     -P1 'HCU' -P2 'HEL' -P3 'HMA' -O 'MGR' \
     --popsFile "${DIR}/data/groups.tsv" \
     -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
     -o "fd_HCU-HEL-HMA-MGR-50kb.csv" \
     -f 'phased' \
     -T "${SLURM_CPUS_PER_TASK}" \
     --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
     -w 100000 \
     -m 1 \
     -P1 'HEL' -P2 'HCU' -P3 'HMA' -O 'MGR' \
     --popsFile "${DIR}/data/groups.tsv" \
     -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
     -o "fd_HEL-HCU-HMA-MGR.csv" \
     -f 'phased' \
     -T "${SLURM_CPUS_PER_TASK}" \
     --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
     -w 50000 \
     -m 1 \
     -P1 'HEL' -P2 'HCU' -P3 'HMA' -O 'MGR' \
     --popsFile "${DIR}/data/groups.tsv" \
     -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
     -o "fd_HEL-HCU-HMA-MGR-50kb.csv" \
     -f 'phased' \
     -T "${SLURM_CPUS_PER_TASK}" \
     --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 100000 \
    -m 1 \
    -P1 'HMA' -P2 'HEL' -P3 'HCU' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HMA-HEL-HCU-MGR.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 50000 \
    -m 1 \
    -P1 'HMA' -P2 'HEL' -P3 'HCU' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HMA-HEL-HCU-MGR-50kb.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 100000 \
    -m 1 \
    -P1 'HEL' -P2 'HMA' -P3 'HCU' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HEL-HMA-HCU-MGR.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 50000 \
    -m 1 \
    -P1 'HEL' -P2 'HMA' -P3 'HCU' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HEL-HMA-HCU-MGR-50kb.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 100000 \
    -m 1 \
    -P1 'HCU' -P2 'HMA' -P3 'HEL' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HCU-HMA-HEL-MGR.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 50000 \
    -m 1 \
    -P1 'HCU' -P2 'HMA' -P3 'HEL' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HCU-HMA-HEL-MGR-50kb.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 100000 \
    -m 1 \
    -P1 'HMA' -P2 'HCU' -P3 'HEL' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HMA-HCU-HEL-MGR.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

python "${GG}/ABBABABAwindows.py" \
    -w 50000 \
    -m 1 \
    -P1 'HMA' -P2 'HCU' -P3 'HEL' -O 'MGR' \
    --popsFile "${DIR}/data/groups.tsv" \
    -g 'HCU-HEL-HMA-MGR-biallelic.geno.gz' \
    -o "fd_HMA-HCU-HEL-MGR-50kb.csv" \
    -f 'phased' \
    -T "${SLURM_CPUS_PER_TASK}" \
    --writeFailedWindows

micromamba deactivate