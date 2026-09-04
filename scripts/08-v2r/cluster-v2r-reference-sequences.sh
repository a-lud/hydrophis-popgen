#!/usr/bin/env bash

# Collapse similar V2R proteins
# Useful information: https://mmseqs.com/latest/userguide.pdf
#   - FAQ secion: How to set the right alignment coverage to cluster

DIR="/home/${USER}/hpcfs/analysis/popgen/data/v2r-proteins"
PROT="${DIR}/v2r-clean.faa"

mkdir -p "${DIR}/tmp"

mmseqs easy-cluster \
    ${PROT} ${DIR}/v2r-clean-clustered ${DIR}/tmp \
    -c 0.9 \
    --cov-mode 1 \
    --cluster-mode 2 \
    --threads 16 \
    --local-tmp ${DIR}/tmp

