#!/usr/bin/env bash

GFF_DIR="/home/${USER}/al-biohub/sequence-db/gffs"
OUT_DIR="/home/${USER}/hpcfs/analysis/popgen/results/v2r/intron-length-dist"

mkdir -p "${OUT_DIR}"

for i in notechis_scutatus.gff3 pseudonaja_textilis.gff3 protobothrops_mucrosquamatus.gff3; do
    BN=$(basename $i .gff3)
    echo "${BN}"

    agat_sp_filter_feature_from_keep_list \
        --gff ${GFF_DIR}/$i \
        --kl "${OUT_DIR}/${BN}-v2r-gene-ids.txt" \
        -o "${OUT_DIR}/${BN}-v2r-genes.gff3"

    agat_sp_add_introns \
        --gff "${OUT_DIR}/${BN}-v2r-genes.gff3" \
        -o "${OUT_DIR}/${BN}-v2r-genes-with-introns.gff3"
done

