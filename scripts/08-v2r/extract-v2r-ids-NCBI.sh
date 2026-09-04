#!/usr/bin/env bash

GFF_DIR="/home/${USER}/al-biohub/sequence-db/gffs"
OUT_DIR="/home/${USER}/hpcfs/analysis/popgen/results/v2r/intron-length-dist"

mkdir -p "${OUT_DIR}"

for i in notechis_scutatus.gff3 pseudonaja_textilis.gff3 protobothrops_mucrosquamatus.gff3; do
    BN=$(basename $i .gff3)
    echo "${BN}"

    grep "vomeronasal type-2" ${GFF_DIR}/$i | 
    grep "mRNA" | 
    cut -f 9 | 
    sed 's/.*;gene=\(.*\);model_evidence.*/\1/' | 
    grep -v "ID=" | 
    sort | 
    uniq | 
    sed 's/\(.*\)/gene-\1/' > "${OUT_DIR}/${BN}-v2r-gene-ids.txt"
done

