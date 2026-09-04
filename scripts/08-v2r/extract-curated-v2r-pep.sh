#!/usr/bin/env bash

cd "$HOME/hpcfs/analysis/popgen/results/v2r/curtating-miniprot-alignments" || exit 1

for i in *curated.gff3; do
    BN=$(basename ${i%-miniprot-curated.gff3})
    agat_sp_extract_sequences \
        --clean_final_stop \
        -f "../../../data/genomes-miniprot-v2r/${BN}.fa" \
        -g "$i" \
        -t cds \
        -p \
        -o "${i%.gff3}.pep"
done

