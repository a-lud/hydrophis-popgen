#!/usr/bin/env bash

DIR="/home/${USER}/hpcfs/analysis/popgen/results/v2r/curtating-miniprot-alignments"
cd $DIR || exit 1

for i in *-miniprot-curated.gff3; do 
    awk '$3 == "gene"' $i | 
    bedtools merge -i - -c 9 -o collapse > ${i%.gff3}.bed
done
