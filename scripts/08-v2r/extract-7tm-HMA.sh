#!/usr/bin/env bash

# Grab 7TM domains from C. adamanteus genome and save to BED file
grep -v '^#' HMA-merged.domtblout | 
grep "7tm" | 
sed 's/ \+/\t/g' | 
cut -f 1-4,12,16-19 | 
awk -F'\t' 'BEGIN{OFS="\t"} $5 <= 1e-5 { $10 = $7 - $6 + 1; $11 = $9 - $8 + 1; if($10 == $3) print}' | 
cut -f 4,8,9 | 
bedtools getfasta -fi ../curtating-miniprot-alignments/HMA-miniprot-curated-merged.pep -fo ../tblastn-7tm-genome-search/HMA.7tm.fa -bed -

cd ../.tblastn-7tm-genome-search || exit 1

cd-hit -i HMA.7tm.fa -o HMA.7tm.cdhit.fa