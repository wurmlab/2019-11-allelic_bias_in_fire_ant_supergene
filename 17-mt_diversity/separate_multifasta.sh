#!/bin/bash

module load samtools
module load parallel
module load R

#Generate multiple fasta files
awk 'BEGIN {O="";} /^>/ { O=sprintf("tmp/%s.fa",substr($0,2));} {if(O!="") print >> O;}' input/multifasta.fa
#Rename the files to remove whitespaces
rename ' ' '_' tmp/*

#Replace whitespace in the headers of the fasta files
ls tmp/*fa | parallel 'sed -i "s/ /_/" {}'  :::

#Remove all fasta files except the ones from Argentina (AR) and the S.invicta MT reference
find tmp/* | egrep -v 'AR43-|AR3-|AR28-|AR118-|AR117-|AR112-|AR111-|AR104-' | grep -v 'MTinv_' | xargs rm

#Extract the sequence for COX1 in each fasta file
SAMPLES=($(ls tmp/*fa | cut -d '/' -f 2 | sed 's/\.fa$//'))
parallel 'samtools faidx tmp/{}.fa {}:0-1528' ::: ${SAMPLES} > tmp/all_cox1.fa

echo 'library(Biostrings)
      dna = readDNAStringSet("tmp/all_cox1.fa")
      cox1_dists <- stringDist(dna, method="hamming")
      save(cox1_dists, file = "tmp/cox1_dists.RData")' | Rscript -
