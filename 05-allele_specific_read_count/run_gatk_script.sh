#!/bin/bash
module load samtools/1.9
#Generate a file with the path to all neutral bam files for North America
ls input/wasp_bam/*bam > path_to_bam.txt

#Generate a file with only the smaple ids
cat path_to_bam.txt | cut -d '/' -f 3 | cut -d '.' -f 1 > samples.ids

gunzip -c input/Si_gnG.fna.gz > tmp/Si_gnG.fna
#GATK needs both an index and a dict file to run:
samtools faidx tmp/Si_gnG.fna
java -jar picard CreateSequenceDictionary REFERENCE=tmp/Si_gnG.fna OUTPUT=tmp/Si_gnG.fna.dict

#GATK needs a sorted vcf to run
java -jar picard SortVcf \
      I=input/subset.vcf \
      O=tmp/sorted_subset.vcf

qsub script_gatk_neutral_ref_north_america.array.sh
qsub script_gatk_neutral_ref_south_america.array.sh

