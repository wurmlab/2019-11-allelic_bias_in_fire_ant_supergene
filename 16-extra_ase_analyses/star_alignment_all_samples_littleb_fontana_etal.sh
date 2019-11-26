#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_rt=12:0:0
#$ -l h_vmem=4G
module load star/2.6.1a

mkdir tmp/star_reference_fontana_etal_lb
#Use little-b reference based on North American samples (Taiwan pops come from North American ones)
gunzip -c input/reference_lb.fna.gz > tmp/reference_lb.fna

#Generate an index for STAR using refseq's gnG assembly and annotations
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir tmp/star_reference_fontana_etal_lb \
     --genomeFastaFiles tmp/reference_lb.fna \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFfile input/GCF_000188075.1_Si_gnG_genomic.gff --sjdbOverhang 125 \
     --limitGenomeGenerateRAM 60000000000

#Align the  reads
mkdir tmp/bam_files_lb
#Generate two variables with all the fastq files
R1_READS=$(ls tmp/trimmed_reads_fontana_etal/*_1.fastq.gz | cat | tr '\n' ',' | sed 's/,$//')
R2_READS=$(ls tmp/trimmed_reads_fontana_etal/*_2.fastq.gz | cat | tr '\n' ',' | sed 's/,$//')

#Run STAR in automatic two pass mode (STAR will run twice, the second time inserting splice junctions detected
# during the first run).
STAR --runThreadN 16 --genomeDir tmp/star_reference_fontana_etal_lb \
     --readFilesIn ${R1_READS} ${R2_READS} \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files_lb/all_samples_fontana_etal_lb\
     --outSAMattributes All --twopassMode Basic --readFilesCommand zcat
