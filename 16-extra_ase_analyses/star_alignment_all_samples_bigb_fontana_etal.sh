#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_rt=12:0:0
#$ -l h_vmem=4G
module load star/2.6.1a

mkdir tmp/star_reference_fontana_etal
gunzip -c input/reference.fna.gz > tmp/reference.fna

#Generate an index for STAR using refseq's gnG assembly and annotations
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir tmp/star_reference_fontana_etal \
     --genomeFastaFiles tmp/reference.fna \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFfile input/GCF_000188075.1_Si_gnG_genomic.gff --sjdbOverhang 125 \
     --limitGenomeGenerateRAM 60000000000

#Align the  reads
mkdir tmp/bam_files
#Generate two variables with all the fastq files
R1_READS=$(ls tmp/trimmed_reads_fontana_etal/*_1.fastq.gz | cat | tr '\n' ',' | sed 's/,$//')
R2_READS=$(ls tmp/trimmed_reads_fontana_etal/*_2.fastq.gz | cat | tr '\n' ',' | sed 's/,$//')

#Run STAR in automatic two pass mode (STAR will run twice, the second time inserting splice junctions detected
# during the first run).
STAR --runThreadN 16 --genomeDir tmp/star_reference_fontana_etal \
     --readFilesIn ${R1_READS} ${R2_READS} \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files/all_samples_fontana_etal\
     --outSAMattributes All --twopassMode Basic --readFilesCommand zcat
