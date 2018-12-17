#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_rt=2:0:0
#$ -l h_vmem=4G
module load star/2.6.1a

mkdir tmp/star_reference_north_america
gunzip -c input/reference.fna.gz > tmp/reference.fna

#Generate index using refseq's gnG assembly and annotations
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir tmp/star_reference_north_america \
     --genomeFastaFiles tmp/reference.fna \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFfile input/GCF_000188075.1_Si_gnG_genomic.gff --sjdbOverhang 74 \
     --limitGenomeGenerateRAM 60000000000

#Align the North American reads
mkdir tmp/bam_files
#Generate two variables with all the fastq files from North America. Sample SRR619956 and SRR619959 seem to be the same
#by looking at the SRA metadata on NCBI. Remove one of them
READS=$(ls input/wurm_2011_reads/*fastq.gz | cat | sed '/SRR619959/d' | tr '\n' ',' | sed 's/,$//')

#Run STAR in automatic two pass mode (STAR will run twice, the second time inserting splice junctions detected
# during the first run).
STAR --runThreadN 16 --genomeDir tmp/star_reference_north_america \
     --readFilesIn ${READS} \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files/all_samples_north_america\
     --outSAMattributes All --twopassMode Basic --readFilesCommand zcat
