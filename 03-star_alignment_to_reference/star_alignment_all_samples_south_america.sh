#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_rt=6:0:0
#$ -l h_vmem=4G
module load star/2.6.1a

mkdir tmp/star_reference_south_america
gunzip -c input/reference.fna.gz > tmp/reference.fna

#Generate an index for STAR using refseq's gnG assembly and annotations, specific for South American populations
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir tmp/star_reference_south_america \
     --genomeFastaFiles tmp/reference.fna \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFfile input/GCF_000188075.1_Si_gnG_genomic.gff --sjdbOverhang 149 \
     --limitGenomeGenerateRAM 60000000000

#Align the South American reads
mkdir tmp/bam_files
#Generate two variables with all the fastq files from South America
R1_READS=$(ls input/ar_trimmed_reads/*R1*fastq.gz | cat | tr '\n' ',' | sed 's/,$//')
R2_READS=$(ls input/ar_trimmed_reads/*R2*fastq.gz | cat | tr '\n' ',' | sed 's/,$//')

#Run STAR in automatic two pass mode (STAR will run twice, the second time inserting splice junctions detected
# during the first run).
STAR --runThreadN 16 --genomeDir tmp/star_reference_south_america \
     --readFilesIn ${R1_READS} ${R2_READS} \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files/all_samples_south_america\
     --outSAMattributes All --twopassMode Basic --readFilesCommand zcat
