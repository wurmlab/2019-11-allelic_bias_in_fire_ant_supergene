#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=4:00:00
#$ -l h_vmem=8G
#$ -t 1-6
module load star/2.6.1a
module load samtools/1.9
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/samples.txt)
#Run all samples in parallel. Multiple mapping reads are not allowed. Only one random alignment of the multiply mapped reads
#will be considered
STAR --runThreadN 4 --genomeDir tmp/star_reference \
     --readFilesIn input/wurm_2011_reads/${INPUT_FILE}.fastq.gz \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files/${INPUT_FILE} --outSAMattributes All \
     --sjdbFileChrStartEnd input/all_samples_north_americaSJ.out.tab --outMultimapperOrder Random --outSAMmultNmax 1 \
     --readFilesCommand zcat

samtools index tmp/bam_files/${INPUT_FILE}Aligned.sortedByCoord.out.bam
