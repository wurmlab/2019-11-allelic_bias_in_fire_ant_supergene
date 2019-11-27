#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -t 1-8
module load star/2.6.1a
module load samtools/1.9
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/samples.ids)

#Run all samples in parallel. Multiple mapping reads are not allowed. Only one random alignment of the multiply mapped reads
#will be considered

STAR --runThreadN 4 --genomeDir tmp/star_reference_fontana_etal \
     --readFilesIn tmp/trimmed_reads_fontana_etal/trimmed_${INPUT_FILE}_1.fastq.gz tmp/trimmed_reads_fontana_etal/trimmed_${INPUT_FILE}_2.fastq.gz \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files/${INPUT_FILE} --outSAMattributes All \
     --sjdbFileChrStartEnd tmp/bam_files/all_samples_fontana_etalSJ.out.tab --outMultimapperOrder Random --outSAMmultNmax 1 \
     --readFilesCommand zcat

samtools index tmp/bam_files/${INPUT_FILE}Aligned.sortedByCoord.out.bam
