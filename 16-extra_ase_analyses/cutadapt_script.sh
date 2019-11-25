#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=2:00:00
#$ -t 1-8

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/samples.ids)
#Script to remove universal Illumina adapters (based on Priyam's pipeline:
#(https://raw.githubusercontent.com/wurmlab/reads-qc/master/Rakefile?token=AJfX3aTIB4hskD8DguRqqLVdgaev78veks5a2GQOwA%3D%3D)

#Run cutadapt (v1.18)------------------------------------------------------
module load python
source ~/cutadapt/bin/activate
module unload python
#Find at least 4 bases of that sequence common to allthe Universal Ilumina Adapters and Tru Seq Adapters in reverse and forward,
#remove empty reads
cutadapt -O 4 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
-a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGTTAACTATCTCGTATGCCGTCTTCTGCTTGA \
-A GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGTTAACTATCTCGTATGCCGTCTTCTGCTTGA \
-m 1 \
-o tmp/trimmed_reads_fontana_etal/trimmed_${INPUT_FILE}_1.fastq.gz \
-p tmp/trimmed_reads_fontana_etal/trimmed_${INPUT_FILE}_2.fastq.gz \
input/reads_fontana_etal/${INPUT_FILE}_1.fastq.gz input/reads_fontana_etal/${INPUT_FILE}_2.fastq.gz

deactivate
