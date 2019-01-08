#!/bin/bash
#Load python virtual environment to run WASP
module load python/3.6.3
source ~/wasp/bin/activate
module load samtools/1.9
module load parallel/20170422
#Generate samples.txt with the name of the bam files. Make sure that the names are the same in both directories
ls input/bigB_bam/*\.bam | cut -d '/' -f 3 | sed 's/Aligned.*//g' > tmp/samples_bb.txt
ls input/littleb_bam/*\.bam | cut -d '/' -f 3 | sed 's/Aligned.*//g' > tmp/samples_lb.txt
cmp tmp/samples_bb.txt tmp/samples_lb.txt
status=$?
if [[ $status != 0 ]]; then
    echo "ERROR: Files are different"
    exit 1
fi
cat tmp/samples_lb.txt tmp/samples_bb.txt| sort -u > tmp/samples.txt


#Step 6 of the WASP pipeline for mapping, merge bam files aligned to
parallel 'samtools merge - input/bigB_bam/{}Aligned.sortedByCoord.out.bam\
          input/littleb_bam/{}Aligned.sortedByCoord.out.bam | \
          samtools sort -o tmp/{}_merged_sorted.bam -' :::: tmp/samples.txt

parallel 'samtools index tmp/{}_merged_sorted.bam' :::: tmp/samples.txt

#Step 7 remove duplicated reads using the rmdup script from WASP
parallel 'python rmdup tmp/{}_merged_sorted.bam results/{}ref_bias_free.bam' :::: tmp/samples.txt

#Index the newly generated bam files
ls results/*bam | parallel 'samtools index'

deactivate
#Add reading groups (RGs) to the wasp generated bam files. This will generated a new set of bam
# files. This is so that GATK can work for extracting the allele specific counts in the enxt step.

parallel 'java -jar picard AddOrReplaceReadGroups I=results/{}ref_bias_free.bam \
          O=results/with_rg/{}.ref_bias_free.rg.bam \
          RGID=1 RGLB=lib{} RGPL=illumina RGPU=unit1 RGSM={}' :::: tmp/samples.txt

#Index the newly generated bam files
ls results/with_rg/*bam | parallel 'samtools index'
