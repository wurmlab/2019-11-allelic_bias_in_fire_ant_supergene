#!/bin/sh
module load parallel
module load genometools/1.5.9
module load python/3.6.3
module load ruby/2.4.1

#Unzip the ensembl reference
gunzip -c input/Si_gnG_ensembl.fa.gz > tmp/Si_gnG_ensembl.fa
#Use flo to generate a chain file between gnG and Ensembl, then use it to lift over the gnG vcf file.
#Flo was created by Anurag Priyam and the explanation on how to use it is here:
#https://github.com/wurmlab/flo
#Copy yaml file
cp ~/software/flo/opts_example.yaml ./flo_opts.yaml
#Modify .yaml file using a text editor, in this case, vim. Comment out akk the 'add to path' commands
#source_fa: 'input/Si_gnG.fa'
#target_fa: 'input/Si_gnG_ensembl.fa'
#processes: '24'
#Comment out the two last lines of the yaml file, as we do not need gff outputs

#Run flo
rake -f ~/software/flo/Rakefile

#Move the run directory into tmp
mv run tmp/

#Run CrossMap to liftover the vcf from gnGA to gnG
source ~/crossmaps/bin/activate
CrossMap.py vcf tmp/run/liftover.chn input/subset_south_america.vcf tmp/Si_gnG_ensembl.fa tmp/ensembl_sa.vcf
deactivate
