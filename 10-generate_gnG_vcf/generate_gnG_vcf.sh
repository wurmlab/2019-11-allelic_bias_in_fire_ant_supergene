module load parallel
module load genometools/1.5.9
module load python/3.6.3

#Use flo to generate a chain file between gnGA and gnG, then use it to lift over the gnGA vcf file.
#Flo was created by Anurag Priyam and the explanation on how to use it is here:
#https://github.com/wurmlab/flo
#Copy yaml file
cp ~/software/flo/opts_example.yaml ./flo_opts.yaml
#Modify .yaml file using a text editor, in this case, vim. Comment out akk the 'add to path' commands
#source_fa: 'input/Si_gnGA.fa'
#target_fa: 'input/Si_gnG.fa'
#processes: '24'
#Comment out the two last lines of the yaml file, as we do not need gff outputs

#Run flo
rake -f ~/software/flo/Rakefile

#Move the run directory into tmp
mv run tmp/

#Unzip the vcf file
gunzip -c input/gnGA_sa.vcf.gz > tmp/gnGA_sa.vcf
#Run CrossMap to liftover the vcf from gnGA to gnG
source ~/crossmaps/bin/activate
CrossMap.py vcf tmp/run/liftover.chn tmp/gnGA_sa.vcf input/Si_gnG.fa tmp/gng_sa.vcf
deactivate
