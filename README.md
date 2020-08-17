# Scripts for "Genomic architecture and evolutionary conflict drive allele-specific expression in the social supergene of the red fire ant"

Scripts for the generation and analysis of the data used for "Genomic architecture and evolutionary conflict drive allele-specific expression in the social supergene of the red fire ant". The scripts listed here go through all the steps from raw RNAseq data, to read count and analysis of differences in gene expression between social forms and between allele-specific expression in the red fire ant *Solenopsis invicta*.

The raw data used for this analysis is available in NCBI, specifically:
* RNAseq and DNAseq data from South American populations: [PRJNA542606](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA542606).
* RNAseq data from North American populations: [PRJDB4088](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB4088) and [PRJNA49629](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA49629).
* RNAseq data from Taiwanese populations: [PRJNA542500](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA542500).

All the subdirectories are numbered according to the order in which they were run for the analysis:
1. **subsetting_vcf:** Extracting variants from the VCF that fall within the supergene region and with fixed differences between SB and Sb. The same script was run in the VCFs based in both South and North American populations.
2. **littleb_ified_reference:** Transform the gnG version of the red fire ant genome reference by replacing the SB variants from the reference by Sb variants from the VCFs generated in the previous step. The same script was run in the VCFs based in both South and North American populations.
3. **star_alignment_to_reference:** Align the raw RNAseq reads to the reference of the red fire ant. The same script was run for data from both South and North America populations and for the normal and Sb-transformed references.
4. **wasp_reference_bias_free:** Generate a reference-bias free alignment BAM file using the alignments generated in the previous step. The same script was run for data from both South and North America populations.
5. **allele_specific_read_count:** Generate allele-specific read counts for the supergene region of the red fire ant based on the VCFs and reference-bias free BAM files generated in previous steps. The same script was run for data from both South and North America populations.
6. **social_form_comparison:** Generate read counts from raw RNAseq data for North American populations. Test for significant differences in gene expression between social forms using DESeq2.
7. **differential_expression_analyses_americas:** Test for allele-specific expression differences between variants in each population independently using DESeq2 for North and South American populations.
8. **taiwan_ase_analysis:** Perform allele specific expression analyses and test for expression differences between variants in Taiwanese populations.
9. **combined_ase_analysis:** Test for allele-specific expression difference between variants using Tawianese, North and South American populations using a linear mixed effects model. Generate plots to visualise allele-specific expression differences between variants. Only the script using all populations was used for the manuscript.
10. **social_form_ase_analysis:** Analyse together the allele-specific expression data and social form gene expression differences from North America
11. **generate_gnG_vcf:** Transform the original whole-genome VCF for North American populations from a gnGA assembly to a gnG assembly.
12. **whole_genome_ase:** Generate allele-specific read counts for the whole genome of *Solenopsis invicta*.
13. **whole_genome_ase_analysis:** Use DESeq2 to test for allele-specific expression differences in the whole genome of *Solenopsis invicta*.
14. **snp_effects:** Run SNPEff on the variants detected in the supergene to asses their potential impact in protein coding sequences.
15. **snp_effects_analysis:** Analyse the results generated in the previous step to test for Sb lower expression with increased nonsynonymous mutations. Check for effects in specific gene lists.
16. **indel_subsetting:** Select indels in the supergene region with fixed differences between SB and Sb in both populations of *Solenopsis invicta*. Not used in the manuscript.
17. **mt_diversity:** Analyze mitocondrial diversity between different South American populations of *Solenopsis invicta*. Not used in the manuscript.
