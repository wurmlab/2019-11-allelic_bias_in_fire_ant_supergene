#Package VariantAnnotation version 1.24.5
library(VariantAnnotation)

#Load vcf
vcf <- readVcf(file = "tmp/fixed_diffs/supergene.vcf", genome = "Si_gnG")
#Extract the name of the samples used for generating the vcf and separate them in Sb vs SB
total_pop <- samples(header(vcf))
bb <- total_pop[grepl("bigB", total_pop)]
lb <- total_pop[grepl("littleb", total_pop)]

#Get the genotypes for all individuals within each group
gt <- geno(vcf)$GT
bb_gt <- gt[,bb]
lb_gt <- gt[,lb]

#Extract only those SNPs which are fixed within either SB or Sb
fixed_bb <- apply(bb_gt, 1, function(x) all(x == x[1]))
fixed_lb <- apply(lb_gt, 1, function(x) all(x == x[1]))

#Extract the fixed SNPs which are different between SB and Sb
bb1_diff_lb1 <- sapply(1:nrow(bb_gt), function(i) bb_gt[i,1] != lb_gt[i,1])

#Generate the new vcf file
fixed_diff_index <- fixed_bb & fixed_lb & bb1_diff_lb1
fixed_diff_vcf <- vcf[fixed_diff_index, ]
writeVcf(fixed_diff_vcf, filename = "tmp/supergene_fixed_diffs_north_america.vcf")
