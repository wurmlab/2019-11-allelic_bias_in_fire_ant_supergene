#!/usr/bin/env Rscript

library(VariantAnnotation)

vcf <- readVcf("tmp/fixed_diffs/supergene.vcf", "Si_gnGA")

total_pop <- samples(header(vcf))
total_pop <- total_pop[grepl("AR", total_pop)]

bb <- total_pop[grepl("bigB", total_pop)]
lb <- total_pop[grepl("littleb", total_pop)]


gt <- geno(vcf)$GT

bb_gt <- gt[,bb]
lb_gt <- gt[,lb]

fixed_bb <- apply(bb_gt, 1, function(x) all(x == x[1]))
fixed_lb <- apply(lb_gt, 1, function(x) all(x == x[1]))

bb1_diff_lb1 <- sapply(1:nrow(bb_gt), function(i) bb_gt[i,1] != lb_gt[i,1])


fixed_diff_index <- fixed_bb & fixed_lb & bb1_diff_lb1

fixed_diff_vcf <- vcf[fixed_diff_index, ]

writeVcf(fixed_diff_vcf, filename = "results/fixed_differences/2018-09-01-AR_fixed_differences.vcf")

# a <- load("filtered_cds.rda")
# length(cds_ranges[cds_ranges %over% fixed_diff_vcf])
