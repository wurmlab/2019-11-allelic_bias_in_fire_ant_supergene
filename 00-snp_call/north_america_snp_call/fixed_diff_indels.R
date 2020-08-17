#!/usr/bin/env Rscript

library(VariantAnnotation)

vcf <- readVcf("results/indels/unfiltered_indels.vcf.gz", "Si_gnGA")

total_pop <- samples(header(vcf))

bb <- total_pop[grepl("bigB", total_pop)]
lb <- total_pop[grepl("littleb", total_pop)]

stopifnot(total_pop %in% c(bb, lb))

gt <- geno(vcf)$GT

bb_gt <- gt[,bb]
lb_gt <- gt[,lb]

fixed_bb <- apply(bb_gt, 1, function(x) all(x == x[1]))
fixed_lb <- apply(lb_gt, 1, function(x) all(x == x[1]))

bb1_diff_lb1 <- sapply(1:nrow(bb_gt), function(i) bb_gt[i,1] != lb_gt[i,1])


fixed_diff_index <- fixed_bb & fixed_lb & bb1_diff_lb1

fixed_diff_vcf <- vcf[fixed_diff_index, ]

writeVcf(fixed_diff_vcf, filename = "results/indels/unfiltered_fixed_different_indels.vcf")

# a <- load("filtered_cds.rda")
# length(cds_ranges[cds_ranges %over% fixed_diff_vcf])
