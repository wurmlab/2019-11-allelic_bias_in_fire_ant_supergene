#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-g", "--gt"), type="character", default=NULL,
              help="GT matrix", metavar="character"),
  make_option(c("-a", "--ao"), type="character", default=NULL,
              help="AO matrix", metavar="character"),
  make_option(c("-r", "--ro"), type="character", default=NULL,
              help="RO matrix", metavar="character"),
  make_option(c("-d", "--mean_dp"), type="numeric", default=30,
              help="Maximum median DP value [default= %default]", metavar="character"),
  make_option(c("-c", "--min_cov"), type="numeric", default=1,
              help="Minimum individual coverage value [default= %default]", metavar="character"),
  make_option(c("-p", "--called_prop"), type="numeric", default=0.6,
              help="Minimum proportion of reads supporting the called allele [default= %default]", metavar="character"),
  make_option(c("-s", "--outsummary"), type="character", default="out.summary",
              help="Name of table with summary [default= %default]", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="out.summary",
              help="Name of file to write filter column [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Inputs
# Matrices with several parameters taken from VCF

# Matrices with the read count supporting:
#      - ro: the reference allele
#      - ao: the alternative allele
geno_ao <- read.table(opt$ao)
geno_ao <- as.matrix(geno_ao)

geno_ro <- read.table(opt$ro)
geno_ro <- as.matrix(geno_ro)

# Genotype matrix
geno_gt <- read.table(opt$gt)
geno_gt <- as.matrix(geno_gt)

# Filter thresholds
mean_dp_threshold           <- opt$mean_dp
lower_cov_threshold         <- opt$min_cov
non_called_support_threhold <- opt$called_prop

# Outputs
filter_summary_out <- opt$outsummary
filter_out         <- opt$outfile

# ---------------------------------------------------------------------------- #
#  "." denote absence of reads -> turn to "0"

geno_ao[geno_ao == "."] <- "0"
geno_ro[geno_ro == "."] <- "0"

geno_ro <- apply(geno_ro, 2, as.numeric)
geno_ao <- apply(geno_ao, 2, as.numeric)

# ---------------------------------------------------------------------------- #
# Checks

# First check: all sites need to be polymorphic
not_polymorphic <- apply(geno_gt, 1, function(x) all(x == x[1]))

if (any(not_polymorphic)) {
  stop("All sites need to be polymorphic")
}

if (!((nrow(geno_ao) == nrow(geno_ro)) & (nrow(geno_ao) == nrow(geno_gt)))) {
  stop("All supplied matrices need to have the same number of sites")
}

if (!((ncol(geno_ao) == ncol(geno_ro)) & (ncol(geno_ao) == ncol(geno_gt)))) {
  stop("All supplied matrices need to have the same number samples")
}

# ---------------------------------------------------------------------------- #
# Make a matrix of the called allele (co):
#      The called allele is either the reference allele or the alternative
#      allele, depending on which was called.

# co == called allele
# Make empty matrix
co <- matrix(-9, nrow = nrow(geno_gt), ncol = ncol(geno_gt))
# Populate with the reference alleles
co[geno_gt == "0"] <- geno_ro[geno_gt == "0"]
# Populate with the alternative alleles
co[geno_gt == "1"] <- geno_ao[geno_gt == "1"]

# Populate with the fields where no allele was Called
# When no genotype is called, make sure there are no reads for either allele
stopifnot(geno_ro[geno_gt == "."] == 0)
stopifnot(geno_ao[geno_gt == "."] == 0)
co[geno_gt == "."] <- 0

# Make sure all is populated
stopifnot(co >= 0)

# Total read depth per individual
dp <- geno_ro + geno_ao

# cf == called fraction = number of reads supporting the called allele / total
cf <- co / dp

# ---------------------------------------------------------------------------- #
# Checks
stopifnot(!is.na(dp))

# ---------------------------------------------------------------------------- #
# Filter: filter by maximum read depth (or maybe by median read depth)
mean_dp <- apply(dp, 1, mean)

no_individual_has_too_much_cov <- mean_dp <= mean_dp_threshold

# ---------------------------------------------------------------------------- #
# Filter: sites where at least one individual has no read coverage
all_individuals_have_cov <- dp >= lower_cov_threshold
all_individuals_have_cov <- apply(all_individuals_have_cov, 1, all)

# ---------------------------------------------------------------------------- #
# Filter: where one or more individuals have too many reads supporting the
#        non-called allele
noIndividualsHaveAmbigousSupport <- function(x, non_called_support_threhold) {
  if (all(is.na(x))) {
    return(TRUE)
  } else {
    x <- x[!is.na(x)]
    return(all(x >= non_called_support_threhold))
  }
}

all_individuals_support_called <- apply(cf,
                                        1,
                                        function(x)
              noIndividualsHaveAmbigousSupport(x, non_called_support_threhold))

# ---------------------------------------------------------------------------- #
# Filter: site QUALITY (if necessary)
# quality_threshold <- 25

# ---------------------------------------------------------------------------- #
# Paste together the different filters to make a column of filters

high_dp                                                <- rep("NA", nrow(geno_gt))
high_dp[!no_individual_has_too_much_cov]               <- "HIGH_DP"

no_cov                                                 <- rep("NA", nrow(geno_gt))
no_cov[!all_individuals_have_cov]                      <- "NO_COV"

ambiguous_read_support                                  <- rep("NA", nrow(geno_gt))
ambiguous_read_support[!all_individuals_support_called] <- "AMBIGUOUS_SUPPORT"

# Paste together
total_filter <- paste(high_dp, no_cov, ambiguous_read_support, sep=",")
total_filter[no_individual_has_too_much_cov &
             all_individuals_have_cov &
             all_individuals_support_called] <- "PASS"

total_filter <- gsub("(\\,*NA\\,*)+", "", total_filter)

stopifnot(!grepl(",,", total_filter))

# ---------------------------------------------------------------------------- #
# Write out
total_filter_df <- data.frame(total_filter)
write.table(total_filter_df, file = filter_out, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# ---------------------------------------------------------------------------- #
# Summary table
filter_summary <- data.frame(
  filter = c("pre-filter", "HIGH_DP", "NO_COV", "AMBIGUOUS_SUPPORT", "multiple_filters", "total_filtered_out", "total_remaining"),
  number = c(nrow(dp),
             sum(!no_individual_has_too_much_cov),
             sum(!all_individuals_have_cov),
             sum(!all_individuals_support_called),
             sum(grepl(",", total_filter)),
             sum(total_filter != "PASS"),
             sum(total_filter == "PASS"))
)

filter_summary$percentage <- round(100 * filter_summary$number / filter_summary$number[1], 1)

write.table(filter_summary, file = filter_summary_out, quote = FALSE,
            col.names = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------- #
