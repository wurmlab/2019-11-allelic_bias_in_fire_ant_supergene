common_snps <- read.table("common_snps.txt")

all_info <- as.character(common_snps$V8)
all_info <- gsub(".+;ANN=[A-Z]*|", replacement = "", all_info)

all_info <- strsplit(x = all_info, split = "|", fixed = TRUE)

snp_position <- sapply(all_info, function(x) x[2])
snp_effect <- sapply(all_info, function(x) x[3])
snp_gene <- sapply(all_info, function(x) x[4])

snp_info <- as.character(common_snps$V3)

scaffold <- gsub(pattern = ":.+", replacement = "", x = snp_info)
position <- gsub(pattern = "(.+:)([0-9]+)(_.+)", replacement = "\\2", x = snp_info)
reference <- gsub(pattern = "(.+_)([A-Z])(/[A-Z])", replacement = "\\2", x = snp_info)
alternative <- gsub(pattern = "(.+_)([A-Z]/)([A-Z])", replacement = "\\3", x = snp_info)

snp_data <- data.frame("Scaffold" = scaffold, "Position_in_scaffold" = position, "Reference_allele" = reference,
                       "Alternative_allele" = alternative, "Position_in_gene" = snp_position, "Gene-s-_affected" = snp_gene,
                       "SNP_effect" = snp_effect)


colnames(snp_data) <- gsub(pattern = "_",
                           replacement = " ",
                           x = colnames(snp_data))
snp_data$`Position in gene` <-  gsub(pattern = "_",
                                     replacement = " ",
                                     x = snp_data$`Position in gene`)


write.csv(file = "snp_info.csv", x = snp_data, quote = FALSE, row.names = FALSE)



