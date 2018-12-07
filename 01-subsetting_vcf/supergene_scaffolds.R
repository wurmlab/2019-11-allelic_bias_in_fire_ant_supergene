#Load the table with the positions of all scaffolds in the gnG assembly of S.invicta.
gng_regions <- read.table("input/gng_regions.tab", header = TRUE)

#Subset the table by selecting only scaffolds in the supergene region
supergene_gng_regions <- subset(x = gng_regions, subset = region == "supergene")

#Extract only the name of the scaffolds
supergene_scaffolds <- as.character(supergene_gng_regions$scaffold)

#Output
write.table(supergene_scaffolds, "tmp/supergene_scaffolds.ids", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
#18 scaffolds in total
