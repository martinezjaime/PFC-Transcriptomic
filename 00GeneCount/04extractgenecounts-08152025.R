##################################################################
# Extract Gene counts and annotation data
# Author: Jose Jaime Martinez-Magana
# Day: 07272025
##################################################################
# Load libraries
library(tidyr)
library(data.table)
##################################################################
# Set working directory
wkdir <- 'C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01transcriptome/00data/'
setwd(wkdir)
##################################################################
# Load datasets
uthbP <- "matrix_counts_featurecounts-07252025-03.rds"
vabbP <- "matrix_counts_featurecounts-07252025-02.rds"

uthb <- readRDS(uthbP)
vabb <- readRDS(vabbP)
##################################################################
# Fix sample names
# uthb
colnames(uthb$counts) <- sub("_", "", sub("_Aligned.*", "", colnames(uthb$counts)))

# vabb
new_colnames <- sub("^((Sample_[^_]+)).*", "\\1", colnames(vabb$counts))
new_colnames <- sub("Aligned.sortedByCoord.out.bam", "", new_colnames)
new_colnames <- sub("_r", "", new_colnames)
new_colnames <- sub("_", "", new_colnames)
colnames(vabb$counts) <- new_colnames
##################################################################
# Keep them separate
expr_uthb  <- uthb$counts
annot_uthb <- uthb$annotation

expr_vabb  <- vabb$counts
annot_vabb <- vabb$annotation
##################################################################
# If you want to combine:
# Make sure genes (rows) are the same and in same order
if (identical(rownames(expr_uthb), rownames(expr_vabb))) {
  expr_combined  <- cbind(expr_uthb, expr_vabb)
  annot_combined <- rbind(annot_uthb, annot_vabb)
}
##################################################################
# Save counts and annotation tables to CSV
##################################################################
# Save uthb counts as gzipped CSV
write.csv(expr_uthb, file = gzfile("00uthb_genecounts-08152025.csv.gz"), row.names = TRUE)

# Save uthb annotation as gzipped CSV
write.csv(annot_uthb, file = gzfile("00uthb_annotation-08152025.csv.gz"), row.names = TRUE)

# vabb counts
write.csv(expr_vabb, file = gzfile("00vabb_genecounts-08152025.csv.gz"), row.names = TRUE)

# vabb annotation
write.csv(annot_vabb, file = gzfile("00vabb_annotation-08152025.csv.gz"), row.names = TRUE)

# Optional combined datasets
if (exists("expr_combined") && exists("annot_combined")) {
  write.csv(expr_combined, file = gzfile("00combined_genecounts-08152025.csv.gz"), row.names = TRUE)
  write.csv(annot_combined, file = gzfile("00combined_annotation-08152025.csv.gz"), row.names = TRUE)
}
