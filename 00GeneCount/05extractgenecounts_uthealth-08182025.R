##################################################################
# Extract Gene counts and annotation data
# Author: Jose Jaime Martinez-Magana
# Day: 07272025
##################################################################
# Load libraries
library(tidyr)
library(data.table)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
##################################################################
# Set working directory
wkdir <- 'C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01transcriptome/00data/'
setwd(wkdir)
##################################################################
# Load datasets
expressionPath <- "matrix_counts_featurecounts-07252025-03.rds"
vabb = readRDS(expressionPath)

# Apply to your actual column names vector
new_colnames <- sub("Aligned.sortedByCoord.out.bam", "", colnames(vabb$counts))
new_colnames <- sub("_", "", new_colnames)
new_colnames <- sub("_", "", new_colnames)
# Optionally assign them back
colnames(vabb$counts) <- new_colnames
expressionData <- vabb$counts
expressionData = as.data.frame(expressionData)

# loading phenotype data
phenoPath = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/00phenotype/01uthealth/00MetaData_UTHealth_RNAseq-07292025.csv"
phenoData <- read.csv(phenoPath, header = T)

# Subsetting to complete cases
phenoData$SampleID = phenoData$ID
# Step 1: Get intersecting sample IDs
shared_samples <- intersect(colnames(expressionData), phenoData$SampleID)
# Step 2: Subset expressionData to shared samples
expressionData <- expressionData[, shared_samples]
# Step 3: Subset and reorder phenoData to match expressionData columns
phenoData <- phenoData[phenoData$SampleID %in% shared_samples, ]
phenoData <- phenoData[match(shared_samples, phenoData$SampleID), ]
# Step 4: Double-check alignment
stopifnot(all(phenoData$SampleID == colnames(expressionData)))
# Step 5: Reassign `zz_nr` and `intro_age` as before
phenoData$zz_nr <- phenoData$SampleID
phenoData$intro_age <- phenoData$Age
################################################################################
# mapping Entrez IDs to Ensembl IDs
mapEntrezIDtoEnsembl <- function(entrez_id) {
  ensembl <- mapIds(org.Hs.eg.db,
                    keys = entrez_id,
                    column = "ENSEMBL",
                    keytype = "ENTREZID",
                    multiVals = "first")  # choose first if multiple mappings
  return(ensembl)
}

# sum observations by column
sumByColumn <- function(x, column) {
  x <- as.data.table(x)
  x2 <- x[, lapply(.SD, sum), by = column]
  return(as.data.frame(x2))
}

# map and summarize gene expression data from EntrezID to Ensembl
summarizeToEnsembl <- function(data) {
  data$ensembl <- mapEntrezIDtoEnsembl(data$gene_id)
  data$gene_id <- NULL
  data <- sumByColumn(na.omit(data), "ensembl")
  rownames(data) <- data$ensembl
  data$ensembl <- NULL
  return(data)
}
##################################################################
# Adjust script to remove the symbol row
expressionData$gene_id = rownames(expressionData)
expressionData <- summarizeToEnsembl(expressionData)
expressionData <- as.matrix(expressionData)
##################################################################
# Save counts and annotation tables to CSV
##################################################################
# Subsetting to complete cases
phenoData$SampleID = phenoData$ID
# Step 1: Get intersecting sample IDs
shared_samples <- intersect(colnames(expressionData), phenoData$SampleID)
# Step 2: Subset expressionData to shared samples
expressionData <- expressionData[, shared_samples]
# Step 3: Subset and reorder phenoData to match expressionData columns
phenoData <- phenoData[phenoData$SampleID %in% shared_samples, ]
phenoData <- phenoData[match(shared_samples, phenoData$SampleID), ]
# Step 4: Double-check alignment
stopifnot(all(phenoData$SampleID == colnames(expressionData)))
# Step 5: Reassign `zz_nr` and `intro_age` as before
phenoData$zz_nr <- phenoData$SampleID
phenoData$intro_age <- as.numeric(phenoData$Age)
phenoDataSub = phenoData[,c('SampleID', 'Age')]
colnames(phenoDataSub) = c('SampleID', 'Age')
write.csv(file='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/00Age_uthealth-08182025.csv',
          phenoDataSub,
          quote=FALSE,
          row.names = FALSE)
# Save as gzipped CSV
write.csv(
  expressionData,
  file = gzfile("C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01GeneMatrix_uthealth-08182025.csv.gz"),
  quote = FALSE,
  row.names = TRUE
)
##################################################################