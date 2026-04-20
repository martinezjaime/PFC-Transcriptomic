################################################################################
# Deep Neural Age ensemble (DNAe) - training
# 2020, NH
# Training of pathway-based neural networks using the hallmark gene sets for ensemble construction
################################################################################
# Load libraries
# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(patchwork)
  library(plot3D)
  library(pals)
  library(ggExtra)
})
################################################################################
## This section is UNIQUE for your dataset
# setting working directory
# Directory for output
# Set working directory
wkdir='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/04geocohorts-08032025/'
setwd(wkdir)
# loading expression database
data = readRDS("00gse102556-08052025.rds")
expressionData = data$expressionData
phenoData = data$phenoData

# Subsetting to complete cases
phenoData$SampleID = phenoData$sample
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
phenoData$intro_age <- as.numeric(phenoData$age)
phenoDataSub = phenoData[,c('SampleID', 'age')]
colnames(phenoDataSub) = c('SampleID', 'Age')
write.csv(file='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/00Age_gse102556-08052025.csv',
          phenoDataSub,
          quote=FALSE,
          row.names = FALSE)
# Save as gzipped CSV
write.csv(
  expressionData,
  file = gzfile("C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01GeneMatrix_gse102556-08052025.csv.gz"),
  quote = FALSE,
  row.names = TRUE
)
##################################################################