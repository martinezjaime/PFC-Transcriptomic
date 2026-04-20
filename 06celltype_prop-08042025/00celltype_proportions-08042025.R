###############################################################################
# Cell-type proportions and surrogate variables in RNASeq data
# Day 08042025
# Author: Jose Jaime Martinez-Magana
###############################################################################
# Set working directory
wkdir='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/05celltype-08042025/'
setwd(wkdir)
##################################################################
# load libraries
library(RUVSeq)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(BRETIGEA)
library(knitr) #only for visualization
##################################################################
# Load datasets
uthbP = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01transcriptome/00data/matrix_counts_featurecounts-07252025-03.rds"
vabbP = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01transcriptome/00data/matrix_counts_featurecounts-07252025-02.rds"

uthb = readRDS(uthbP)
vabb = readRDS(vabbP)

# Set the path for the output object
vabb_cells = "00vabb_cellprop_svas-08042025.csv"
uthh_cells = "01uthealth_cellprop_svas-08042025.csv"
##################################################################
# Fixing names of samples
### VABB
# Apply to your actual column names vector
new_colnames <- sub("^((Sample_[^_]+)).*", "\\1", colnames(vabb$counts))
new_colnames <- sub("Aligned.sortedByCoord.out.bam", "", new_colnames)
new_colnames <- sub("_r", "", new_colnames)
new_colnames <- sub("_", "", new_colnames)
new_colnames <- sub("^Sample(\\d+)[A-Za-z]+(\\d{1,2})$", "Sample\\1_\\2", new_colnames)
# Optionally assign them back
colnames(vabb$counts) <- new_colnames
expressionData1 <- vabb$counts
expressionData1 = as.data.frame(expressionData1)


### UTHealth
# Apply to your actual column names vector
new_colnames <- sub("Aligned.sortedByCoord.out.bam", "", colnames(uthb$counts))
new_colnames <- sub("_", "", new_colnames)
new_colnames <- sub("_", "", new_colnames)
# Optionally assign them back
colnames(uthb$counts) <- new_colnames
expressionData2 <- uthb$counts
expressionData2 = as.data.frame(expressionData2)
################################################################################
"%ni%" <- Negate("%in%")

# mapping emsembl ids to gene symbols
# Changed the mapping function on 07292025
mapEnsemblIDtoSymbol <- function(ensembl_id) {
  symbol <- mapIds(org.Hs.eg.db, ensembl_id, "SYMBOL", "ENSEMBL")
  return(symbol)
}

mapEntrezIDtoSymbol <- function(entrez_id) {
  symbol <- mapIds(org.Hs.eg.db, entrez_id, "SYMBOL", "ENTREZID")
  return(symbol)
}

# sum observations by column
sumByColumn <- function(x, column) {
  library(data.table)
  x <- as.data.table(x)
  x2 <- x[, lapply(.SD, sum), by = column]
  return(as.data.frame(x2))
}
# map and summarize gene expression data from ensembl to symbols
summarizeToSymbolsEnsem <- function(data) {
  # mapping to symbols
  data$symbol <- mapEnsemblIDtoSymbol(data$gene_id)
  #data$symbol <- mapEntrezIDtoSymbol(data$gene_id)
  data$gene_id <- NULL
  # summarizing by gene symbols
  data <- sumByColumn(na.omit(data), "symbol")
  rownames(data) <- data$symbol
  data$symbol <- NULL
  return(data)
}

# map and summarize gene expression data from ensembl to symbols
summarizeToSymbolsEntrez <- function(data) {
  # mapping to symbols
  #data$symbol <- mapEnsemblIDtoSymbol(data$gene_id)
  data$symbol <- mapEntrezIDtoSymbol(data$gene_id)
  data$gene_id <- NULL
  # summarizing by gene symbols
  data <- sumByColumn(na.omit(data), "symbol")
  rownames(data) <- data$symbol
  data$symbol <- NULL
  return(data)
}
################################################################################
# Adjust section for your data
# Adjust script to remove the symbol row
expressionData1$gene_id = rownames(expressionData1)
expressionData1 <- summarizeToSymbolsEnsem(expressionData1)
expressionData1$gene_id = NULL
#expressionData <- as.matrix(expressionData)


# Adjust section for your data
# Adjust script to remove the symbol row
expressionData2$gene_id = rownames(expressionData2)
expressionData2 <- summarizeToSymbolsEntrez(expressionData2)
expressionData2$gene_id = NULL
#expressionData <- as.matrix(expressionData)
################################################################################
ct_res1 = brainCells(expressionData1, nMarker = 50)
ct_res2 = brainCells(expressionData2, nMarker = 50)
################################################################################
# Using RUVseq for surrogate variables
## VABB
# Create SeqExpressionSet
set1 <- newSeqExpressionSet(as.matrix(expressionData1))
# Normalize between lanes (e.g., for GC or library size biases)
set1 <- betweenLaneNormalization(set1, which = "upper")
# Compute row-wise variance
gene_variances <- apply(normCounts(set1), 1, var)
nonzero_var_genes <- which(gene_variances > 0)
least_variable_genes <- order(gene_variances[nonzero_var_genes])[1:1000]
least_variable_genes_symbols <- rownames(set1)[nonzero_var_genes][least_variable_genes]
ruv1 <- RUVg(set1, least_variable_genes_symbols, k=10)
ruv1 <- pData(ruv1)

## UTHealth
# Create SeqExpressionSet
set2 <- newSeqExpressionSet(as.matrix(expressionData2))
# Normalize between lanes (e.g., for GC or library size biases)
set2 <- betweenLaneNormalization(set2, which = "upper")
# Compute row-wise variance
gene_variances <- apply(normCounts(set2), 1, var)
nonzero_var_genes <- which(gene_variances > 0)
least_variable_genes <- order(gene_variances[nonzero_var_genes])[1:1000]
least_variable_genes_symbols <- rownames(set2)[nonzero_var_genes][least_variable_genes]
ruv2 <- RUVg(set2, least_variable_genes_symbols, k=10)
ruv2 <- pData(ruv2)
################################################################################
# Building output objects
out1 = as.data.frame(ct_res1)
out1$name = row.names(ct_res1)
out2 = as.data.frame(ct_res2)
out2$name = row.names(ct_res2)

# SVA factors
ruv1$name = row.names(ruv1)
ruv2$name = row.names(ruv2)

# Merge
out1 = merge(out1, ruv1, by='name')
out2 = merge(out2, ruv2, by='name')

# Saving data
write.csv(file = vabb_cells,
          out1, quote = FALSE,
          row.names = FALSE)
write.csv(file = uthh_cells,
          out2, quote = FALSE,
          row.names = FALSE)