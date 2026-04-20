##################################################################
# Estimate RNAAgeCalc age in GEO cohorts
# Author: Jose Jaime Martinez-Magana
# Day: 08052025
##################################################################
# Load libraries
library(RNAAgeCalc)
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
##################################################################
# Set working directory
wkdir='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/04geocohorts-08032025/'
setwd(wkdir)
##################################################################
# Load datasets
# GSE102556
fpkm_102556 <- read.delim("GSE102556_HumanMDD_fpkmtab.txt.gz", stringsAsFactors = FALSE)
pheno_102556 <- read.csv("GSE102556_phenoData.csv", stringsAsFactors = FALSE)
# Subset rows with "BA" in the title
pheno_ba <- pheno_102556[grepl("BA", pheno_102556$title), ]

# Extract number before colon
pheno_ba$sample_num <- sub(":.*", "", pheno_ba$title)
# Determine region and assign BA label
pheno_ba$region_id <- ifelse(
  grepl("OFC", pheno_ba$title), 
  paste0(pheno_ba$sample_num, ".BA11"),
  ifelse(
    grepl("dlPFC", pheno_ba$title),
    paste0(pheno_ba$sample_num, ".BA8_9"),
    NA  # fallback if neither matches
  )
)
# Determine region and assign BA label
pheno_ba$region_id <- paste0("S", pheno_ba$region_id)
##################################################################
# Add expression
# Add "S" to all column names except "gene_id"
colnames(fpkm_102556)[-1] <- paste0("S", colnames(fpkm_102556)[-1])
colnames(fpkm_102556) <- gsub("^SX", "S", colnames(fpkm_102556))
fpkm_102556 <- fpkm_102556[, -((ncol(fpkm_102556) - 2):ncol(fpkm_102556))]
# Set gene_id as row names
rownames(fpkm_102556) <- fpkm_102556$gene_id
# Remove gene_id column
fpkm_102556$gene_id <- NULL
exprdataframe <- fpkm_102556[, grepl("BA", colnames(fpkm_102556))]

# Change the annotation and datatype if needed
# Make sure it matches your input data
#annot = "ENTREZID"
annot = "ENSEMBL"

# Make sure it matches your input data
#exprtype = 'counts'
exprtype = 'FPKM'

# Warning DO NOT modify the following code
### Estimate RNAAgeCalc
# This part of the script takes the counts of a feature counts object and 
# estimates all the signatures and tissues in the RNAAgeCalc and gives back a csv file
# Define your tissue and signature options
tissues <- c("adipose_tissue", "adrenal_gland", "blood", "blood_vessel", "brain",
             "breast", "colon", "esophagus", "heart", "liver", "lung", "muscle",
             "nerve", "ovary", "pancreas", "pituitary", "prostate", "salivary_gland",
             "skin", "small_intestine", "spleen", "stomach", "testis", "thyroid",
             "uterus", "vagina")

signatures <- c("DESeq2", "Pearson", "Dev", "deMagalhaes", "GenAge", "GTExAge", "Peters", "all")

# Store results
age_results <- list()

# Loop through tissues and signatures
for (tissue in tissues) {
  for (signature in signatures) {
    message(paste("Predicting for:", tissue, signature))
    
    pred <- tryCatch({
      predict_age(exprdata = exprdataframe,
                  tissue = tissue,
                  exprtype = exprtype,
                  idtype = annot,
                  signature = signature)
    }, error = function(e) {
      warning(paste("Failed for", tissue, signature, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(pred)) {
      pred_numeric <- unlist(pred)
      df <- data.frame(
        sample = names(pred_numeric),
        predicted_age = as.numeric(pred_numeric),
        tissue = tissue,
        signature = signature,
        stringsAsFactors = FALSE
      )
      age_results[[paste(tissue, signature, sep = "_")]] <- df
    }
  }
}

# Combine into one dataframe
all_predictions <- do.call(rbind, age_results)

# Wide format (optional)
age_matrix <- pivot_wider(all_predictions,
                          id_cols = sample,
                          names_from = c(signature, tissue),
                          values_from = predicted_age,
                          names_sep = "_")
# Preview
age_matrix$sample = colnames(exprdataframe)

# Extract relevant metadata from characteristics
pheno_clean <- pheno_ba %>%
  mutate(
    age = str_extract(characteristics_ch1.3, "\\d+"),
    gender = str_extract(characteristics_ch1.1, "(?<=gender: )\\w+"),
    rin = as.numeric(str_extract(characteristics_ch1.13, "[0-9.]+")),
    phenotype = str_extract(characteristics_ch1.6, "(?<=phenotype: )\\w+"),
    cause_of_death = str_extract(characteristics_ch1.2, "(?<=Cause of death: )[^;]+"),
    pmi = as.numeric(str_extract(characteristics_ch1.4, "[0-9.]+"))
  ) %>%
  select(X, title, geo_accession, age, gender, rin, phenotype, cause_of_death, pmi)

# 
pheno_clean$sample <- pheno_ba$region_id
pheno_clean$X = NULL
pheno_clean$title = NULL

# Merge
merged = merge(age_matrix, pheno_clean, by='sample')

# Save results
outfile = "00gse102556-08052025.csv"
fwrite(merged,
       file = outfile,
       row.names = FALSE,
       quote = FALSE)

# Create a list to save expression and phenodata
out = list()
out$expressionData = exprdataframe
out$phenoData = pheno_clean
outRDS = "00gse102556-08052025.rds"
saveRDS(out, file = outRDS)
##################################################################
