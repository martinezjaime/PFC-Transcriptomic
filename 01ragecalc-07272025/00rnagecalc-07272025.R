##################################################################
# Estimate RNAAgeCalc age
# Author: Jose Jaime MArtinez-Magana
# Day: 07272025
##################################################################
# Load libraries
library(RNAAgeCalc)
library(tidyr)
library(data.table)
##################################################################
# Set working directory
wkdir='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/01scripts/01ragecalc-07272025/'
setwd(wkdir)
##################################################################
# Load datasets
uthbP = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01transcriptome/00data/matrix_counts_featurecounts-07252025-03.rds"
vabbP = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/01transcriptome/00data/matrix_counts_featurecounts-07252025-02.rds"

uthb = readRDS(uthbP)
vabb = readRDS(vabbP)
##################################################################
# Fixing names of samples
colnames(uthb$counts) <- sub("_", "", sub("_Aligned.*", "", colnames(uthb$counts)))

# Apply to your actual column names vector
new_colnames <- sub("^((Sample_[^_]+)).*", "\\1", colnames(vabb$counts))
new_colnames <- sub("Aligned.sortedByCoord.out.bam", "", new_colnames)
new_colnames <- sub("_r", "", new_colnames)
new_colnames <- sub("_", "", new_colnames)

# Optionally assign them back
colnames(vabb$counts) <- new_colnames
##################################################################
# Set the input object
exprdataframe = uthb$counts
exprdataframe = vabb$counts
# Set the path for the output object
gzfile_path = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/02trans_age-07272025/00uthhealth_rnaagecalc-07272025.csv"
gzfile_path = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/02trans_age-07272025/01vabb_rnaagecalc-07272025.csv"

# Change the annotation if needed
# Make sure it matches your input data
annot = "ENTREZID"
annot = "ENSEMBL"

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
                  exprtype = "counts",
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
# Save results
fwrite(age_matrix, file = gzfile_path, row.names = FALSE)
##################################################################
# End of Script
