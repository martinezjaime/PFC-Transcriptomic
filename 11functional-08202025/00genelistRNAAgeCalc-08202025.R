##################################################
# Functional Gene Analysis
##################################################
# Load libraries
# Load package
library(RNAAgeCalc)
library(org.Hs.eg.db)
library(AnnotationDbi)
##################################################
# Set working directory
wkdir = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/08functional"
setwd(wkdir)
##################################################
# Access internal objects
genelist_all <- RNAAgeCalc:::genelist_all
genelist_cau <- RNAAgeCalc:::genelist_cau

# Now run your extraction function
extract_signatures <- function(genelist) {
  df_list <- list()
  
  for (tissue in names(genelist)) {
    for (signature in names(genelist[[tissue]])) {
      coefs <- genelist[[tissue]][[signature]]
      df <- data.frame(
        Tissue = tissue,
        Signature = signature,
        Gene = names(coefs)[-1],      # exclude intercept
        Coefficient = coefs[-1],      # corresponding weights
        stringsAsFactors = FALSE
      )
      df_list[[paste(tissue, signature, sep = "_")]] <- df
    }
  }
  do.call(rbind, df_list)
}

# Extract all
df_all <- extract_signatures(genelist_all)
df_cau <- extract_signatures(genelist_cau)

# Map SYMBOL -> ENSEMBL
ensembl_ids <- mapIds(
  org.Hs.eg.db,
  keys = unique(df_all$Gene),
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"   # in case multiple mappings
)

# Add as a new column
df_all$Ensembl <- ensembl_ids[df_all$Gene]
##############################################################
# Save output
write.csv(file="00geneRNAAgeCalc-08202025.csv",
          df_all,
          quote = FALSE,
          row.names = FALSE)