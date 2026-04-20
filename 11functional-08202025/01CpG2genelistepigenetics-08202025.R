##################################################
# Functional of CpG Gene Analysis
##################################################
### Load libraries
library(dplyr)
library(stringr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

### Set working directory
wkdir = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/08functional/00epigenetics/"
setwd(wkdir)

# -------------------------------
# Load annotation packages
# -------------------------------
annos450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annosEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# -------------------------------
# List all CSVs in your directory
# -------------------------------
data_dir <- "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/03references/03epiclocks"
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
cpg_files <- csv_files

# -------------------------------
# Start Annotation
# -------------------------------
all_annotations <- data.frame()
for (file in cpg_files) {
  clock_name <- tools::file_path_sans_ext(basename(file))
  
  # Read CpG IDs
  df <- read.csv(file, stringsAsFactors = FALSE)
  cpgs <- df[[1]]   # assumes CpGs are in the first column
  
  # Try 450K first
  gene_info <- annos450k[rownames(annos450k) %in% cpgs,
                         c("Name", "UCSC_RefGene_Name", "chr", "pos")]
  
  # If no matches, try EPIC
  if (nrow(gene_info) == 0) {
    gene_info <- annosEPIC[rownames(annosEPIC) %in% cpgs,
                           c("Name", "UCSC_RefGene_Name", "chr", "pos")]
  }
  
  if (nrow(gene_info) > 0) {
    # For CpGs with no gene symbol, replace with chr:pos
    gene_symbols <- gene_info$UCSC_RefGene_Name
    missing_genes <- gene_symbols == "" | is.na(gene_symbols)
    gene_symbols[missing_genes] <- paste0(gene_info$chr[missing_genes], ":", gene_info$pos[missing_genes])
    
    res <- data.frame(
      Clock = rep(clock_name, nrow(gene_info)),
      CpG   = rownames(gene_info),
      Gene  = gene_symbols,
      stringsAsFactors = FALSE
    )
    all_annotations <- rbind(all_annotations, res)
  } else {
    message("⚠️ No CpGs matched annotation for clock: ", clock_name)
  }
}

# -------------------------------
# Read Ensembl mapping
# -------------------------------
map_file <- file.path(data_dir, "ensembl_to_gene.csv")  # adjust path/filename if needed
map_df <- read.csv(map_file, stringsAsFactors = FALSE)
colnames(map_df) <- c("EnsemblID", "GeneSymbol")  # ensure consistent names

# -------------------------------
# Merge Ensembl IDs into all_annotations
# -------------------------------
# Function to map multiple gene symbols to Ensembl IDs
map_genes_to_ensembl <- function(gene_str, map_df) {
  # If the gene looks like chr:pos, return NA
  if (grepl("^chr\\:", gene_str)) return(NA_character_)
  
  # Split multiple genes by ";" and remove duplicates
  genes <- unique(unlist(strsplit(gene_str, ";")))
  
  # Lookup Ensembl IDs
  ensembl_ids <- map_df$EnsemblID[match(genes, map_df$GeneSymbol)]
  
  # Collapse back into semicolon-separated string
  paste(ensembl_ids, collapse = ";")
}
# Apply to all_annotations
all_annotations$EnsemblID <- sapply(all_annotations$Gene, map_genes_to_ensembl, map_df = map_df)

# -------------------------------
# Save combined results
# -------------------------------
write.csv(all_annotations, "01EpigeneticCLocks_CpG_GeneAnnotations-08202025.csv", row.names = FALSE)
