# Set working directory
wkdir = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/04geocohorts-08032025/"
setwd(wkdir)

# Load libraries
library(GEOquery)

# Download the GEO dataset (ExpressionSet)
gse <- getGEO("GSE211792", GSEMatrix = TRUE)
# Sample metadata (phenotypic data)
gse_data <- gse[[1]]
phenoData <- pData(gse_data)
write.csv(phenoData, "GSE211792_phenoData.csv")

# Download expression counts
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE211792&format=file&file=GSE211792%5Freadcount%5Fmatrix%2Etxt%2Egz",
  destfile = "GSE211792_featurecounts_rawcounts.txt.gz"
)
# Load them
counts <- read.delim(gzfile("GSE211792_featurecounts_rawcounts.txt.gz"), row.names = 1)
# Check shapes
dim(counts)     # Should be genes x samples
dim(phenoData)   # Should be samples x covariates


# Download the GEO dataset (ExpressionSet)
gse <- getGEO("GSE102556", GSEMatrix = TRUE)
# Sample metadata (phenotypic data)
gse_data <- gse[[1]]
phenoData <- pData(gse_data)
write.csv(phenoData, "GSE102556_phenoData.csv")
# Download expression counts
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102556&format=file&file=GSE102556%5FHumanMDD%5Ffpkmtab%2Etxt%2Egz",
  destfile = "GSE102556_HumanMDD_fpkmtab.txt.gz"
)
# Load them
counts <- read.delim(gzfile("GSE102556_HumanMDD_fpkmtab.txt.gz"), row.names = 1)
# Check shapes
dim(counts)     # Should be genes x samples
dim(phenoData)   # Should be samples x covariates
