################################################################################
# Deep Neural Age ensemble (DNAe) - training
# 2020, NH
# Training of pathway-based neural networks using the hallmark gene sets for ensemble construction
################################################################################
# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})
################################################################################
## This section is UNIQUE for your dataset
# setting working directory
# Directory for output
# Set working directory
wkdir='C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/02results/04geocohorts-08032025/'
setwd(wkdir)

# loading pathway database
pathwaysPath = "C:/Users/jjm262/OneDrive - Yale University/Documents/Documents/00yale/04fourthyear/01projects/03tclock/00databases/03references/00gmt/h.all.v7.0.symbols.gmt"

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
##################################################################
# not in
"%ni%" <- Negate("%in%")

# mapping emsembl ids to gene symbols
# Changed the mapping function on 07292025
mapEnsemblIDtoSymbol <- function(ensembl_id) {
  symbol <- mapIds(org.Hs.eg.db, ensembl_id, "SYMBOL", "ENSEMBL")
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
summarizeToSymbols <- function(data) {
  # mapping to symbols
  data$symbol <- mapEnsemblIDtoSymbol(data$gene_id)
  data$gene_id <- NULL
  # summarizing by gene symbols
  data <- sumByColumn(na.omit(data), "symbol")
  rownames(data) <- data$symbol
  data$symbol <- NULL
  return(data)
}
# function for mean centering and scaling
scale_features <- function(x) {
  x <- apply(x, 2, function(y) (y - mean(y)) / sd(y))
}
# function for parsing .gmt files
parseGMT <- function (filepath) {
  res <- list(genes = list(), desc = list())
  gmt <- file(filepath)
  gmt_lines <- readLines(gmt)
  close(gmt)
  gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, split = "\t")))
  gmt_names <- sapply(gmt_list, "[", 1)
  gmt_desc <- lapply(gmt_list, "[", 2)
  gmt_genes <- lapply(gmt_list, function(x) {
    x[3:length(x)]
  })
  names(gmt_desc) <- names(gmt_genes) <- gmt_names
  res <- do.call(rbind, lapply(names(gmt_genes), function(n) cbind.data.frame(term = n, 
                                                                              gene = gmt_genes[[n]], 
                                                                              stringsAsFactors = F)))
  res$term <- as.factor(res$term)
  return(res)
}
################################################################################
# Adjust section for your data
# Adjust script to remove the symbol row
expressionData$gene_id = rownames(expressionData)
expressionData <- summarizeToSymbols(expressionData)
expressionData <- as.matrix(expressionData)


#################################################################################
# CRUCIAL: enabling compatibility with tf1
# This requieres a conda environment with the following packages
# python=3.9; tensorflow==2.12, numpy==1.24.3
library(reticulate)
use_condaenv("tf-07292025",
             conda = "C:/Users/jjm262/AppData/Local/miniconda3/Scripts/conda.exe",
             required = TRUE)
library(keras)
library(tensorflow)
#tf$compat$v1$enable_eager_execution()
# We changed this behavior
tf$compat$v1$disable_v2_behavior()
################################################################################
# Initial data loading and setup
pathways <- parseGMT(pathwaysPath)
min_n_of_genes <- 5 # not applicable for hallmark collection

# Quality check: Ensure pathways object is not empty after parsing
stopifnot("Pathways object is empty after parsing." = nrow(pathways) > 0)
message(paste0("Initial pathways loaded: ", length(unique(pathways$term)), " terms, ", nrow(pathways), " gene-term associations."))

# cleaning pathway descriptions
pathways$term <- make.names(pathways$term)
# Quality check: Ensure no NA terms introduced after make.names
stopifnot("NA terms found after make.names." = !anyNA(pathways$term))
message("Pathway terms cleaned.")

# preparing data
x <- t(expressionData)
# Quality check: Ensure x has dimensions and no NA rows before log transformation
stopifnot("Expression data 'x' is empty or invalid." = !is.null(x) && all(dim(x) > 0))
message(paste0("Expression data dimensions after transpose: ", paste(dim(x), collapse = "x")))

# log transformation
# Note: na.omit is applied row-wise, ensuring no NA rows proceed
x <- t(na.omit(apply(x, 1, function(x) log(x+1))))
# Quality check: Ensure x still has dimensions and no NA after log transform and na.omit
stopifnot("Expression data 'x' is empty or invalid after log transformation and NA omission." = !is.null(x) && all(dim(x) > 0))
message(paste0("Expression data dimensions after log transformation and NA omission: ", paste(dim(x), collapse = "x")))

# filtering genes
x <- x[,colMeans(x) > 0]
# Quality check: Ensure x still has columns after filtering based on colMeans
stopifnot("No genes remaining after filtering by colMeans > 0." = ncol(x) > 0)
message(paste0("Genes filtered by colMeans > 0. Remaining genes: ", ncol(x)))
x <- x[,colnames(x) %in% pathways$gene]
# Quality check: Ensure x still has columns after filtering based on presence in pathways
stopifnot("No genes remaining after filtering by presence in pathways." = ncol(x) > 0)
message(paste0("Genes filtered by presence in pathways. Remaining genes: ", ncol(x)))

# final shape
message(paste0("Final shape of expression data 'x' before scaling: ", paste(dim(x), collapse = "x")))
# Quality check: Confirm final dimensions of x
stopifnot("Final expression data 'x' has invalid dimensions." = all(dim(x) > 0))

# scaling the data
x <- scale_features(x)
# Quality check: Ensure x remains numeric and has no NaNs/Infs after scaling
stopifnot("Expression data 'x' is not numeric or contains NaNs/Infs after scaling." = is.numeric(x) && !anyNA(x) && !any(is.infinite(x)))
message("Expression data scaled.")

# subsetting to include only pathways with more than n genes
pathways <- pathways[pathways$gene %in% colnames(x),]
# Quality check: Ensure pathways object is not empty after subsetting based on genes in x
stopifnot("Pathways object is empty after filtering by genes present in 'x'." = nrow(pathways) > 0)
message(paste0("Pathways filtered by genes present in expression data. Remaining gene-term associations: ", nrow(pathways)))

sizes <- pathways %>% group_by(term) %>% count(term)
pathways <- pathways[pathways$term %in% sizes$term[sizes$n > min_n_of_genes],]
# Quality check: Ensure pathways object is not empty after filtering by min_n_of_genes
# Corrected Quality Check:
if (nrow(pathways) == 0) {
  stop(paste0("No pathways remaining after filtering for > ", min_n_of_genes, " genes. The 'pathways' object became empty."))
}
message(paste0("Pathways filtered by minimum number of genes (>", min_n_of_genes, "). Remaining unique pathways: ", length(unique(pathways$term))))

# preparing pathway filter matrix
num_unique_pathways <- length(unique(pathways$term))
num_x_cols <- ncol(x)
pathwayFilter <- as.data.frame(matrix(nrow = num_unique_pathways, ncol = num_x_cols))
rownames(pathwayFilter) <- unique(pathways$term)
colnames(pathwayFilter) <- colnames(x)
# Quality check: Confirm dimensions of pathwayFilter matrix
stopifnot("Pathway filter matrix has incorrect dimensions." = nrow(pathwayFilter) == num_unique_pathways && ncol(pathwayFilter) == num_x_cols)
message(paste0("Initialized pathwayFilter matrix with dimensions: ", paste(dim(pathwayFilter), collapse = "x")))

# filling the filter matrix
index <- 1
for(i in rownames(pathwayFilter)) {
  pathwayFilter[i,] <- as.numeric(colnames(pathwayFilter) %in% pathways$gene[pathways$term == i])
  message(paste0(index, " / ", length(rownames(pathwayFilter)), " pathways built"))
  
  # Corrected Quality check inside loop: Ensure at least one gene is set to 1 for each pathway
  if (sum(pathwayFilter[i,]) == 0) {
    stop(paste0("Pathway '", i, "' has no genes marked (all zeros) in filter matrix. This pathway might be empty or its genes are not in the expression data."))
  }
  index <- index + 1
}
message("Pathway filter matrix filled.")

# subsetting and checking gene orders
pathwayFilter <- pathwayFilter[,colnames(pathwayFilter) %in% colnames(x)]
# Quality check: Ensure pathwayFilter still has columns
stopifnot("pathwayFilter lost all columns after subsetting by colnames(x)." = ncol(pathwayFilter) > 0)
x <- x[,colnames(x) %in% colnames(pathwayFilter)]
# Quality check: Ensure x still has columns
stopifnot("x lost all columns after subsetting by colnames(pathwayFilter)." = ncol(x) > 0)
pathwayFilter <- pathwayFilter[,order(match(colnames(pathwayFilter), colnames(x)))]
# Quality check: Ensure sorting didn't introduce issues
stopifnot("Sorting pathwayFilter columns resulted in empty matrix." = ncol(pathwayFilter) > 0)
message("Pathway filter matrix and expression data columns aligned.")

# final sanity check
final_sanity_check <- all(colnames(pathwayFilter) == colnames(x))
stopifnot("Final sanity check failed: Column names of pathwayFilter and x do not match exactly." = final_sanity_check)
message("Final sanity check passed: Column names of pathwayFilter and x match.")

# setting seed for reproducibility of train/test split
set.seed(2)
message("Random seed set for reproducibility.")

# creating train/test split
# Quality check: Ensure phenoData is not empty and has zz_nr column
stopifnot("phenoData is empty or missing 'zz_nr' column." = !is.null(phenoData$zz_nr) && length(phenoData$zz_nr) > 0)

train_ids <- sample(phenoData$zz_nr, size = floor(nrow(phenoData)*0.7))
test_ids <- phenoData$zz_nr[phenoData$zz_nr %ni% train_ids]
# Quality check: Ensure train and test IDs are not empty and don't overlap
stopifnot("Train IDs are empty." = length(train_ids) > 0)
stopifnot("Test IDs are empty." = length(test_ids) > 0)
stopifnot("Train and test IDs overlap." = length(intersect(train_ids, test_ids)) == 0)
stopifnot("Total train and test IDs do not sum to total phenoData rows." = (length(train_ids) + length(test_ids)) == nrow(phenoData))
message(paste0("Train/test split created. Train samples: ", length(train_ids), ", Test samples: ", length(test_ids)))

# splitting data matrices
x_train <- x[rownames(x) %in% train_ids,]
x_test <- x[rownames(x) %in% test_ids,]
y_train <- phenoData$intro_age[phenoData$zz_nr %in% train_ids]
y_test <- phenoData$intro_age[phenoData$zz_nr %in% test_ids]

# Quality checks for split data
stopifnot("x_train is empty or has incorrect dimensions." = all(dim(x_train) > 0) && nrow(x_train) == length(train_ids) && ncol(x_train) == ncol(x))
stopifnot("x_test is empty or has incorrect dimensions." = all(dim(x_test) > 0) && nrow(x_test) == length(test_ids) && ncol(x_test) == ncol(x))
stopifnot("y_train is empty or has incorrect length." = length(y_train) == length(train_ids) && !anyNA(y_train))
stopifnot("y_test is empty or has incorrect length." = length(y_test) == length(test_ids) && !anyNA(y_test))

# Quality check: Ensure rownames of x_train/x_test match corresponding phenoData zz_nr and y_train/y_test
stopifnot("Row names of x_train do not match train_ids." = all(sort(rownames(x_train)) == sort(as.character(train_ids))))
stopifnot("Row names of x_test do not match test_ids." = all(sort(rownames(x_test)) == sort(as.character(test_ids))))
stopifnot("Length of y_train does not match nrow(x_train)." = length(y_train) == nrow(x_train))
stopifnot("Length of y_test does not match nrow(x_test)." = length(y_test) == nrow(x_test))

message("Data successfully split into train and test sets.")
################################################################################
# Setting model options
ensemble_size <- 10
input <- ncol(x)
regularization <- 0.01
dropout <- 0.1
# training options
alpha <- 0.4
loss <- "mean_squared_error"
metric <- "mean_absolute_error" 
batch_size <- 16
batch_size_prediction <- 512
epochs <- 200
lr <- 1e-3

################################################################################
# function performing model construction with quality checks
build_model <- function(inputTensor, model_id) {
  
  # Ensure global variables are defined
  required_globals <- c("pathwayFilter", "dropout", "regularization", "lr", "alpha", "metric")
  missing_globals <- required_globals[!sapply(required_globals, exists, inherits = TRUE)]
  if (length(missing_globals) > 0) {
    stop("Missing global variables: ", paste(missing_globals, collapse = ", "))
  }
  
  # Check input tensor
  if (is.null(inputTensor)) stop("inputTensor cannot be NULL")
  
  # Check pathwayFilter
  if (!is.data.frame(pathwayFilter) && !is.matrix(pathwayFilter)) {
    stop("pathwayFilter must be a data.frame or matrix")
  }
  
  # defining function for targeted connections
  customConnection <- function(tensor, indices) {
    if (length(indices) == 0) stop("Empty index list in customConnection()")
    tensor <- tensorflow::tf$gather(tensor, indices, axis = as.integer(1))
    return(tensor)
  }
  
  # defining loss function with return
  customLossWrappper <- function(alpha) {
    force(alpha)
    return(function(y_true, y_pred) {
      main_loss <- k_mean(k_square(y_pred[[1]] - y_true[[1]]))
      aux_loss <- k_mean(k_square(y_pred[[2]] - y_true[[2]]))
      loss <- (1 - alpha) * main_loss + alpha * aux_loss
      return(loss)
    })
  }
  
  # Initialize variables
  pind <- 1
  pnumb <- nrow(pathwayFilter)
  layers <- list()
  
  # build model
  for (i in rownames(pathwayFilter)) {
    selected <- pathwayFilter[pind, ] %>% as.logical()
    
    # skip empty or all-TRUE pathways
    if (all(selected) || !any(selected)) {
      message(paste0("Skipping ", i, ": empty or full gene set"))
      pind <- pind + 1
      next
    }
    
    gind <- which(selected)
    gind <- as.integer(gind - 1)  # zero-based indexing
    message(paste0(i, ": ", length(gind), " genes"))
    
    # Create lambda layer for sparse connection
    tensor <- layer_lambda(
      f = function(x) customConnection(tensor = x, indices = as.list(gind)),
      name = paste0(i, "_filter_", model_id)
    )(inputTensor)
    
    for (d in 1:4) {
      tensor <- layer_dense(units = 5 + ceiling(length(gind) * 0.5),
                            kernel_initializer = initializer_he_uniform(),
                            kernel_regularizer = regularizer_l2(regularization),
                            name = paste0(i, "_dense", d, "_", model_id))(tensor)
      tensor <- layer_activation_elu(name = paste0(i, "_activation", d, "_", model_id))(tensor)
      tensor <- layer_dropout(rate = dropout, name = paste0(i, "_dropout", d, "_", model_id))(tensor)
    }
    
    # Final pooling
    tensor <- layer_dense(units = 1,
                          kernel_initializer = initializer_he_uniform(), 
                          kernel_regularizer = regularizer_l2(regularization),
                          name = paste0(i, "_pooling_", model_id))(tensor)
    
    layers[[pind]] <- tensor
    message(paste0(pind, " / ", pnumb, " pathways connected"))
    pind <- pind + 1
  }
  
  if (length(layers) == 0) stop("No valid pathway layers were created.")
  
  # Auxiliary output
  AuxOutput <- layer_concatenate(layers, name = paste0("aux_output_", model_id))
  
  # Final model output
  OutputTensor <- layer_concatenate(layers, name = paste0("concatenate_", model_id))
  OutputTensor <- layer_batch_normalization(name = paste0("normalization_", model_id))(OutputTensor)
  OutputTensor <- layer_dropout(rate = dropout, name = paste0("dropout_", model_id))(OutputTensor)
  OutputTensor <- layer_dense(units = 1, 
                              kernel_regularizer = regularizer_l2(regularization),
                              kernel_initializer = initializer_zeros(),
                              kernel_constraint = constraint_nonneg(),
                              name = paste0("output_", model_id))(OutputTensor)
  
  # Create model
  model <- keras_model(inputTensor, output = list(OutputTensor, AuxOutput))
  model %>% summary()
  
  # Compile model
  model %>% compile(
    loss = customLossWrappper(alpha),
    optimizer = optimizer_adam(learning_rate = lr),
    metrics = metric
  )
  
  message(paste0("Model ", model_id, " construction and compilation finished"))
  return(model)
}

################################################################################
# creating the input layer
inputTensor <- layer_input(shape = c(input), name = "input")

if (ncol(pathwayFilter) != input) {
  stop("Mismatch: pathwayFilter columns (", ncol(pathwayFilter), ") != input shape (", input, ")")
}

# training models for the ensembl
for(model_id in 1:ensemble_size) {
  # buidling the model
  model <- build_model(inputTensor, model_id)
  # training the model
  history <- model %>% fit(x = x_train,
                           y = list(y_train, y_train), 
                           epochs = epochs, 
                           batch_size = batch_size, 
                           validation_data = list(x_test, list(y_test, y_test)),
                           verbose = 1)
  message(paste0("model ", model_id, " training finished"))
  # saving model weights
  model %>% save_model_weights_tf(paste0("model_", model_id))
  message(paste0("model ", model_id, " weights saved"))
  # saving training history
  history_df <- as.data.frame(history)
  write.table(history_df, paste0("training_history_model_", model_id, ".tsv"), row.names = F, sep = "\t", quote = F)
  # saving model parameters
  model_parameters <- data.frame(n_pathways = nrow(pathwayFilter),
                                 n_input_params = input,
                                 n_parameters = model$count_params(),
                                 activation_function = "elu",
                                 weight_initialization = "he",
                                 regularization = regularization,
                                 dropout = dropout,
                                 alpha = alpha,
                                 loss = loss,
                                 metric = metric,
                                 lr = lr,
                                 batch_size = batch_size,
                                 epochs = epochs) %>% t() %>% as.data.frame()
  write.table(model_parameters, file = paste0("parameters_model_", model_id, ".tsv"), row.names = T, col.names = F, sep = "\t", quote = F)
}