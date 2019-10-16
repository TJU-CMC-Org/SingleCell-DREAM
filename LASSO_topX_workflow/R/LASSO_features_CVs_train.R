########################################
# Make sure you change the paths 
# of the input files to rerun the script 
########################################

rm(list = ls())
set.seed(1234)
require(glmnet)
require(data.table)
source("LASSO_topX_workflow/R/glmnetExtensionLibrary.R")

# ***************
# USER input ####
# ***************
numFeaturesIn <- c(20, 40, 60)
repCV <- 20
nfoldsCV <- 5

# Data to use
# NOTE: Change the paths of the input files to rerun the script
singleCellRNAseq <- fread(input = "Data/dge_normalized.txt.gz", sep = "\t")
colnames(singleCellRNAseq)[1] <- "Genes"

# insitu data
# NOTE: Change the paths of the input files to rerun the script
insitu.matrix = read.csv("Data/binarized_bdtnp.csv.gz",check.names = FALSE)

# Coordinates to predict 
load(file = "Results_Common/geometry.RNAseq.RData")
# Keep cells which are uniquelly mapped
rm(geometry.RNAseq.all)

# Use inSitu genes only
useOnlyInSitu <- TRUE

# Read provided training index - 10 fold CV
folds_train <- fread(input = "Data/CV_folds/folds_train.csv")
folds_test <- fread(input = "Data/CV_folds/folds_test.csv")
cell_ids <- fread(input = "Data/CV_folds/cell_ids.csv")

# Manually set lambda grid
grid = 10^seq(from = 1.3, to = -1, length.out = 300)

# Path to folder to save results
path2saveResults <- "Modified_LASSO_workflow/Results_UniquelyMapped_cells_inSituRNAseq"
# path2saveResults <- "Modified_LASSO_workflow/Results_UniquellyMapped_cells_allRNASeq"

# *****************
# CALCULATIONS ####
# *****************

# Filter with insitu genes
if(useOnlyInSitu){
    cat("Using inSitu genes...\n")
    singleCellRNAseq <- singleCellRNAseq[Genes %in% colnames(insitu.matrix)]    
}else{
    cat("Using all genes...\n")    
}


# transpose singleCellRNAseq
singleCellRNAseq.matrix <- as.matrix(singleCellRNAseq[,-1])
rownames(singleCellRNAseq.matrix) <- singleCellRNAseq$Genes
singleCellRNAseq.matrix <- t(singleCellRNAseq.matrix)

singleCellRNAseq.dt <- data.table(singleCellRNAseq.matrix, keep.rownames = TRUE)
colnames(singleCellRNAseq.dt)[1] <- "Cell"
# Align features with outcomes
all.data.dt <- merge(x = geometry.RNAseq.unique, y = singleCellRNAseq.dt, by = "Cell")
geometry.RNAseq <- all.data.dt[, c("xcoord", "ycoord", "zcoord", "Cell")]

# Reformat singleCell matrix to use with LASSO
singleCellRNAseq.matrix <- as.matrix(all.data.dt[, -c(1:4)])
dim(singleCellRNAseq.matrix)
rownames(singleCellRNAseq.matrix) <- all.data.dt$Cell


# Execute cv folds
for(fold_i in 1:nrow(folds_train)){
    
    cat("CV", fold_i, "out of", nrow(folds_train), "...\n")
    
    # Cells to train
    cells_train_fold_i <- cell_ids[cell_ids$id %in% as.numeric(folds_train[fold_i, ])]
    # Cells to test
    cells_test_fold_i <- cell_ids[!cell_ids$id %in% as.numeric(folds_train[fold_i, ])]
    
    # Training Sets
    singleCellRNAseq.train <- singleCellRNAseq.matrix[
        rownames(singleCellRNAseq.matrix) %in% cells_train_fold_i$cell, ]
    geometry.RNAseq.train <- geometry.RNAseq[
        match(x = rownames(singleCellRNAseq.train), table = geometry.RNAseq$Cell), ]
    if(!all(rownames(singleCellRNAseq.train) == geometry.RNAseq.train$Cell)){
        stop("training sets are not aligned")
    }
    geometry.RNAseq.train <- as.matrix(geometry.RNAseq.train[, -4])
    
    # Test set
    singleCellRNAseq.test <- singleCellRNAseq.matrix[
        rownames(singleCellRNAseq.matrix) %in% cells_test_fold_i$cell, ]
    geometry.RNAseq.test <- geometry.RNAseq[
        match(x = rownames(singleCellRNAseq.test), table = geometry.RNAseq$Cell), ]
    if(!all(rownames(singleCellRNAseq.test) == geometry.RNAseq.test$Cell)){
        stop("test sets are not aligned")
    }
    geometry.RNAseq.test <- as.matrix(geometry.RNAseq.test[, -4])
    
    # Identify most stable / important inSitu genes
    cv.lasso.mse <- LASSO.topX(x = singleCellRNAseq.train, y = geometry.RNAseq.train,
                                 numFeatures = numFeaturesIn, lambda = grid, nfolds = nfoldsCV, 
                                 reps = repCV, ncores = detectCores()/2)
    
    # Save results
    filename_i <- paste0(path2saveResults, "/FeaturesInformation_CV_",
                         fold_i, 
                         "_Data_UniquellyMappedCells_LASSO_Reps20_Folds5.RData")
    save(singleCellRNAseq.train, geometry.RNAseq.train, 
         singleCellRNAseq.test, geometry.RNAseq.test,
         cv.lasso.mse, file = filename_i)
    
    # Remove variables
    rm(singleCellRNAseq.train, geometry.RNAseq.train, 
       singleCellRNAseq.test, geometry.RNAseq.test,
       cv.lasso.mse)
    
    gc()

}

