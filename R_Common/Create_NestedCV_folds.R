# Create random Nested CV folds #
rm(list = ls())
library(data.table)
library(foreach)
cat("Creating Random Outer folds for Nested CV...\n")

format_train <- fread(file = "../../src_myPubl/Data/CV_folds/folds_train.csv")
format_train[, (ncol(format_train)-3) : ncol(format_train), with = FALSE]

format_test <- fread(file = "../../src_myPubl/Data/CV_folds/folds_test.csv")
format_cell_ids <- fread(file = "../../src_myPubl/Data/CV_folds/cell_ids.csv")


# ***************
# USER input ####
# ***************
numFolds <- 10
folder2save <- "Data/CV_folds"

# *****************
# CALCULATIONS ####
# *****************
set.seed(123)

# Read Normalized Expression to get cells names
expression <- fread(file = "Data/dge_normalized.txt.gz")
cells <- colnames(expression)[-1]

# Randomly assign each cell to a fold
folds <- sample(x = 1:numFolds, size = length(cells), replace = TRUE)
max(table(folds)) - min(table(folds))

# Extract train test IDs
trainIDs <- foreach(fold_i = 1:numFolds) %do% {
    return(as.data.table(t(matrix(which(folds != fold_i)))))
}
trainIDs <- rbindlist(l = trainIDs, fill = TRUE)
testIDs <- foreach(fold_i = 1:numFolds) %do% {
    return(as.data.table(t(matrix(which(folds == fold_i)))))
}
testIDs <- rbindlist(l = testIDs, fill = TRUE)

# Cell IDs
cellIDs <- data.table(id = 1:length(cells), 
                      cell = cells)
# Write files
fwrite(x = trainIDs, file = paste0(folder2save, "/folds_train.csv"))
fwrite(x = testIDs, file = paste0(folder2save, "/folds_test.csv"))
fwrite(x = cellIDs, file = paste0(folder2save, "/cell_ids.csv"))





















