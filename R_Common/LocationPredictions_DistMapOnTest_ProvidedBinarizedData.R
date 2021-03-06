#################################################
# SUMMARY
# This script subsets the insitu matrix using the 
# genes selected by LASSO or Deep Neural Networks 
# and predicts the cell locations using distMap. 
# It performs this task in the Cross Validation folds
# provided.

# Make sure you change the paths 
# of the input files to rerun the script 
#################################################


#################################################
# Set everything up to run distMap 
# This part of the script was provided online 
#################################################

rm(list = ls())
library(DistMap)
library(foreach)
library(doMC)
library(data.table)
set.seed(1234)

# Import raw and normalised data decompress the downloaded dataset and import to data.frame.
# Important: quote is set to empty string, this way we can prevent the mis-interpretation of a gene name, that contains an apostrophe

# We convert the data.frame to matrix: each column represents a cell and each row corresponds to a gene.
raw.data = read.table("Data/dge_raw.txt.gz", sep = "\t", row.names = NULL,
                      stringsAsFactors = F, quote = "")
raw.data.genes = raw.data$V1
raw.data$V1 = NULL

# Let’s fix the gene names that contains apostrophe – this would generate issues
# gene names with apostrophes
# print(grep("'",raw.data.genes,value = T,fixed = T))
# [1] "beta'COP" "PP2A-B'"
raw.data.genes = gsub("'","",raw.data.genes,fixed = T)

raw.data = as.matrix(raw.data)
rownames(raw.data) = raw.data.genes

# Repeat for the normalised data
normalized.data = read.table("Data/dge_normalized.txt.gz", sep = "\t",
                             row.names = NULL, stringsAsFactors = F, quote = "")
normalized.data.genes = normalized.data$row.names
normalized.data$row.names = NULL

# gene names with apostrophes
# print(grep("'",normalized.data.genes,value = T,fixed = T))

normalized.data.genes = gsub("'","",normalized.data.genes,fixed = T)

normalized.data = as.matrix(normalized.data)
rownames(normalized.data) = normalized.data.genes

# MY ASSUMPTION #
colnames(raw.data) <- colnames(normalized.data)

# Check that the gene names are identical in the raw and normalised dataset
stopifnot(all(normalized.data.genes == raw.data.genes))

# Import in situ datasets
insitu.matrix = read.table("Data/binarized_bdtnp.csv.gz", sep = ",",header = T)

## Warning in read.table(gzfile("binarized_bdtnp.csv.gz", "rt"), sep = ",", :
## seek on a gzfile connection returned an internal error

## Warning in read.table(gzfile("binarized_bdtnp.csv.gz", "rt"), sep = ",", :
## seek on a gzfile connection returned an internal error

insitu.genes_orig <- colnames(insitu.matrix)

# Match the gene names across datasets

# The following few lines of code checks if there is any mismatch between the gene names.

# We will find that 2 gene names were changed during the import of the dataset and the - (dash) character was changed to .(dot).
# 
# Further, one gene names contained brackets (), which were replaced by dashes .

# 2 gene names are not matched:
missingGenes = insitu.genes_orig[which(!insitu.genes_orig %in% normalized.data.genes)]
# print(missingGenes)

## [1] "Blimp.1"      "E.spl.m5.HLH"

# this was reported by Nikos
# lets fix this by changing the . characters in the gene names to -
insitu.genes = gsub(".","-",insitu.genes_orig,fixed = T)
# also replace .spl. --> (spl)
insitu.genes = gsub("-spl-","(spl)",insitu.genes,fixed = T)

# assert that all institu genes appear in the gene names
stopifnot(all(insitu.genes %in% raw.data.genes))

# Now we can rename the genes in the institu.matrix with the correct names:

insitu.matrix = as.matrix(insitu.matrix)
colnames(insitu.matrix) = insitu.genes

# Read geometry data

# The column naming of the geometry is not consistent with distMap expectation, xcoord, ycoord and zcoord must be renamed to x, y, z

geometry = read.csv("Data/geometry.txt.gz",sep = " ")

## Warning in read.table(file = file, header = header, sep = sep, quote =
## quote, : seek on a gzfile connection returned an internal error

## Warning in read.table(file = file, header = header, sep = sep, quote =
## quote, : seek on a gzfile connection returned an internal error

colnames(geometry) = c("x","y","z")


####################################################
# Run DistMap for the several features identified by 
# - LASSO
# - Neural Nets
# - Baseline
####################################################

cat("Reading important features from all methods...\n")

# Load important features 
# - Lasso - 
# 20, 40,60 genes
files_cv_results_LASSO <- list.files(path = "LASSO_topX_workflow/Results_UniquelyMapped_cells_inSituRNAseq/", 
                                     pattern = "CV", 
                                     full.names = TRUE)

# - Neural Networks -
# 20, 40, 60 genes
files_cv_results_NeuralNets <- 
    list.files(path = "NeuralNetworks/PostChallenge/", 
               full.names = TRUE, pattern = "Results.NeuralNetworks.PostChallenge10CV.*FeatureSelection.InsituGenes.txt", recursive = TRUE)

# - Random - 
# 20, 40, 60 genes
files_cv_results_Random <- list.files(path = "Results_Common/Baseline_Method", 
                                      full.names = TRUE)

# Configuration for comparison
featuresInfo <- data.table(numFeatures = rep(x = c(20, 40, 60), each = 10, times = 3), 
                           Method = rep(c("Lasso", "NeuralNets", "Random"), each = 30), 
                           CV = c(rep(x = 1:10, times = 9)))

# Read provided training index - 10 fold CV
folds_train <- fread(input = "Data/CV_folds/folds_train.csv")
folds_test <- fread(input = "Data/CV_folds/folds_test.csv")
cell_ids <- fread(input = "Data/CV_folds/cell_ids.csv")

# Read Binarized data
load(file = "Results_Common/my_mapCells_run.RData")
binarizedExpression <- as.data.frame(dm@binarized.data)

# Run DistMap
cat("Running DistMap...\n")
registerDoMC(cores = 4)
distMaps <- foreach(i = 1:nrow(featuresInfo)) %dopar% {
    
    cat(i, "out of", nrow(featuresInfo), "\n")
    
    featuresInfo_i <- featuresInfo[i, ]
    
    if(featuresInfo_i$Method == "Random"){
        # Read file
        importantGenes <- read.table(file = grep(pattern = paste0(featuresInfo_i$numFeatures, "genes",
                                                                  "_CV_", featuresInfo_i$CV, ".txt"), 
                                                 x = files_cv_results_Random, value = TRUE), 
                                     as.is = TRUE)$V1
        
        
    }else if(featuresInfo_i$Method == "Lasso"){
        # Read File
        load(file = files_cv_results_LASSO[grepl(pattern = paste0("CV_", featuresInfo_i$CV, "_"),
                                                 x = files_cv_results_LASSO)])
        
        # Select features using lasso results
        importantGenes <- cv.lasso.mse$ImportantFeatures[
            numFeatures == featuresInfo_i$numFeatures, topGenes]
        
        # Remove
        rm(cv.lasso.mse, singleCellRNAseq.test, singleCellRNAseq.train, 
           geometry.RNAseq.train, geometry.RNAseq.test)
        gc()
        
    }else if(featuresInfo_i$Method == "NeuralNets"){

	 	  subchallenge = "error";
	     if (featuresInfo_i$numFeatures == 20) {
		     subchallenge = "Subchallenge3";
		  } else if (featuresInfo_i$numFeatures == 40) {
		     subchallenge = "Subchallenge2";
		  } else if (featuresInfo_i$numFeatures == 60) {
		     subchallenge = "Subchallenge1";
		  }

        # Read features using NeuralNets method
        importantGenes <- fread(files_cv_results_NeuralNets[
            grepl(pattern = paste0("Slice_",(featuresInfo_i$CV - 1), ".", subchallenge, ".FeatureSelection.InsituGenes.txt"), 
                  x = files_cv_results_NeuralNets)], 
            header = FALSE)$V1
        
    }else if (featuresInfo_i$Method == "all"){
        # use all insitu genes
        importantGenes <- colnames(insitu.matrix)
    }
    
    # Cells to test using provided folds
    cells_test_fold_i <- cell_ids[cell_ids$id %in% na.omit(as.numeric(folds_test[featuresInfo_i$CV, ])), ]
    
    # Normalized data test
    normalized.data.test <- normalized.data[, cells_test_fold_i$cell]
    
    # Raw data test
    raw.data.test <- raw.data[, cells_test_fold_i$cell]
    
    # Checks
    stopifnot(all(colnames(normalized.data.test) == colnames(raw.data.test)))
    stopifnot(all(rownames(normalized.data.test) == rownames(raw.data.test)))
    
    # FILTERING for the important features first
    insitu.matrix.2use <- insitu.matrix[, colnames(insitu.matrix) %in% importantGenes]
    
    # Evaluate if all selected genes are in insitu genes
    stopifnot(identical(x = which(!importantGenes %in% colnames(insitu.matrix)), 
                        y = integer(0)))
    
    # Run DistMap
    distMap_i = new("DistMap",
                    raw.data = raw.data.test,
                    data = normalized.data.test,
                    insitu.matrix = insitu.matrix.2use,
                    geometry = as.matrix(geometry))
    
    # # Evaluate quantiles effect
    distMap_i <- binarizeSingleCellData(object = distMap_i, quantiles = seq(from = 0.15, to = 0.5,
                                                                            by = 0.01))

    # Use information from Provided binarized expression table
    binarizedExpression.part <- binarizedExpression[match(x = colnames(insitu.matrix.2use), 
                                                          table = rownames(binarizedExpression)), 
                                                    match(x = colnames(normalized.data.test), 
                                                          table = colnames(binarizedExpression))]
    
    # Check - TO REMOVE
    stopifnot(all(rownames(binarizedExpression.part) == rownames(distMap_i@binarized.data)) & 
                  all(colnames(binarizedExpression.part) == colnames(distMap_i@binarized.data)))
    
    # Add to DistMap Object
    distMap_i@binarized.data <- as.matrix(binarizedExpression.part)
    
    # Run distMap
    distMap_i <- mapCells(distMap_i)
    
    # Check the input data dimensions
    stopifnot(ncol(distMap_i@insitu.matrix) == featuresInfo_i$numFeatures)
    
    return(distMap_i)
}
names(distMaps) <- paste(featuresInfo$Method, featuresInfo$numFeatures, "CV", featuresInfo$CV, 
                         sep = ".")

############################################################
# Top 10 locations predictions
# For each distMap object generate the top 10 predictions by extracting the “mcc.scores” object from within the distMap object and identifying the bins corresponding to the top 10 highest mcc scores per cell.
############################################################

cat("Predicting top 10 locations...\n")

topX <- 10
topXbins.All <- foreach(i = 1:length(distMaps), .combine = rbind) %dopar% {
    distMap_i <- distMaps[[i]]
    mccScores_i <- distMap_i@mcc.scores
    colnames(mccScores_i) <- colnames(distMap_i@data)
    
    # Top X bins per cell
    topXbins <- foreach(col_i = 1:ncol(mccScores_i), .combine = rbind) %do% {
        bins <- sort.int(x = mccScores_i[,col_i], decreasing = TRUE, 
                         index.return = TRUE, na.last = TRUE)$ix[1:topX]
        topXbins_i <- geometry[bins, ]
        topXbins_i$Cell <- colnames(mccScores_i)[col_i]
        topXbins_i$bin <- bins
        topXbins_i$id <- which(colnames(normalized.data) == colnames(mccScores_i)[col_i])
        return(topXbins_i)
    }
    
    # Add number of features information
    topXbins$numFeatures <- ncol(distMap_i@insitu.matrix)
    topXbins$type <- unlist(strsplit(x = names(distMaps)[i], split = ".", fixed = TRUE))[1]
    topXbins$CV <- unlist(strsplit(x = names(distMaps)[i], split = ".", fixed = TRUE))[4]
    
    # return
    return(topXbins)
}
topXbins.All <- data.table(topXbins.All[, c("id","Cell", "bin", "x", "y", "z", 
                                            "numFeatures", "type", "CV")])
topXbins.All[, type := paste(type, "DistMap", sep = ".")]

# Check that cells in CVs are as expected
# - No cells in training set are in test
# - All cells in test sets are used in predictions
for(fold_i in 1:nrow(folds_test)){
    stopifnot(all(unique(topXbins.All[CV == fold_i, id]) %in% na.omit(as.numeric(folds_test[fold_i, ]))))
    stopifnot(!any(unique(topXbins.All[CV == fold_i, id]) %in% na.omit(as.numeric(folds_train[fold_i, ]))))
}

# Save predictions
cat("Saving results...\n")
save(topXbins.All, 
     file = "Results_Common/cell_Locations_top10_MCC_DistMapOnTest_10foldCV_UsingProvidedBinaryData.RData")














