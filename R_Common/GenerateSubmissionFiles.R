# Generate Submissions files
rm(list = ls())

# ***************
# USER input ####
# ***************

# # @ Test Cells - provided Binarized Table ####
# load(file = "Results_Common/cell_Locations_top10_MCC_DistMapOnTest_10foldCV_UsingProvidedBinaryData.RData")
# path2save <- "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/"

# # @ TrainTest mode of DistMap ####
# load(file = "Results_Common/cell_Locations_top10_MCC_DistMapTrainTest_10foldCV.RData")
# path2save <- "Results_Common/SubmissionFiles_CV_DistMapTrainTest/"


# *****************
# CALCULATIONS ####
# *****************
require(data.table)
require(foreach)

# Function to write submission files
write.submittion.file <- function(mypredictions, features, inSituGenes, filename){
    
    # Checks # 
    # All genes are in the insitu data
    if(!all(features %in% inSituGenes)){
        stop("Not all features in insitu genes")
    }
    
    mypredictions$tmp <- rep(x = paste("Prediction", c(letters[c(1:10)]), sep = "."), 
                             times = (nrow(mypredictions) / 10))
    out2csv.predictions <- dcast.data.table(data = mypredictions, formula = id ~ tmp, 
                                            value.var = "bin")
    
    if(any(is.na(out2csv.predictions))){
        stop("NAs in predictions, recalculate predictions!")
    }
    colnames(out2csv.predictions)[-1] <- paste0("V", c(1:10))
    out2csv.features <- data.table(matrix(data = features, nrow = length(features) / 10, 
                                          ncol = 10))
    out2csv <- rbindlist(list(out2csv.features, out2csv.predictions), fill = TRUE)
    
    # Count NAs
    if(sum(is.na(out2csv)) > (length(features) / 10) ){
        stop("More NAs than expented!!")
    }
    
    # Reorder columns
    out2csv <- out2csv[, c("id", colnames(out2csv)[1:10]), with = FALSE]
    
    # Write file
    fwrite(x = out2csv, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE, 
           na = "")
    
}

# Keep 20 / 40 / 60 files
topXbins.All[, tmpName := paste(type, numFeatures, "CV", CV, sep = ".")]
topXbins.All <- topXbins.All[numFeatures %in% c(20, 40, 60)]

# Important features, files information
# - Lasso - 
files_cv_results_LASSO <- 
    list.files(path = "Modified_LASSO_workflow/Results_UniquelyMapped_cells_inSituRNAseq", 
               pattern = "FeaturesInformation_CV", 
               full.names = TRUE)
# - Neural Nets - 
files_cv_results_NeuralNets <- 
    list.files(path = "NeuralNetworks/Results_UniquelyMapped_cells_inSituRNAseq",
               recursive = TRUE, full.names = TRUE)

# - Random - 
# 20, 40, 60 genes
files_cv_results_Random <- list.files(path = "Results_Common/Baseline_Method", 
                                      full.names = TRUE)


# Insitu genes information
insitu <- fread("Data/bdtnp.txt.gz")

# Subset on test folds
# Load Test sets
folds_test <- fread(input = "Data/CV_folds/folds_test.csv")
topXbins.All2submit <- foreach(foldi = 1:nrow(folds_test), .combine = rbind) %do% {
    topXbins_Foldi <- topXbins.All[CV == foldi & id %in% as.numeric(folds_test[foldi, ]), ]
}

counter <- 0
for(tmpName_i in unique(topXbins.All2submit$tmpName)){
    counter <- counter + 1    
    
    topXbins.All_i <- topXbins.All2submit[tmpName_i == tmpName]
    
    # Message
    cat("Writing predictions, ", counter, "out of ", 
        length(unique(topXbins.All2submit$tmpName)), "\n")
    
    # load features information
    if(grepl(pattern = "Lasso", x = unique(topXbins.All_i$type))){
        
        # Read features using Lasso method
        load(file = files_cv_results_LASSO[grepl(pattern = paste0("CV_", unique(topXbins.All_i$CV),"_"),
                                                 x = files_cv_results_LASSO)])    
        importantGenes <- cv.lasso.mse$ImportantFeatures[numFeatures == unique(topXbins.All_i$numFeatures), 
                                                         topGenes]
        
        rm(cv.lasso.mse)
    }else if(grepl(pattern = "NeuralNets", x = unique(topXbins.All_i$type))){
        
        # Read features using NeuralNets method
        importantGenes <- fread(files_cv_results_NeuralNets[
            grepl(pattern = paste0(strsplit(x = tmpName_i, split = ".", fixed = TRUE)[[1]][3], 
                                   "_WithCut_",
                                   (as.numeric(strsplit(x = tmpName_i, split = ".", 
                                                        fixed = TRUE)[[1]][5]) - 1), 
                                   ".txt"), 
                  x = files_cv_results_NeuralNets)], 
            header = FALSE)$V1
        
        
    }else if(grepl(pattern = "Random", x = unique(topXbins.All_i$type))){
        
        # Read file
        importantGenes <- read.table(file = grep(pattern = paste0(strsplit(x = tmpName_i, split = ".", 
                                                                           fixed = TRUE)[[1]][3], 
                                                                  "genes_CV_",
                                                                  strsplit(x = tmpName_i, split = ".", 
                                                                           fixed = TRUE)[[1]][5], 
                                                                  ".txt"),
                                                 x = files_cv_results_Random, value = TRUE), 
                                     as.is = TRUE)$V1
        
    }

    # Write file
    write.submittion.file(mypredictions = topXbins.All_i, 
                          features = importantGenes,
                          inSituGenes = colnames(insitu), 
                          filename = paste0(path2save, tmpName_i, ".csv"))

}






