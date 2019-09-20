##########################################
# SUMMARY
# The script evaluates different methods predictions,
# - LASSO
# - Neural Nets
# - Baseline
# It calculates the mean Euclidean distance of the top 10 predictions per cell from its "true" position and produces the results plots
#########################################

# Load dependencies ##
rm(list = ls())
library(data.table)
library(ggplot2)
library(ggfortify)

# ***************
# USER INPUT ####
# ***************

folder2save <- "Results_Common/Figures/"


# *****************
# CALCULATIONS ####
# *****************

# Euclidean Distance two points
euclidean_distance <- function(myPred_i, true_position_i){
    return(sqrt(sum((myPred_i - true_position_i)^2)))
}


# Load necessary data ##

# Load Cell locations ####

# @ Using provided Binarized Data ####
load("Results_Common/cell_Locations_top10_MCC_DistMapOnTest_10foldCV_UsingProvidedBinaryData.RData")
topXbins.All_ProvidedBinaryData <- topXbins.All
topXbins.All_ProvidedBinaryData[, type := paste(type, "testCells_PB", sep = "_")]
topXbins.All_ProvidedBinaryData[, type := sub(pattern = ".DistMap", replacement = "", x = type)]
topXbins.All_ProvidedBinaryData[, Cell_test_CV := paste0(Cell, "_CV", CV)]

# @ Binarize data - train cell, selected genes ####
load("Results_Common/cell_Locations_top10_MCC_DistMapTrainTest_10foldCV.RData")
topXbins.All_trainCells_selGenes <- topXbins.All
topXbins.All_trainCells_selGenes[, Cell_test_CV := paste0(Cell, "_CV", CV)]
topXbins.All_trainCells_selGenes[ 
    , CV_test := Cell_test_CV %in% topXbins.All_ProvidedBinaryData$Cell_test_CV]
topXbins.All_trainCells_selGenes[CV_test == FALSE, type := paste(type, "trainCells_selGenes", sep = "_")]
topXbins.All_trainCells_selGenes[CV_test == TRUE, type := paste(type, "testCells_selGenes", sep = "_")]
topXbins.All_trainCells_selGenes[, type := sub(pattern = ".DistMap", replacement = "", x = type)]

rm(topXbins.All)


# Combine everything
topXbins.All <- rbindlist(l = list(topXbins.All_ProvidedBinaryData, 
                                   topXbins.All_trainCells_selGenes[, -"CV_test"]))
topXbins.All[, Cell_test_CV := NULL]

colnames(topXbins.All)[8] <- "Method"

# Get original RNAseq bins x, y, z 
load("Results_Common/geometry.RNAseq.RData")

# Add true positions to predictions
# Note: here we add the predictions of all cells and not only those that are 
# uniquely mapped
topXbins.All <- merge(x = topXbins.All, y = geometry.RNAseq.all, by = "Cell")
topXbins.All[, unique_line_ID := paste("line", 1:nrow(topXbins.All), sep = ".")]

# Calculate Euclidean Distances Per Prediction
topXbins.All[, EuclideanDistancePerPrediction := 
                 euclidean_distance(myPred_i = c(x, y, z), 
                                    true_position_i = c(xcoord, 
                                                        ycoord, 
                                                        zcoord)), 
             by = "unique_line_ID"]

# Calculate Mean Euclidean Distances 
# - Per Cell 
# - Per number of Features 
# - Per Method
topXbins.All[, MeanEuclideanPerCellPerFold := mean(EuclideanDistancePerPrediction), 
             by = c("Cell", "numFeatures", "Method", "CV")]


topXbins.All <- unique(topXbins.All[, c("Cell", "numFeatures", "Method", "CV", 
                                        "MeanEuclideanPerCellPerFold")])


# Calculate Mean Euclidean Distance across all cells
# - Per number of Features 
# - Per method
topXbins.All[, MeanEuclideanAllCellPerFold := mean(MeanEuclideanPerCellPerFold), 
             by = c("numFeatures", "Method", "CV")]


# Set methods order for plots
topXbins.All$Method <- factor(x = topXbins.All$Method, 
                              levels = c("Lasso_testCells_PB",
                                         "Lasso_trainCells_selGenes", 
                                         "Lasso_testCells_selGenes",
                                         "NeuralNets_testCells_PB",
                                         "NeuralNets_trainCells_selGenes",
                                         "NeuralNets_testCells_selGenes", 
                                         "Random_testCells_PB",
                                         "Random_trainCells_selGenes",    
                                         "Random_testCells_selGenes"))


# Compare all methods
myplots <- list()
myplots[[1]] <- 
    ggplot(data = topXbins.All, 
       mapping = aes(x = as.factor(numFeatures), y = MeanEuclideanPerCellPerFold,
                     fill = Method)) +
    geom_boxplot(notch = TRUE) +
    ggtitle("Methods Comparison") +
    ylab(label = "Mean Euclidean Distance (per cell)") +
    xlab(label = "Number of Features")

# Evalute mean values #
data2plot <- unique(topXbins.All[, c("numFeatures", "Method", "MeanEuclideanAllCellPerFold")])
# Mean performance across nested cv
data2plot[, MeanPerformance := mean(MeanEuclideanAllCellPerFold), by = c("numFeatures", "Method")]

myplots[[2]] <- 
    qplot(x = as.factor(numFeatures), y = MeanEuclideanAllCellPerFold, 
          data = data2plot, #[!grepl(pattern = "trainCells", x = data2plot$Method)], 
          fill = Method, geom = "blank") + 
    geom_boxplot() +
    ylab(label = "MeanEuclDistPerFold, (Nested CV)") + 
    xlab(label = "Number of Features") +
    ggtitle("Methods Comparison") + 
    theme_bw()

myplots[[3]] <- 
    qplot(x = as.factor(numFeatures), y = MeanPerformance, 
      data = data2plot, #[!grepl(pattern = "trainCells", x = data2plot$Method), ], 
      fill = Method,
      geom = "blank") + 
    xlab(label = "Number of Features") + 
    ylab(label = "MeanEuclDistAllFold, (Nested CV)") +
    geom_bar(position = position_dodge(), stat = "identity") + 
    ggtitle("Methods Comparison") + 
    theme_bw()

print(autoplot(myplots[c(2,3)]))

ggsave(filename = paste0(folder2save, "MethodsComparison_NestedCV_Distribution.jpg"), 
       plot = myplots[[2]], 
       width = 6, height = 4)

ggsave(filename = paste0(folder2save, "MethodsComparison_NestedCV_mean.jpg"), 
       plot = myplots[[3]], 
       width = 6, height = 4)




