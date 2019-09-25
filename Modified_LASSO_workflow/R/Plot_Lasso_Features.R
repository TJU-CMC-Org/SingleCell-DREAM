# plot LASSO feature selection
rm(list = ls())

# ***************
# USER INPUT ####
# ***************

# File with CV results to load
fileCV <- "Modified_LASSO_workflow/Results_UniquelyMapped_cells_inSituRNAseq/FeaturesInformation_CV_1_Data_UniquellyMappedCells_LASSO_Reps20_Folds5.RData"

# Num of Features 
numFeatures <- c(20, 40, 60)

# Folder to save plots
folder2save <- "Modified_LASSO_workflow/Results_UniquelyMapped_cells_inSituRNAseq/"

# Folds
nfolds <- 5

# Reps
reps <- 20


# CALCULATIONS ####
library(data.table)
library(ggplot2)
library(ggfortify)

# Load information
load(fileCV)

# Extract Coefficient information
lasso.coefs.cv.dt.m <- cv.lasso.mse$LassoCoeffient

# Get Features information
lambdaFeatures <- unique(lasso.coefs.cv.dt.m[
    numFeaturesPerLambdaPerFoldPerCoord %in% numFeatures,
    c("numFeaturesPerLambdaPerFoldPerCoord", "lambda")])


# @ Plot performance evolutions ####
cv.errors.m <- cv.lasso.mse$CrossValidationErrors
cv.errors.m <- merge(x = cv.errors.m, 
                     y = unique(lasso.coefs.cv.dt.m[,c("RepFold", 
                                                   "lambda",
                                                   "numFeaturesPerLambdaPerFoldPerCoord")]), 
                     by = c("RepFold", "lambda"))

data2plot.tmp <- unique(cv.errors.m[, c("lambda", "MeanEuclideanDistance", 
                                        "numFeaturesPerLambdaPerFoldPerCoord")])

data2plot.tmp[, Number_Of_Features := ifelse(test = numFeaturesPerLambdaPerFoldPerCoord %in% numFeatures,
                                   yes = numFeaturesPerLambdaPerFoldPerCoord, 
                                   no = "No")]

if(all(unique(data2plot.tmp$Number_Of_Features) == "No")){
    errorMessage <- paste("No models with the desired number of features.
                          Available number of features are:",  
                          paste(sort(as.numeric(unique(data2plot.tmp$numFeaturesPerLambdaPerFoldPerCoord))), collapse = " "))
    stop(errorMessage)
}



plotLambdaSelection <- 
    ggplot(data = data2plot.tmp, mapping = aes(x = log(lambda), 
                                           y = MeanEuclideanDistance)) +  
    geom_line() +
    geom_point(mapping = aes(x = log(lambda), y = MeanEuclideanDistance,
                             colour = Number_Of_Features),
               data = data2plot.tmp[Number_Of_Features != "No"]) +
    ylab(label = "Mean Euclidean Distance") +
    ggtitle("Error Evolution") + 
    theme_bw() + theme(legend.position = "bottom")
print(plotLambdaSelection)

ggsave(filename = paste0(folder2save, "/ErrorEvolution.jpg"), 
       plot = plotLambdaSelection, width = 4.2, height = 4)


data2plot <- cv.errors.m[lambda %in% lambdaFeatures$lambda, c(1, 2, 4, 5)]
data2plot <- merge(x = data2plot, y = lambdaFeatures, by = "lambda")
data2plot[, log.lambda := as.factor(round(log(lambda), digits = 3))]
data2plot[, numFeaturesPerLambda := as.factor(numFeaturesPerLambdaPerFoldPerCoord)]


# @ Plots features selection process ####

# Select best features
for(i in 1:length(numFeatures)){
    numFeatures_i = numFeatures[i]
    
    # Filter on the number of features
    # - The Selected genes are predictive for ALL Coordinates
    # - Are important in at least ONE cross validation RepFold
    # Stability: Features which are present in more cross validation folds 
    # are more stable than features in only one cross validation RepFold
    filterFeatures <- lasso.coefs.cv.dt.m[ numFeaturesPerLambdaPerFoldPerCoord == numFeatures_i & 
                                               numFoldsPerLambdaPerGenePerCoord > 1, ]
    
    # If there are more than one lambda keep the features from the best lambda
    # Best lambda - Lambda with the minimum euclidean distance
    if(length(unique(filterFeatures$lambda)) > 1){
        
        # Add mean euclidean distance info
        filterFeatures <- merge(x = filterFeatures, 
                                y = unique(cv.errors.m[, c("MeanEuclideanDistance", "lambda")]), 
                                by = "lambda")
        # Select best lambda
        best_lambda <- filterFeatures$lambda[which.min(filterFeatures$MeanEuclideanDistance)]
        
        # Keep features only for best lambda
        filterFeatures <- filterFeatures[lambda == best_lambda, ]
    }

    # Sort genes based on STABILITY and absolute value of coefficients
    # Using RankSUM statistic
    filterFeatures[, meanCoef := mean(abs(Coefficiens)), by = c("Genes")]
    filterGenes <- unique(filterFeatures[, c("Genes", "numFoldsPerLambdaPerGenePerCoord", 
                                             "meanCoef")])
    # Calculate RankSum
    filterGenes[, rankCoef := frank(x = meanCoef, ties.method = "average")]
    filterGenes[, rankStability := frank(x = numFoldsPerLambdaPerGenePerCoord, 
                                         ties.method = "average")]
    filterGenes[, RankSum := (rankStability + rankCoef)]
    
    filterGenes[, numFolds := mean(numFoldsPerLambdaPerGenePerCoord), by = "Genes"]
    # filterGenes <- unique(filterGenes[,c("Genes", "numFolds")])
    filterGenes[, Genes := factor(x = Genes, levels = Genes[
        sort.int(RankSum, index.return = TRUE)$ix])]
    
    # Extract number of features in csv file
    fwrite(x = filterGenes, file = paste0(folder2save, numFeatures_i, "_Features.csv"))
    
    topGenes <- as.character(filterGenes$Genes)[
        sort.int(filterGenes$RankSum, decreasing = TRUE, 
                 index.return = TRUE)$ix][1:numFeatures_i]
    
    filterGenes[, SelectedGenes := Genes %in% topGenes]
    filterFeatures[, Genes.f := factor(Genes, levels = levels(filterGenes$Genes))]
    filterFeatures[, SelectedGenes := Genes %in% topGenes]
    
    # Initialize plots
    myplots <- list()
    myplots[[1]] <- 
        ggplot(data = filterGenes, 
               mapping = aes(x = Genes, y = numFolds, fill = SelectedGenes)) + 
        geom_bar(stat = "identity") + 
        # scale_fill_manual(values = c("grey80", "grey60")) +
        # scale_color_manual(values = c("black", "white")) +
        ylab(label = "Stability (Number of Folds)") + 
        scale_y_continuous(breaks = seq(from = 0, to = (nfolds * reps), by = 25),
                           limits = c(0, (nfolds * reps))) +
        coord_flip() + 
        ggtitle(paste("Selecting",numFeatures_i, "features")) + 
        theme_bw()
    
    
    myplots[[2]] <- 
        ggplot(data = filterFeatures, mapping = aes(x = Genes.f, y = abs(Coefficiens), 
                                                    fill = SelectedGenes)) + 
        geom_boxplot() + 
        geom_point(mapping = aes(x = Genes.f, y = meanCoef), colour = "black") +
        # scale_fill_manual(values = c("grey80", "grey60")) +
        coord_flip() + 
        ggtitle(paste("Selecting",numFeatures_i, "features")) +
        xlab(label = "") + theme_bw()

    # Print plots
    print(autoplot(myplots))
    
    ggsave(filename = paste0(folder2save, "Lasso_Stability_", numFeatures_i, "_features.jpg"), 
           plot = myplots[[1]], width = 5, height = 4 * i)
    
    ggsave(filename = paste0(folder2save, "Lasso_absCoeff_", numFeatures_i, "_features.jpg"), 
           plot = myplots[[2]], width = 5, height = 4 * i)
    
    
    
    
}



