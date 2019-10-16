# Load Dependencies
library(glmnet)
library(doMC)
library(ggplot2)
library(ggfortify)

########################################
# LASSO.topX function Documentation ##
########################################
# LASSO.topX() performs a repeated cross validation scheme in order to identify a number of desired important features. 
# The function builds (nfols * reps * lambda) LASSO models.

# For each model, LASSO.topX() extracts the following information
# - Its performance. The loss function that is used to calculate the performance it employs the Euclidean of the predicted xyz location to the real location.  
# - The number of features which were used.
# - Their corrensponding coefficients

# Then for same size models - models with the same number of features - but with different lambda values the function keeps the models with the smaller lambda. As evident from figure2, models with smaller lambdas produce smaller meanEuclidean distance error. 

# After, for same size, same lambdas models we extact the features and for each one we calculate two metrics: 
# - Stability: the number of times a feature was selected as important  
# - Mean Coefficient: the mean value of the coefficients that a feature was assigned.

# Finally LASSO.topX() function utilizes the rankSum statistic to combine these two metrics. In cases where there are more than desired number of features in the final list the function returns the genes with the higher rankSum. 

# - Arguments -
# x : matrix, rows correspond to cells and columns to genes.
# y : matrix, rows correspond to cells and columns to coordinates
# lambda : numeric, the lambdas for which to fit the LASSO model.
# nfolds : integer, number of fold to be created during cross validation
# numFeatures : integer, vector of integers, corresponding to the number of features we are interested in identifying. 
# reps : integer, corresponding to the number of times to repeat the cross validation procedure.
########################################
LASSO.topX <- function(x, y, lambda, nfolds = 10, numFeatures = c(20,40,60), reps, ncores = 32){
    
    # Euclidean Distance Matrix
    euclidean_distance <- function(myPred_i, test_i){
        out <- apply(X = (myPred_i - test_i), MARGIN = 1, 
                     FUN = function(li){sqrt(sum(li^2))})
        return(out)
    }
    
    # parallel backend
    registerDoMC(cores = ncores)
    
    # Build LASSO models
    models.cvErrors <- foreach(rep_i = 1:reps) %dopar%{
        
        # Fold indexes
        foldid <- sample(x = c(1:nfolds), size = nrow(x), replace = TRUE)
        
        # train cv models
        models.cv <- foreach(fold_i = 1:nfolds) %do% {
            
            # train
            fit.glmnet <- glmnet(x = x[foldid != fold_i, ], y = y[foldid != fold_i, ], 
                                 alpha = 1, lambda = lambda, family = "mgaussian")
            
            return(fit.glmnet)
        }
        
        # Calculate cv errors
        cv.errors <- foreach(fold_i = 1:nfolds) %do% {
            
            # Get Model
            model_i = models.cv[[fold_i]]
            
            # test
            myPred <- predict(object = model_i, newx = x[foldid == fold_i,], 
                              s = model_i$lambda)
            
            # cv errors
            cvErrors.tmp <- apply(X = myPred, MARGIN = 3, FUN = euclidean_distance, 
                                  test_i = y[foldid == fold_i, ])
            cvErrors.tmp <- data.table(cvErrors.tmp)
            
            # Add fold Information
            cvErrors.tmp$fold <- fold_i
            return(cvErrors.tmp)
            
        }
        cv.errors <- rbindlist(l = cv.errors, fill = TRUE)
        
        # Add Rep Information
        cv.errors$Rep <- rep_i
        
        # return
        return(list(cv.errors = cv.errors, models.cv = models.cv))
        
    }
    
    # Extract Cross Validation Errors
    cv.errors.out <- rbindlist(l = lapply(X = models.cvErrors, FUN = function(li){
        li$cv.errors
    }), fill = TRUE)
    cv.errors.out[, RepFold := paste("Rep", Rep, "Fold", fold, sep = ".")]
    
    # Extract Models
    models.Out <- lapply(X = models.cvErrors, FUN = function(li){
        models.cvErrors[[1]]$models.cv
    })
    models.Out <- unlist(models.Out, recursive = FALSE)
    names(models.Out) <- unique(cv.errors.out$RepFold)
    
    # Melt data
    cv.errors.m <- melt.data.table(cv.errors.out[, -c("fold", "Rep")], id.vars = c("RepFold"), 
                                   variable.name = "lambdaid", value.name = "EuclideanDistance")
    
    # Calculate Mean Euclidean Distance Error
    cv.errors.m[, MeanEuclideanDistance := mean(EuclideanDistance), by = "lambdaid"]
    
    # Add Lambda information
    cv.errors.m[, lambda := lambda[lambdaid]]
    
    
    # Extract important features #
    
    # Extract Lasso coefficients - Stability
    lasso.coefs.cv <- foreach(RepFold_i = unique(cv.errors.m$RepFold)) %dopar%{
        
        # Select Model
        model_i <- models.Out[[which(names(models.Out) == RepFold_i)]]
        
        # Extract Coefficients across all coordinates for each Lambda
        coefs = coef(model_i, s = model_i$lambda)
        
        coefs.tmp <- foreach(coord_i = 1:length(coefs), .combine = rbind) %do%{
            coef_out <- data.table(as.matrix(coefs[[coord_i]])[-1, ], 
                                   keep.rownames = TRUE)
            colnames(coef_out)[1] <- c("Genes")
            coef_out$Coord <- names(coefs)[coord_i]
            return(coef_out)
        }
        coefs.tmp$RepFold <- RepFold_i
        return(coefs.tmp)
    }
    lasso.coefs.cv.dt <- rbindlist(l = lasso.coefs.cv, use.names = TRUE)
    lasso.coefs.cv.dt.m <- melt.data.table(data = lasso.coefs.cv.dt, 
                                           id.vars = c("Genes", "Coord", "RepFold"), 
                                           variable.name = "lambdaid", 
                                           value.name = "Coefficiens")
    lasso.coefs.cv.dt.m[, lambda := lambda[lambdaid]]
    
    # Keep important genes for all coordinates
    lasso.coefs.cv.dt.m[, numCoordsPerLambdaPerGenePerFold := sum(Coefficiens != 0), 
                        by = c("lambdaid", "Genes", "RepFold")]
    # Filtering
    lasso.coefs.cv.dt.m <- lasso.coefs.cv.dt.m[numCoordsPerLambdaPerGenePerFold != 0]
    
    # Number of Features Per Lambda Per RepFold Per Coordinate
    # This corresponds also to the
    # MINIMUM Number of Features Per Lambda Per RepFold ACROSS ALL COORDINATES 
    # due to the previous filtering
    lasso.coefs.cv.dt.m[, numFeaturesPerLambdaPerFoldPerCoord := sum(Coefficiens != 0), 
                        by = c("lambdaid", "RepFold", "Coord")]
    
    # Check if the desired number of features in results
    if(!any(numFeatures %in% unique(lasso.coefs.cv.dt.m$numFeaturesPerLambdaPerFoldPerCoord))){
        cat("Stopping..\n-Information-\n")
        cat("Desired number of features:", numFeatures, "\n")
        cat("Resulting number of features:\n", 
            sort(unique(lasso.coefs.cv.dt.m$numFeaturesPerLambdaPerFoldPerCoord)), "\n")
        stop("The desired number of features not in results\nSuggestions: pick a different lambda range")
    }
    if(!all(numFeatures %in% unique(lasso.coefs.cv.dt.m$numFeaturesPerLambdaPerFoldPerCoord))){
        cat("-Information-\n")
        cat("Desired number of features:", numFeatures, "\n")
        cat("Resulting number of features:\n", 
            unique(lasso.coefs.cv.dt.m$numFeaturesPerLambdaPerFoldPerCoord), "\n")
        warning("Some of the desired number of features not in results\nSuggestions: pick a different lambda range\nThe computations will continue.")
    }

    # Number of folds each Gene is important Per lambda Per Coord
    lasso.coefs.cv.dt.m[, numFoldsPerLambdaPerGenePerCoord := sum(Coefficiens != 0), 
                        by = c("lambdaid", "Genes", "Coord")]
    
    # Select best features
    topGenesAll <- foreach(numFeatures_i = numFeatures, .combine = rbind) %do% {
        
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
        filterGenes[, Genes := factor(x = Genes, levels = Genes[
            sort.int(RankSum, index.return = TRUE)$ix])]

        topGenes <- as.character(filterGenes$Genes)[
            sort.int(filterGenes$RankSum, decreasing = TRUE, 
                     index.return = TRUE)$ix][1:numFeatures_i]
        
        return(data.table(topGenes, numFeatures = numFeatures_i))
        
    }
    
    # Generate output
    out <- list(CrossValidationErrors = cv.errors.m, 
                LassoCoeffient = lasso.coefs.cv.dt.m,
                ImportantFeatures = topGenesAll)
    
    return(out)
}




