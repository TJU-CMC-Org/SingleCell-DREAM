Authors
-------
- Phillipe Loher (http://cm.jefferson.edu/phillipe-loher/) 
- Nestoras Karathanasis (http://cm.jefferson.edu/nestoras-karathanasis/)

We are part of the Computational Medicine Center at Thomas Jefferson University.  To find out more about our research, please visit https://cm.jefferson.edu. 


General
-------

This is a technical report on the techniques and protocols used to
perform the analysis on the Single Cell Transcriptomics challenge, see
<https://www.synapse.org/#!Synapse:syn15665609/wiki/>. Please follow the
steps below to rerun our analysis.

### Software Dependencies

* Unix System
* Python 3.5+, packages PyTorch, numpy, pandas, pickle, csv, sys, random.
* R version 3.6.0 (2019-04-26) -- "Planting of a Tree"

###### R packages


    # Install required packages for Modified LASSO workflow 
    # Start R
    # Set working directory to SingleCell-DREAM_FOLDER
    setwd("PATH_TO_SingleCell-DREAM_FOLDER/")
    # Installing required dependencies
    source("Modified_LASSO_workflow/R/Dependencies.R") 


### Setting things up

    # Setting Things up
    # - Create folders tree
    # - Download data from: 
    # https://shiny.mdc-berlin.de/DVEX/
    source("R_Common/SettingThingsUp.R")
    
    # Create NestedCV folds
    # The created 10 folds are generated randomly and are used only for illustration purposes.
    # The results in our paper are based on the 10 cross validation folds provided by the challenge's organiners and are available here.TBA
    source("R_Common/Create_NestedCV_folds.R")


### Feature Selection 
We employed three methods to perform the feature selection step, namely Random, a modified LASSO workflow and Deep Neural Nets.

#### Baseline - Random
We randomly selected the desired number of genes to baseline our feature selection algorithms.

    # Randomly select genes
    source(file = "Modified_LASSO_workflow/R/Baseline.R")


#### Modified LASSO workflow

##### Generating the 3d positions of the Cells

In order to generate the cell's 3d positions, which we used for labels, we run DistMap using the code provided online. Then we extracted the “mcc.scores” object from within distMap’s output and the bins corresponding to the maximum mcc.score per cell were identified. There are 287 cells that are not uniquely mapped to only one bin. We used only the uniquely mapped cells, 1100 out of 1297, in our feature selection process.

    # Run Distmap to identify 3d cell locations
    source(file = "R_Common/ReproducePaperResults.R")
    
    # Identify Cell positions
    source(file = "Modified_LASSO_workflow/R/ExtractRNAseqLocations.R")


##### Feature selection step
We modified LASSO workflow as described in our publication, TBA.
Our code and respective documentation can be found in 
`Modified_LASSO_workflow/R/glmnetExtensionLibrary.R`
and is provided as an extension of glmnet package.

    # Reproduce our gene selection using the modified LASSO workflow using 
    # the nested cross validation folds.
    # - To select only from the inSitu genes set: useOnlyInSitu <- TRUE. 
    # - To select across all genes set: useOnlyInSitu <- FALSE. 
    # Run feature selection process
    source(file = "Modified_LASSO_workflow/R/LASSO_features_CVs_train.R")
    
    # Plot error evolution and selected features of one of the provided Nested cross validation folds. 
    # Variable `fileCV` in the USER INPUT section of the script can be used to specify the results file that you want to use. 
    source(file = "Modified_LASSO_workflow/R/Plot_Lasso_Features.R")


#### Instructions for running Neural Network codes used for variable selection
Please visit the sub-directory named "NeuralNetworks/DuringChallenge_Subchallenge2/" and then run the below 3 steps.

#### Pre-process/generate a table that will be used for training:
* Run: Rscript Step1_GenTrainingDataMatrix.R

#### Create the models and perform inference:
* Run: python EvalData.py

#### Combine variable importance scores across all 200 models:
* Run: bash Step3_GenerateRankedLists.sh

#### Run the PostChallenge 10CV
* Go to directory: NeuralNetworks/PostChallenge/
* Run: python EvalData.py
* Run: bash Step3_GenerateRankedLists.sh

### Location prediction
After selecting the most informative genes, using Random, the modified version of LASSO and Deep Neural Nets, we predicted the 10 locations per cell using a modified version of DistMap, as described in our publication, add link TBA. The modified version of DistMap employs only the cells in the training set to calculate all DistMap parameters and predicts the cell locations in the both the training and test sets. 
The modified version of DistMap can be found here
`R_Common/distmap/R/myDistMap.R`


    # Predict Locations using the modified version of DistMap
    source(file = "R_Common/LocationPredictions_DistMap_TrainTest.R")
    
    # Predict Locations using the provided binarized table
    source(file = "R_Common/LocationPredictions_DistMapOnTest_ProvidedBinarizedData.R")


### Score predictions 
We score our prediction using our in-house blind metric and the challenge's organizers score functions, which are available online at `https://github.com/dream-sctc`.

For details see publication - TBA

    # Score predictions - blind metric
    source(file = "R_Common/Evaluate_Predictions_Blind_Metric.R")

    # Score predictions - organizers score functions
    # Open R script "R_Common/GenerateSubmissionFiles.R"
    # - Uncomment lines 9 and 10 to generate the submission files needed using the Provided Binarized table.
    source(file = "R_Common/GenerateSubmissionFiles.R")
    # - Uncomment lines 13 and 14 to generate the submission files needed using the location predictions from the modified DistMap.
    source(file = "R_Common/GenerateSubmissionFiles.R")
    # Finally generate scores and figures for both approaches by running.
    source(file = "R_Common/Evaluate_Predictions_Organizers_Scores.R")
