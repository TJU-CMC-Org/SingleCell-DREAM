Authors
-------
- Phillipe Loher (http://cm.jefferson.edu/phillipe-loher/) 
- Nestoras Karathanasis (http://cm.jefferson.edu/nestoras-karathanasis/)

We are part of the Computational Medicine Center at Thomas Jefferson University.  To find out more about our research, please visit https://cm.jefferson.edu .  


General
-------

The source code accompanies our publication (PMID and link: TBA) describing
our approach to the DREAM Single Cell Transcriptomics challenge <https://www.synapse.org/#!Synapse:syn15665609/wiki/>. 

Notes

-   Parallel processing. There are several code chunks where we run
    computations in parallel using all available cores.


### Sub-Challenge 1 and 3

Before running the following steps you need to install all the
dependencies from the file `R/Dependencies`.

    # Setting things up to rerun analysis #
    
    # Results_TEST folder:
    # You will need the Results_TEST folder as there is where the algorithm will save its outputs. Use the Results_TEST folder provided so you can generate also the comparison figures of our writeUp. 
    
    # Set working directory to src folder
    setwd("~/Documents/Projects/DREAM_Challenge_Upload/src/")
    
    # Installing required dependencies
    source(file = "SubChallenge_1_3/R/Dependencies.R")


#### Generating the true positions of the Cells

In order to identify the “true” cell positions we run distMap using the
code provided online. Then we extracted the “mcc.scores” object from
within distMap’s output and the bins corresponding to the maximum
mcc.score per cell were identified. In contrast to the webinar and as
shown in the generated Figure there are 287 cells that are not uniquely mapped to
only one bin. We used only the uniquely mapped cells, 1100 out of 1297,
in our feature selection process

    # Run DistMap
    source(file = "SubChallenge_1_3/R/ReproducePaperResults.R")
    
    # Identify Cell positions
    source("SubChallenge_1_3/R/ExtractRNAseqLocations_SUBMIT.R")


#### LASSO

We trained LASSO as it is described in our write up report.

#### Method Generalization

In order to evaluate the generalization of our method on an unseen
dataset, we split the full dataset to

-   training set (770 cells), 70% of the cells which are mapped
    uniquely, and
-   hold out set (330 cells), 30% of the cells which are mapped
    uniquely.

We run LASSO using only the training set to select the important
features and we predict the location of the cells to both training and
hold out set using distMap. As evident from the last figure of this
section, our method is able to generalize to the hold out set and to
avoid overfitting as the performance of the hold out set and the
training sets are very similar. Intermediate figures show,

-   euclidean distance error in relation to the values of log lambda,
-   stability, the number of times a feature was selected as important across the repeated cross validation procedure and
-   the distribution of the coefficients per selected feature across the
    repetitive cross validation process.

How to run:

    # Extract Features using LASSO employing 70% of the cells which are mapped uniquely
    # Note: You should change the paths of the input files to rerun the script
    source("SubChallenge_1_3/R/LASSO_features_07train_SUBMIT.R")


#### Final Run

Similar to above we run lasso using the 1100 uniquely mapped cells to
identify the important features and we fed these features to distMap in
order to predict the top 10 locations of the 1297 cells.

    # Extract Features using LASSO employing 100% of the cells which are mapped uniquely (1100 cells). 
    # Note: You should change the paths of the input files to rerun the script
    source("SubChallenge_1_3/R/LASSO_features_all_data_SUBMIT.R")
    
    # Employ selected features and distMap to predict the final locations of the cells
    # using the several methods
    # You should change the paths of the input files to rerun the script
    source("SubChallenge_1_3/R/FinalLocationPredictions_SUBMIT.R")


### Sub-Challenge 2

#### Instructions for running Neural Network codes used for variable selection

#### Pre-process/generate a table that will be used for training:
* Other: Depending on the paths on your system, you may have to change the paths to the input files for files: dge_raw.txt, dge_normalized.txt, binarized_bdtnp.csv, and geometry.txt (these are provided by the challenge organizers at https://www.synapse.org/#!Synapse:syn15665609/wiki/582909)

* Run: Rscript Step1_GenTrainingDataMatrix.R

#### Create the models and perform inference:
* Pre-reqs include: Python 3.5+, packages PyTorch, numpy, pandas, pickle, csv, sys, random
* Other: at top of EvalData.py file, you will see variables that start with 'fn_'.  Change the paths accordingly for your system.
* Run: python EvalData.py

#### Combine variable importance scores across all 200 models:
* Pre-reqs: Unix machine
* Run: bash Step3_GenerateRankedLists.sh



### Methods Comparison

Finally, we evaluated our location predictions across the different
methods by calculating for each cell the mean Euclidean distance of the
top 10 predicted locations from the true location.

    # Comparison of different methods
    source("SubChallenge_1_3/R/Evaluate_Final_Predictions_SUBMIT.R")

#### Post-challenge
After the challenge ended, we evaluated our performance using an outter cross-validation provided by the challenge organizers.  

### Neural Network Approach
The codes for this are included in the 'NeuralNetworks/PostChallenge' folder and utilize the same set of steps as described above.



