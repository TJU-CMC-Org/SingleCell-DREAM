import numpy as np
import random
import pandas as pd
import json;
import operator;
import pickle;
import csv;

import torch.nn as nn ## neural net library
import torch.optim as optim # optimization package
import torch.utils.data

import TripleRegressor as da
fn_trainingdata = "trainingTable.txt"; # generated using script Step1_GenTrainingDataMatrix.R
fn_predictme = "../../Data/dge_normalized.txt" # this file has been supplied to us
fn_geometry = "../../Data/geometry.tsv" # same as the supplied geometry.txt, but as a tab-separated table (instead of space separated)
fn_insitugenes = "../DuringChallenge_Subchallenge2/bdtnp.genenames.n84.txt" # an 84-line file with the names of the 84 insitu genes
fn_10cv = "10cv.txt" # representative example of 10cv's given by challenge organizers

def main():
   random.seed (314159); # set pi seed
   
   cv10 = pd.read_csv (fn_10cv, sep = '\t', header = None, index_col=None); # special NaN processing because a feature is called 'nan' which messed up python.  grr!!!
   for loopindex in range (20): # 20 iterations (total of 1000 models because 20 iterations * 10-static cv * 5-innercv)
      print ("%d time performing k-folds" % (loopindex))
      geomtable = pd.read_csv (fn_geometry, sep = '\t', header = 0);
      
      # train using supplied 10-fold crossvalidation-ish provided to us for paper
      for foldcv10 in range (cv10.shape[0]):
         # in this case we are going to do a cvfold of 5 within a cvfold of 10 (crazy huh?).  But this way I can use the validation set without touching the test set ;-)
         predicttable_full = pd.read_csv (fn_predictme, sep = '\t', header = 0, index_col=0, keep_default_na=False, na_values = 'NaN'); # special NaN processing because a feature is called 'nan' which messed up python.  grr!!!
         list_train = [(int (x)-1) for x in cv10.iloc[foldcv10,].tolist() if not np.isnan(x)] # python is 0-indexes and input file is 1-indexed, convert
      
         predicttable_train = predicttable_full.iloc[:,list_train]
         predicttable_val = predicttable_full.iloc[:,list (set (list (range (predicttable_full.columns.__len__()))) - set (list_train))]
         
         if (len (set (predicttable_train.columns).intersection (set (predicttable_val.columns))) != 0):
            print ("Error, overlap in train and validation - should never get here")
            exit (1);
      
         # train using 5-fold cross validation (base the split on uniq RNASeq predictor variables since there are duplicated labels for the same predictors.  This way there are never values in the training and validation sets that have the exactly same predictor values)
         validationFoldIndexes = (np.arange(predicttable_train.shape[1]) // ((predicttable_train.shape[1] // 5) + 1)); # split data (uniq predictor variables) into k=5 pieces for 5-fold cross-feature_validation
         np.random.shuffle (validationFoldIndexes); # make sure the k-folds are split randomly
      
         # iterate across the folds  
         for g, df_subset in predicttable_train.groupby(validationFoldIndexes, axis=1):
            
            
            observationsInValidationFold = list (set (df_subset.columns.tolist())) # cell names in INNER cross-validation validation fold
            observationsInTrainingFold = list (set (predicttable_train.columns.tolist()) - set (df_subset.columns.tolist())); # cell names in INNER cross-validation training fold
            
            # Load the training table
            dfinput_ExplanatoryVariables, dfinput_ResponseVariable = genDataset (); # load full dataset, not yet subselected according to outter or inner folds
      
            # find/remove correlated feature only using the data from the Training fold
            correlated_features_fromTraining = removeCorrelated (predicttable_train.loc[:,observationsInTrainingFold].transpose(), corthresh=0.6) # get a list of features/genes (from INNER-fold training set only) that were correlated
            dfinput_ExplanatoryVariables.drop (correlated_features_fromTraining, axis=1, inplace=True); # remove correlated genes from full dataset
            print ("Removed %s columns that were correlated with insitu genes" % (len (correlated_features_fromTraining)))
            
            # get attributes needed to 0-mean center and unit-variance the RNAseq data.  This should be specific to the INNER training set (not validation) set only.  It will later be used to normalize the training and validation set
            norm_attributes = GetAttributes_MeanCenter_UnitVariance (predicttable_train.loc[dfinput_ExplanatoryVariables.columns,observationsInTrainingFold].transpose()); # only using the survived features and data from the training set, calculate values needed to nomralization (0-mean center and unit variance) the dataset
            
            
            df_subset_reduced = df_subset.loc[list (dfinput_ExplanatoryVariables.columns),] # remove correlated features from predictions as well
            innertraining = predicttable_train.loc[dfinput_ExplanatoryVariables.columns,observationsInTrainingFold];
            
            if ((set (predicttable_train.columns) == set(list (innertraining.columns) + list (df_subset_reduced.columns))) != True):
              print ("ERROR, wrong columns!  the INNERFOLD training and validation folds should match the OUTTERFOLD training set, otherwise something is wrong") 
              exit (1)
              
            if (df_subset_reduced.index.equals (innertraining.index) != True):
              print ("ERROR, wrong rows/genes!  the INNERFOLD training and validation folds should match the OUTTERFOLD training set, otherwise something is wrong") 
            
            
            mymodel = da.TripleRegressor (dfinput_ExplanatoryVariables.columns, # pass in names of explanatory variables.  the number of response variables is pegged at 1 for this use-case
                                         size_hl1=100, # optional argument, just for demonstration purposes
                                         size_hl2=100, # optional argument, just for demonstration purposes
                                         UseStdout=True); # optional argument, just for demonstration purposes
            mymodel.trainme (dfinput_ExplanatoryVariables,  # full dataset (subselected later) with correlated features dropped
                             dfinput_ResponseVariable,      # labels of full dataset (subselected later)
                             df_subset_reduced.T, # these will be the INNER validation set 
                             innertraining.T,     # these will be the INNER training set 
                             norm_attributes, # mean and std deviation needed to normalize the dataset (based on training split only)
                             num_epochs=20000,
                             earlystop_epochsNoChange=50); # optional argument, just for demonstration purposes
            mymodel.CalcVIP_Gedeon()
            mymodel.SaveModel ('saved_model_data.Iter-%d.Outterfold-%d.Kfold-%d' % (loopindex, foldcv10, g)) # save model for loading latter
            # save order if genes (descending VIP scores) to a file 
            descending_VIPscore_genelist = [x[0] for x in sorted (mymodel.output_logging['VIP_Gedeon']['VIP_Rankings'].items(), key=operator.itemgetter(1), reverse=True)]
            with open('saved_model_data.Iter-%d.Outterfold-%d.Kfold-%d.VIPdescending.txt' % (loopindex, foldcv10, g), 'wt') as fp:
               fp.write ("\n".join (descending_VIPscore_genelist));
         
# look for genes that are correlated with each other, and then remove one of them strategically if they are
# remove them in place
def removeCorrelated (df_explanatory, corthresh=0.7):
   protectedgenes = pd.read_csv (fn_insitugenes, sep = '\t', header = None); # never throw out insitu genes
   protectedgenes_bool = df_explanatory.columns.isin(protectedgenes[0])
   protectedgenes_index = np.where (df_explanatory.columns.isin(protectedgenes[0]))[0]  # get indexes of protected genes
   
   numpymatrix = np.matrix (df_explanatory)
   pearson_coefs = np.corrcoef(numpymatrix, rowvar=False)
  
   removeindexes = set() 
   for row in range (pearson_coefs.shape[0]): # iterate over each row
      if row in protectedgenes_index: # don't remove any gene indexes that are protected
         removeme = np.where ((pearson_coefs[row].__abs__()>=corthresh) & (~protectedgenes_bool))[0] # find indexes that both aren't protected and exceed correlation threshold
         for removeelement in removeme:
            removeindexes.add (int (removeelement));
   
   toremove_colnames = df_explanatory.columns[list (removeindexes)] 
   return toremove_colnames;

def GetAttributes_MeanCenter_UnitVariance (param_dataframe):
   normStats = np.matrix (np.zeros (shape=(param_dataframe.columns.__len__(), 2))) # columns are (mean & std)
   data_means = np.matrix (param_dataframe.apply(np.mean, axis=0)).transpose ()
   data_stds = np.matrix (param_dataframe.apply(np.std, axis=0)).transpose ()
   normStats[:,0] = data_means;
   normStats[:,1] = data_stds;
      
   return normStats

def genDataset ():
   inputtable = pd.read_csv (fn_trainingdata, sep = '\t', header = 0);
   dfinput_explanatory_inputfeatures = inputtable.iloc[:,3:];
   dfoutput_response_transfection =  inputtable.iloc[:,:3];
   
   return dfinput_explanatory_inputfeatures,  dfoutput_response_transfection;

if __name__ == '__main__':
   main()
