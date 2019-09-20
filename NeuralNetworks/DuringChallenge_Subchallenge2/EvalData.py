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
fn_predictme = "dge_normalized.txt" # this file has been supplied to us
fn_geometry = "geometry.tsv" # same as the supplied geometry.txt, but as a tab-separated table (instead of space separated)
fn_insitugenes = "bdtnp.genenames.n84.txt" # an 84-line file with the names of the 84 insitu genes

def main():
   random.seed (314159); # set pi seed
   
   for loopindex in range (40): # do 5-fold cross validation 40 times for a total of 200 models
      print ("%d time performing k-folds" % (loopindex))
      # lets do some predictions on our model 
      predicttable = pd.read_csv (fn_predictme, sep = '\t', header = 0, index_col=0, keep_default_na=False, na_values = 'NaN'); # special NaN processing because a feature is called 'nan' which messed up python.  grr!!!
      geomtable = pd.read_csv (fn_geometry, sep = '\t', header = 0);
      
      # train using 5-fold cross validation (base the split on uniq RNASeq predictor variables since there are duplicated labels for the same predictors.  This way there are never values in the training and validation sets that have the exactly same predictor values)
      validationFoldIndexes = (np.arange(predicttable.shape[1]) // ((predicttable.shape[1] // 5) + 1)); # split data (uniq predictor variables) into k=5 pieces for 5-fold cross-feature_validation
      np.random.shuffle (validationFoldIndexes); # make sure the k-folds are split randomly
   
      # iterate across the folds  
      for g, df_subset in predicttable.groupby(validationFoldIndexes, axis=1):
         observationsInValidationFold = list (set (df_subset.columns.tolist()))
         observationsInTrainingFold = list (set (predicttable.columns.tolist()) - set (df_subset.columns.tolist()));
         
         # Load the training table
         dfinput_ExplanatoryVariables, dfinput_ResponseVariable = genDataset ();
   
         # find/remove correlated feature only using the data from the Training fold
         correlated_features_fromTraining = removeCorrelated (predicttable.loc[:,observationsInTrainingFold].transpose(), corthresh=0.6) # get a list of features (from training set only) that were correlated
         dfinput_ExplanatoryVariables.drop (correlated_features_fromTraining, axis=1, inplace=True);
         print ("Removed %s columns that were correlated with insitu genes" % (len (correlated_features_fromTraining)))
         
         # get attributes needed to 0-mean center and unit-variance the RNAseq data.  This should be specific to the training set (not validation) set only.  It will later be used to normalize the training and validation set
         norm_attributes = GetAttributes_MeanCenter_UnitVariance (predicttable.loc[dfinput_ExplanatoryVariables.columns,observationsInTrainingFold].transpose()); # only using the survived features and data from the training set, calculate values needed to nomralization (0-mean center and unit variance) the dataset
         
         
         df_subset_reduced = df_subset.loc[list (dfinput_ExplanatoryVariables.columns),] # remove correlated features from predictions as well
         mymodel = da.TripleRegressor (dfinput_ExplanatoryVariables.columns, # pass in names of explanatory variables.  the number of response variables is pegged at 1 for this use-case
                                      size_hl1=100, # optional argument, just for demonstration purposes
                                      size_hl2=100, # optional argument, just for demonstration purposes
                                      UseStdout=True); # optional argument, just for demonstration purposes
         mymodel.trainme (dfinput_ExplanatoryVariables, 
                          dfinput_ResponseVariable, 
                          df_subset_reduced.T, # these will be the validation set based on the random fold
                          norm_attributes, # mean and std deviation needed to normalize the dataset (based on training split only)
                          num_epochs=20000,
                          earlystop_epochsNoChange=50); # optional argument, just for demonstration purposes
         mymodel.CalcVIP_Gedeon()
         mymodel.SaveModel ('saved_model_data.Iter-%d.Kfold-%d' % (loopindex, g)) # save model for loading latter
         # save order if genes (descending VIP scores) to a file 
         descending_VIPscore_genelist = [x[0] for x in sorted (mymodel.output_logging['VIP_Gedeon']['VIP_Rankings'].items(), key=operator.itemgetter(1), reverse=True)]
         with open('saved_model_data.Iter-%d.Kfold-%d.VIPdescending.txt' % (loopindex, g), 'wt') as fp:
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
