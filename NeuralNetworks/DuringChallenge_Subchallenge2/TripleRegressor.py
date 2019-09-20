import numpy as np
import time
import copy
import pickle
import sys

import torch.nn as nn ## neural net library
import torch.optim as optim # optimization package
import torch.utils.data

class TripleRegressor:
   output_logging = {};
   model = None;
   
   def __init__(self, names_inputneurons, size_hl1=10, size_hl2=10, device=torch.device ('cpu'), UseStdout=True):
      self.device = device;
      self.names_inputneurons = names_inputneurons;
      self.UseStdout = UseStdout;
      
      self.model =  self.__CreateBasicNN (self.names_inputneurons.__len__(), size_hl1, size_hl2);
      
   def __EuclideanLoss (self, input, target):
      return (((input[:,0]-target[:,0])**2)+((input[:,1]-target[:,1])**2)+((input[:,2]-target[:,2])**2)).sqrt().sum() / input.shape[0]
      
   def trainme (self, input_explanatory, input_response, input_valset, normalization_params, num_epochs=10000, earlystop_epochsNoChange=100, size_minibatch=100, datasplit = { "percent_train": 0.80 }):
      dfinput_ExplanatoryVariables = input_explanatory;
      dfinput_ResponseVariables = input_response;
      self.adjustExplanatory = normalization_params
      shuf_uniqrows_validation = np.array (input_valset);

      indexes_train = []
      indexes_validation = []
      for index, row in dfinput_ExplanatoryVariables.iterrows():
         a = (list (row) == shuf_uniqrows_validation)
         
         if (np.array(row) == shuf_uniqrows_validation).all(1).any() == True:
            indexes_validation.append (index)
         else:
            indexes_train.append (index)

      # transform data: mean center with unit variance using stats from training split
      self.__Normalize_MeanCenter_UnitVariance (dfinput_ExplanatoryVariables);
      # np.mean (dfinput_ExplanatoryVariables.apply(np.mean,axis=0))
      # np.mean (dfinput_ExplanatoryVariables.apply(np.std,axis=0))
      
      # perform validation/test/training split.  then load training data in mini-batches
      train_tensor_x = torch.stack ([torch.Tensor(i) for i in np.array (dfinput_ExplanatoryVariables.loc[indexes_train] )]) # transform to torch tensors
      train_tensor_y = torch.stack ([torch.Tensor(i) for i in np.array (dfinput_ResponseVariables.loc[indexes_train] )])
      train_my_dataset = torch.utils.data.TensorDataset(train_tensor_x, train_tensor_y)
      val_tensor_x = torch.stack ([torch.Tensor(i) for i in np.array (dfinput_ExplanatoryVariables.loc[indexes_validation] )]) # transform to torch tensors
      val_tensor_y = torch.stack ([torch.Tensor(i) for i in np.array (dfinput_ResponseVariables.loc[indexes_validation] )])
      val_my_dataset = torch.utils.data.TensorDataset(val_tensor_x, val_tensor_y)
      dataloaders = {}
       
      dataloaders['train'] = torch.utils.data.DataLoader(train_my_dataset, batch_size=size_minibatch, shuffle=True) # create your dataloader (train)
      dataloaders['train_Eval'] = torch.utils.data.DataLoader(train_my_dataset, batch_size=train_my_dataset.__len__()); #  generate training data in one big batch
      dataloaders['val'] = torch.utils.data.DataLoader(val_my_dataset, batch_size=size_minibatch, shuffle=True) # create your dataloader (validation)
      dataloaders['val_Eval'] = torch.utils.data.DataLoader(val_my_dataset, batch_size=val_my_dataset.__len__()); #  generate validating data in one big batch

      # create and train model
      self.output_logging['TrainingAndEvaluationOutput'] = self.__train_model(dataloaders, self.__EuclideanLoss, optim.Adadelta (self.model.parameters ()), num_epochs, earlystop_epochsNoChange) # at least on the test data Adadelta seems to converge faster than regular SGD
   
   # do forward pass.  Input is a Pandas Series of expression values.  Output is a numpy array for the X, Y, Z predictions 
   def getPrediction (self, expressionData, jiggleScale=None):
      self.model.eval()
      
      # adjust inputs (scale) to same scale as during training (based on training data!).  This is a must
      adjust_inputs_with = self.adjustExplanatory.transpose()
      norm_expression  = ((np.array (expressionData) - adjust_inputs_with[0,:]) / adjust_inputs_with[1,:]);
      norm_expression = np.array (norm_expression)[0]; # likely not needed but wanted to keep the data shape consistent as I had it before
      
      if (jiggleScale is not None): # jiggle the unit-variance inputs by a 0-mean normal distribution of std jiggleScale
         jiggles = [ np.random.normal (loc=0, scale=jiggleScale) for i in range (len (norm_expression))];
         norm_expression = (norm_expression+jiggles);
      
      # get output
      unnorm_predictions = self.model (torch.tensor (norm_expression).float ()).detach().numpy();
      return (unnorm_predictions);
     
   # print variable importance scores 
   def CalcVIP_Gedeon (self):
      self.output_logging['VIP_Gedeon'] = self.__VIP_Gedeon (self.names_inputneurons); 
      
   def SaveModel (self, out_fn):
      if (self.model == None):
         print ("Can't save model as it hasn't been created");
         return;
      
      # saving only model parameters per: https://pytorch.org/docs/stable/notes/serialization.html
      torch.save (self.model.state_dict (), "%s.model" % (out_fn))

      with open ("%s.outputlogging.pickle" % (out_fn), 'wb') as fp:      
         pickle.dump (self.output_logging, fp);
         
      with open ("%s.explanatoryadjust.pickle" % (out_fn), 'wb') as fp:      
         pickle.dump (self.adjustExplanatory, fp);
      
   def LoadModel (self, in_fn):
      # load model parameters per https://pytorch.org/docs/stable/notes/serialization.html
      self.model.load_state_dict (torch.load ("%s.model" % (in_fn)));
      
      with open ("%s.outputlogging.pickle" % (in_fn), 'rb') as fp:      
         self.output_logging = pickle.load (fp);
         
      with open ("%s.explanatoryadjust.pickle" % (in_fn), 'rb') as fp:      
         self.adjustExplanatory = pickle.load (fp);
      
   
   def __PrintAndLog (self, input_message, input_dict, input_dict_key):
      if (input_dict_key in input_dict):
         input_dict[input_dict_key] = "%s%s%s" % (input_dict[input_dict_key], "\n", input_message);
      else:
         input_dict[input_dict_key] = input_message;
         
      if (self.UseStdout):
         print ("%s: %s" % (input_dict_key, input_message))

   # print VIP (variable importance scores) for model 
   def __VIP_Gedeon (self, input_names):
      outdict = {'VIP_Rankings': {}}
      
      ## gedeon approach here
      #a. Gedeon paper here: http://users.cecs.anu.edu.au/~Tom.Gedeon/pdfs/ContribDataMinv2.pdf
      #b. Another paper that used the method and had a useful illustration of it:
      #   - https://pdfs.semanticscholar.org/188d/26a005b6aac1448b9c52529b93a186c33685.pdf
      #   - see Fig. 2 and section "Analyzing Weights" for a decent description
      carryOverC = np.asmatrix (np.ones ([3,1])); # 3 output neurons (x,y,z)
      for layer in [4,2,0]: # recurse-backwards through neural layers
         w = np.absolute (list (self.model.parameters())[layer].data.numpy());
         w_sum = np.sum (np.absolute (list (self.model.parameters())[layer].data.numpy()), axis=1);
         c = np.asmatrix ((w/w_sum[:, np.newaxis])); # neuron-level relative weight normalization
         carryOverC = c.transpose () * carryOverC;
      print ()
      self.__PrintAndLog ("VarImp (Gedeon): Attempt to implement Gedeon's 1997 method", outdict, 'VerboseOutput');
      relImp3 = np.asarray (carryOverC/np.sum(carryOverC)).flatten (); # get relative contribute of input neurons (not sure if necessary like it was above)
     
      for i in np.argsort (relImp3)[::-1]:
         self.__PrintAndLog (relImp3[i], outdict['VIP_Rankings'], input_names[i]);
         
      return outdict;

   # create basic neural network model
   def __CreateBasicNN (self, n_InputNeurons, n_hl1, n_hl2):
      return nn.Sequential(
                nn.Dropout(0.1), # 10% hinton dropout on input layers
                nn.Linear(n_InputNeurons, n_hl1),
                nn.ReLU(),
                nn.Dropout(0.1), # 10% hinton dropout on hidden layer 1
                nn.Linear(n_hl1, n_hl2),
                nn.ReLU(),
                nn.Dropout(0.1), # 10% hinton dropout on hidden layer 2
                nn.Linear(n_hl2, 3), # triple-regression - 3 output neurons
              ).to (self.device);
              
   # For all columns (separately) for the passed in data frame: transform data - mean center with unit variance
   def __Normalize_MeanCenter_UnitVariance (self, param_dataframe):
      thecol = 0; 
      for thecolumn in param_dataframe.columns: # not an efficient loop but no time to optimize at the moment
         param_dataframe[thecolumn] = (param_dataframe[thecolumn] - self.adjustExplanatory[thecol,0]) / self.adjustExplanatory[thecol,1]
         thecol += 1
         
        
         

   # training loop originally derived from tutorial at: https://pytorch.org/tutorials/beginner/transfer_learning_tutorial.html
   def __train_model(self, model_data, criterion, optimizer, num_epochs, earlystop_epochsNoChange):
      outdict = {'LossMetricsAfterTraining': {}};
      since = time.time()

      best_model_wts = copy.deepcopy(self.model.state_dict())
      best_val_loss = -1;
      earlystop_counter = 0;

      for epoch in range(num_epochs):
         self.__PrintAndLog ('Epoch {}/{}'.format(epoch, num_epochs - 1), outdict, 'VerboseTrainingOutput');
         self.__PrintAndLog ('-' * 10, outdict, 'VerboseTrainingOutput');
         
         if (earlystop_counter == earlystop_epochsNoChange):
            self.__PrintAndLog ('Early stopping at beginning of epoch {} because validation loss hasn\'t improved in {} cycles'.format (epoch, earlystop_epochsNoChange), outdict, 'VerboseTrainingOutput');
            break;

         # Each epoch has a training and validation phase
         for phase in ['train', 'val']:
            if phase == 'train':
               self.model.train()  # Set model to training mode
            else:
               self.model.eval()   # Set model to evaluate mode

            running_loss = 0.0

            # Iterate over data.
            for inputs, labels in model_data[phase]: # iterate over the various mini-batches for the phase of interest
               inputs = inputs.to(self.device)
               labels = labels.to(self.device)

               # zero the parameter gradients
               optimizer.zero_grad()

               # forward
               # track history if only in train
               with torch.set_grad_enabled(phase == 'train'):
                  outputs = self.model(inputs)
                  loss = criterion(outputs, labels)

                  # backward + optimize only if in training phase
                  if phase == 'train':
                     loss.backward()
                     optimizer.step()

               # statistics
               running_loss += loss.item() * inputs.size(0)
                  
            epoch_loss = running_loss / (model_data[phase].__len__()*model_data[phase].batch_size); # running_loss / dataset_sizes[phase]
            self.__PrintAndLog ('{} Loss: {:.4f}'.format(phase, epoch_loss), outdict, 'VerboseTrainingOutput');

            # deep copy the model
            if phase == 'val':
               if ((epoch_loss < best_val_loss) or (best_val_loss == -1)):
                  best_val_loss = epoch_loss
                  best_model_wts = copy.deepcopy(self.model.state_dict())
                  earlystop_counter = 0; # reset early-stopping counter because validation set had better accuracy then a prior validation set
               else:
                  earlystop_counter += 1;

         self.__PrintAndLog ('', outdict, 'VerboseTrainingOutput');

      time_elapsed = time.time() - since
      self.__PrintAndLog ('Training complete in {:.0f}m {:.0f}s'.format(time_elapsed // 60, time_elapsed % 60), outdict, 'VerboseTrainingOutput');
      self.__PrintAndLog ('Best val Acc: {:.4f}'.format(best_val_loss), outdict, 'VerboseTrainingOutput');
      
      # load best model weights
      self.model.load_state_dict(best_model_wts)

      # generate evaluation loses in one big batch
      self.model.eval()   # Set model to evaluate mode
      
      self.__PrintAndLog ("Generating loss results (training scores may not match above because batchNorm using runningStats on evaluation mode", outdict, 'VerboseTrainingOutput');
      for phase in ['train_Eval', 'val_Eval']:
         for x, y in model_data[phase]:
            self.__PrintAndLog (criterion (self.model (x), y).data.tolist (), outdict['LossMetricsAfterTraining'], phase);
         
      return outdict;
