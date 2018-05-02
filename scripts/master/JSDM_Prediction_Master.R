###########################################
###########################################
###########################################
###                                     ###
###    JSDM PREDICTION MASTER SCRIPT    ###
###                                     ###
###  This is the script called during   ###
### sbatch job submission on Spartan.   ###
### This script reads in command line   ###
### arguments from the job submission   ###
### process that serve as indices for   ###
### model, dataset, and CV fold. This   ###
### allows the master script to call    ###
### the appropriate combination of      ###
### scripts for the analysis            ###
###                                     ###
###########################################
###########################################
###########################################

##################################
### Command line arguments and ###
###       defining indices     ###
##################################

## Read in the command line arguments

command_args <- commandArgs(trailingOnly = TRUE)

model_index <- command_args[1]
dataset_index <- command_args[2]
fold_index <- command_args[3]

## Possible ID options

model_options <- c("MPR",
                   "HPR",
                   "LPR",
                   "DPR",
                   "HLR_NS",
                   "HLR_S",
                   "SSDM")

dataset_options <- c("frog",
                     "eucalypt",
                     "bird")

fold_options <- c(1:5)

## Extract IDs corresponding to command line argument indices

model_id <- model_options[model_index]

dataset_id <- model_options[model_index]

fold_id <- model_options[model_index]

#####################
### Load packages ###
#####################

### Model script packages

if(model_id == "MPR"){
  
  library(BayesComm)
  library(coda)
  
} else if(model_id == "HPR"){
  
  library(R2jags)
  library(parallel)
  library(random)
  library(abind)
  library(MCMCpack)
  library(MASS)
  library(mclust)
  
} else if(model_id == "LPR"){
  
  library(boral)
  library(coda)
  
} else if(model_id == "DPR"){
  
  library(gjam)
  
} else if(model_id == "HLR_NS" | model_id == "HLR_S"){
  
  library(HMSC)
  library(coda)
  
} else if(model_id == "SSDM"){
  
  library(coda)
  
}

### Prediction packages

library(tmvtnorm)

### Test statistics packages

library(Metrics)
library(caret)
library(vegan)

##################
### Run Status ###
##################

## Define run status. 
## Frog data only has two folds so don't run for folds >2

if(dataset_id == "frog" & fold_id > 2){
  run_status <- FALSE
} else {
  run_status <- TRUE
}

########################
### Run Model Script ###
########################

if(run_status){
  
command <- sprintf("source('scripts/models/%s_model.R')",
                   model_id)

eval(parse(text = command))

} 

##############################
### Run Prediction Scripts ###
##############################

if(run_status){
  
  if(model_id != "SSDM"){
    
    source("scripts/prediction/prediction_functions_JSDM.R")
    source("scripts/prediction/prediction_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source()
    source()
  }
}

#################################
### Run Test Statistic Script ###
#################################

if(run_status){
  
  source("scripts/test_statistics/test_statistics.R")

}