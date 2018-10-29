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

start_time <- Sys.time()

##################################
### Command line arguments and ###
###       defining indices     ###
##################################

## Read in the command line arguments

command_args <- commandArgs(trailingOnly = TRUE)

model_index <- as.numeric(command_args[1])
dataset_index <- as.numeric(command_args[2])
fold_index <- as.numeric(command_args[3])

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

fold_options <- 1:5

## Extract IDs corresponding to command line argument indices

model_id <- model_options[model_index]

dataset_id <- dataset_options[dataset_index]

fold_id <- fold_options[fold_index]

#####################
### Load packages ###
#####################

.libPaths("/home/davidpw/R/lib/3.5")

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
library(psych)

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

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

if(run_status){
  
command <- sprintf("source('scripts/models/%s_model.R')",
                   model_id)

eval(parse(text = command))

} 

##############################
### Run Prediction Scripts ###
##############################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

if(run_status){
  
  if(model_id != "SSDM"){
    
    source("scripts/prediction/prediction_functions_JSDM.R")
    source("scripts/prediction/prediction_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/prediction/prediction_functions_SESAM.R")
    source("scripts/prediction/prediction_SSDM_SESAM.R")
    
  }
}

#################################
### Run Log-Likelihood Script ###
#################################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

if(run_status){
  
  if(model_id != "SSDM"){
    
    source("scripts/test_statistics/likelihood_function_JSDM.R")
    source("scripts/test_statistics/likelihood_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/test_statistics/likelihood_SSDM.R")
    
  }
}

#################################
### Run Test Statistic Script ###
#################################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

if(run_status){
  
  if(model_id != "SSDM"){
    
    source("scripts/test_statistics/test_statistic_function.R")
    source("scripts/test_statistics/test_statistics_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/test_statistics/test_statistic_function.R")
    source("scripts/test_statistics/test_statistics_SSDM_SESAM.R")
    
  }
}

###################################
### Run Species Richness Script ###
###################################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

if(run_status){
  
  if(model_id != "SSDM"){
    
    source("scripts/species_richness/species_richness_function.R")
    source("scripts/species_richness/species_richness_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/species_richness/species_richness_function.R")
    source("scripts/species_richness/species_richness_SSDM_SESAM.R")
    
  }
}

##############
### Output ###
##############

## Model posteriors saved in model fitting script
##   to outputs/posteriors
## Prediction posteriors saved in prediction script
##   to outputs/predictions
## Test statistics saved in test statistic script
##   to outputs/test_statistics

meta_data <- list(model = model_id,
                  dataset = dataset_id,
                  fold = fold_id,
                  run_status = run_status,
                  session_info = sessionInfo(),
                  start_time = start_time,
                  end_time = Sys.time())

filename <- sprintf("outputs/meta_data/%s_%s_fold%s_metadata.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(meta_data,
        filename)
