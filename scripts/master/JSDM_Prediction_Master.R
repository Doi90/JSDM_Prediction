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
                   "HLR-S",
                   "SSDM")

dataset_options <- c("Frog",
                     "Eucalypt",
                     "Bird")

fold_options <- c(1:5)

## Extract IDs corresponding to command line argument indices

model_id <- model_options[model_index]

dataset_id <- model_options[model_index]

fold_id <- model_options[model_index]

## Run the correct model script

command <- sprintf("source('scripts/models/%s_model.R')",
                   model_id)

eval(parse(text = command))

## Run the prediction script


## Run the test statistic script