###########################################
###########################################
###########################################
###                                     ###
###    JSDM PREDICTION MASTER SCRIPT    ###
###                                     ###
###         SIMULATION DATA ONLY        ###
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

message(sprintf("Job started at %s",
                start_time))

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

dataset_options <- c("sim1random",
                     "sim2random",
                     "sim3random",
                     "sim4random",
                     "sim5random",
                     "sim6random",
                     "sim7random",
                     "sim8random",
                     "sim9random",
                     "sim10random",
                     "sim1spatial",
                     "sim2spatial",
                     "sim3spatial",
                     "sim4spatial",
                     "sim5spatial",
                     "sim6spatial",
                     "sim7spatial",
                     "sim8spatial",
                     "sim9spatial",
                     "sim10spatial")

fold_options <- 1:5

## Extract IDs corresponding to command line argument indices

model_id <- model_options[model_index]

dataset_id <- dataset_options[dataset_index]

fold_id <- fold_options[fold_index]

message("Command line arguments finished")

message(sprintf("Combination:\nModel: %s \nDataset: %s \nFold: %s",
                model_id,
                dataset_id,
                fold_id))

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

library(TruncatedNormal)
library(Metrics)
library(caret)
library(vegan)
library(psych)
library(abind)

message("Packages loaded")

##################
### Run Status ###
##################

## Define run status. 
## Frog data only has two folds so don't run for folds >2

# if(dataset_id == "frog" & fold_id > 2){
#   
#   run_status <- FALSE
#   
# } else {
#   
#   run_status <- TRUE
#   
# }

## Frogs now random five fold so all run_status = TRUE

run_status <- TRUE

message(sprintf("Run status: %s",
                run_status))

########################
### Run Model Script ###
########################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

model_start <- Sys.time()

if(run_status){
  
  command <- sprintf("source('scripts/models/%s_model.R')",
                     model_id)
  
  eval(parse(text = command))
  
} 

message(sprintf("Model fitting duration: %s hours",
                round(difftime(Sys.time(),
                               model_start,
                               units = "hours")[[1]],
                      digits = 5)))

message(sprintf("Total time elapsed: %s hours",
                round(difftime(Sys.time(),
                               start_time,
                               units = "hours")[[1]],
                      digits = 5)))

########################################
### Split up the Bird Runs From Here ###
########################################

## The bird dataset is too computationally intensive to run
## all 1000 posterior samples sequentially for prediction/
## species richness/test statistics.
##
## Here we submit more jobs to Spartan's queue to tackle
## the problem in smaller chunks. 100 runs of 10 samples

if(dataset_id == "bird" & model_id != "SSDM"){
  
  model_options <- c("MPR",
                     "HPR",
                     "LPR",
                     "DPR",
                     "HLR_NS",
                     "HLR_S",
                     "SSDM")
  
  model_command_arg <- which(model_options == model_id)
  
  for(i in seq(1, 1000, 10)){
    
    start_sample <- i
    
    end_sample <- i + 9
    
    command <- sprintf("system('sbatch -p cloud,physical --cpus-per-task=1 --mem=51200 --time=10-00 --job-name=B_%1$s_%2$s_%3$s --output=outputs/slurm_outputs/B_%1$s_%2$s_%3$s.out scripts/slurm/bird_submission.slurm %1$s 3 %2$s %3$s %4$s')",
                       model_command_arg,
                       fold_id,
                       start_sample,
                       end_sample)
    
    eval(parse(text = command))
    
  }
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

prediction_start <- Sys.time()

if(run_status & dataset_id != "bird"){
  
  if(model_id != "SSDM"){
    
    source("scripts/prediction/pmvtnorm_new_function.R")
    source("scripts/prediction/prediction_functions_JSDM.R")
    source("scripts/prediction/prediction_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/prediction/prediction_functions_SESAM.R")
    source("scripts/prediction/prediction_SSDM_SESAM.R")
    
  }
}

if(run_status & dataset_id == "bird" & model_id == "SSDM"){
  
  source("scripts/prediction/prediction_functions_SESAM.R")
  source("scripts/prediction/prediction_SSDM_SESAM.R")
  
}

message(sprintf("Total prediction duration: %s hours",
                round(difftime(Sys.time(),
                               prediction_start,
                               units = "hours")[[1]],
                      digits = 5)))

message(sprintf("Total time elapsed: %s hours",
                round(difftime(Sys.time(),
                               start_time,
                               units = "hours")[[1]],
                      digits = 5)))

#################################
### Run Log-Likelihood Script ###
#################################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

likelihood_start <- Sys.time()

if(run_status){
  
  if(model_id != "SSDM"){
    
    source("scripts/test_statistics/likelihood_function_JSDM.R")
    source("scripts/test_statistics/likelihood_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/test_statistics/likelihood_SSDM.R")
    
  }
}

message(sprintf("Likelihood calculation duration: %s hours",
                round(difftime(Sys.time(),
                               likelihood_start,
                               units = "hours")[[1]],
                      digits = 5)))

message(sprintf("Total time elapsed: %s hours",
                round(difftime(Sys.time(),
                               start_time,
                               units = "hours")[[1]],
                      digits = 5)))

#################################
### Run Test Statistic Script ###
#################################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

test_start <- Sys.time()

if(run_status & dataset_id != "bird"){
  
  if(model_id != "SSDM"){
    
    source("scripts/test_statistics/test_statistic_function.R")
    source("scripts/test_statistics/test_statistics_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/test_statistics/test_statistic_function.R")
    source("scripts/test_statistics/test_statistics_SSDM_SESAM.R")
    
  }
}

if(run_status & dataset_id == "bird" & model_id == "SSDM"){
  
  source("scripts/test_statistics/test_statistic_function.R")
  source("scripts/test_statistics/test_statistics_SSDM_SESAM.R")
  
}

message(sprintf("Test statistic calculation duration: %s hours",
                round(difftime(Sys.time(),
                               test_start,
                               units = "hours")[[1]],
                      digits = 5)))

message(sprintf("Total time elapsed: %s hours",
                round(difftime(Sys.time(),
                               start_time,
                               units = "hours")[[1]],
                      digits = 5)))

###################################
### Run Species Richness Script ###
###################################

## Purge environment for memory's sake. Need to keep 5 constants.

rm(list = ls()[-which(ls() %in% c("model_id",
                                  "dataset_id",
                                  "fold_id",
                                  "run_status",
                                  "start_time"))])

richness_start <- Sys.time()

if(run_status & dataset_id != "bird"){
  
  if(model_id != "SSDM"){
    
    source("scripts/species_richness/species_richness_function.R")
    source("scripts/species_richness/species_richness_JSDM.R")
    
  } else if(model_id == "SSDM"){
    
    source("scripts/species_richness/species_richness_function.R")
    source("scripts/species_richness/species_richness_SSDM_SESAM.R")
    
  }
}

if(run_status & dataset_id == "bird" & model_id == "SSDM"){
  
  source("scripts/species_richness/species_richness_function.R")
  source("scripts/species_richness/species_richness_SSDM_SESAM.R")
  
}

message(sprintf("Species richness calculation duration: %s hours",
                round(difftime(Sys.time(),
                               richness_start,
                               units = "hours")[[1]],
                      digits = 5)))

message(sprintf("Total time elapsed: %s hours",
                round(difftime(Sys.time(),
                               start_time,
                               units = "hours")[[1]],
                      digits = 5)))

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

message("Metadata saved to file")

if(dataset_id == "bird" & model_id != "SSDM"){
  
  message("Bird analysis master script complete\nSplit up prediction jobs still in progress")
  
} else {
  
  message("Analysis complete!")
  
}

message(sprintf("Total time elapsed: %s hours",
                round(difftime(meta_data$end_time,
                               start_time,
                               units = "hours")[[1]],
                      digits = 5)))
