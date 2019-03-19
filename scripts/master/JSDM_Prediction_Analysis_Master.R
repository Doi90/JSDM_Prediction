#######################################################
#######################################################
#######################################################
###                                                 ###
###         JSDM PREDICTION ANALYSIS SCRIPT         ###
###                                                 ###
###   This script is called to perform the analysis ###
### of the JSDM predictions and their test          ###
### statistics. This script is called during the    ###
### sbatch job submission process on Spartan.       ###
### It calculates summary statistics including:     ###
### NA percentages, means, medians, 95% prediction  ###
### interval ranges, etc.                           ###
###                                                 ###
#######################################################
#######################################################
#######################################################

start_time <- Sys.time()

message(sprintf("Job started at %s",
                start_time))

#####################
### Load packages ###
#####################

.libPaths("/home/davidpw/R/lib/3.5")

library(abind)
library(ggplot2)
library(RColorBrewer)

message("Packages loaded")

##################################
### Set combination parameters ###
##################################

## Models

model_options <- c("MPR",
                   "HPR",
                   "LPR",
                   "DPR",
                   "HLR_NS",
                   "HLR_S",
                   "SSDM",
                   "SESAM")

model_order <- c("SSDM",
                 "SESAM",
                 "MPR",
                 "HPR",
                 "LPR",
                 "DPR",
                 "HLR_NS",
                 "HLR_S")

JSDM_models <- model_options[1:6]

SSDM_models <- model_options[7:8]

## Datasets

dataset_options <- c(#"frog",
"eucalypt")#,
# "bird",
# "sim1random",
# "sim2random",
# "sim3random",
# "sim4random",
# "sim5random",
# "sim6random",
# "sim7random",
# "sim8random",
# "sim9random",
# "sim10random",
# "sim1spatial",
# "sim2spatial",
# "sim3spatial",
# "sim4spatial",
# "sim5spatial",
# "sim6spatial",
# "sim7spatial",
# "sim8spatial",
# "sim9spatial",
# "sim10spatial")

## Prediction types

prediction_options <- c("marginal_bin",
                        "marginal_prob",
                        "condLOI",
                        "joint",
                        "SESAM",
                        "SSDM_bin",
                        "SSDM_prob")

prediction_levels <- c("marginal_bin",
                       "marginal_prob",
                       "condLOI_low",
                       "condLOI_med",
                       "condLOI_high",
                       "joint",
                       "SESAM",
                       "SSDM_bin",
                       "SSDM_prob")

binary_predictions <- c("marginal_bin",
                        "condLOI_low",
                        "condLOI_med",
                        "condLOI_high",
                        "joint",
                        "SESAM",
                        "SSDM_bin")

probability_predictions <- c("marginal_prob",
                             "SSDM_prob")

## Test statistics

### By Species

ts_species <- c("AUC",
                "bias",
                "MSE",
                "R2",
                "RMSE",
                "SSE",
                "Pearson",
                "Spearman",
                "Kendall",
                "TP",
                "FP",
                "TN",
                "FN",
                "TPR",
                "FPR",
                "TNR",
                "FNR",
                "PLR",
                "NLR",
                "DOR",
                "Prevalence",
                "Accuracy",
                "PPV",
                "FOR",
                "FDR",
                "NPV",
                "F_1",
                "Youden_J",
                "Kappa")

### By Site

ts_site <- c("Binomial",
             "Bray",
             "Canberra",
             "Euclidean",
             "Gower",
             "Gower_alt",
             "Horn",
             "Jaccard",
             "Kulczynski",
             "Mahalanobis",
             "Manhattan",
             "Mountford",
             "Raup")

#############################
### Load Function Scripts ###
#############################

source("scripts/analysis/proportionNA_function.R")
source("scripts/analysis/notin_function.R")

#######################
#######################
### TEST STATISTICS ###
#######################
#######################

#############################################################
### Create empty dataframe to store test statistic output ###
#############################################################

n_row <- length(model_options) * length(dataset_options) * 
  (length(prediction_options) + 2) *                        # +2 for extra condLOI low/med/high
  (length(ts_site) + length(ts_species))                    # add number of ts together

ts_df <- data.frame(model = factor(character(n_row),
                                   levels = model_order[model_order %in% model_options]),
                    dataset = factor(character(n_row),
                                     levels = dataset_options),
                    prediction_type = factor(character(n_row),
                                             levels = prediction_levels),
                    prediction_class = factor(character(n_row),
                                              levels = c("binary",
                                                         "probability")),
                    test_statistic = factor(character(n_row),
                                            levels = c(ts_species,
                                                       ts_site)),
                    mean = numeric(n_row),
                    median = numeric(n_row),
                    lower = numeric(n_row),
                    upper = numeric(n_row),
                    porportion_NA = numeric(n_row))

########################################
### Summarise Test Statistic Metrics ###
########################################

## Set up a row index counter
## Increase AFTER usage

row_index <- 1

## Loop over model

for(model in model_options){
  
  ## Loop over datasets
  
  for(dataset in dataset_options){
    
    ## Loop over prediction type
    
    for(pred_type in prediction_options){
      
      ########################################
      ### Combination Compatibility Checks ###
      ########################################
      
      ## Skip combinations of loop iterators that are non-compatible
      
      ### SSDM
      
      if(model == "SSDM" & pred_type %notin% c("SSDM_bin", "SSDM_prob")){
        
        next()
        
      }
      
      ### SESAM
      
      if(model == "SESAM" & pred_type != "SESAM"){
        
        next()
        
      }
      
      ### JSDMs
      
      if(model %in% JSDM_models & pred_type %notin% c("marginal_bin",
                                                      "marginal_prob",
                                                      "condLOI",
                                                      "joint")){
        
        next()
        
      }
      
      #################
      ### Load Data ###
      #################
      
      ## SSDM
      
      if(model == "SSDM"){
        
        filename1 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold1_ts.rds",
                             pred_type,
                             dataset)
        
        filename2 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold2_ts.rds",
                             pred_type,
                             dataset)
        
        filename3 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold3_ts.rds",
                             pred_type,
                             dataset)
        
        filename4 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold4_ts.rds",
                             pred_type,
                             dataset)
        
        filename5 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold5_ts.rds",
                             pred_type,
                             dataset)
        
        ts_array <- readRDS(filename1)
        
        tmp1 <- readRDS(filename2)
        
        tmp2 <- readRDS(filename3)
        
        tmp3 <- readRDS(filename4)
        
        tmp4 <- readRDS(filename5)
        
        ts_array$test_statistics_species <- abind(ts_array$test_statistics_species,
                                                  tmp1$test_statistics_species,
                                                  tmp2$test_statistics_species,
                                                  tmp3$test_statistics_species,
                                                  tmp4$test_statistics_species,
                                                  along = 1)
        
        ts_array$test_statistics_site <- abind(ts_array$test_statistics_site,
                                               tmp1$test_statistics_site,
                                               tmp2$test_statistics_site,
                                               tmp3$test_statistics_site,
                                               tmp4$test_statistics_site,
                                               along = 1)
        
      }
      
      ## SESAM
      
      if(model == "SESAM"){
        
        filename1 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold1_ts.rds",
                             pred_type,
                             dataset)
        
        filename2 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold2_ts.rds",
                             pred_type,
                             dataset)
        
        filename3 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold3_ts.rds",
                             pred_type,
                             dataset)
        
        filename4 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold4_ts.rds",
                             pred_type,
                             dataset)
        
        filename5 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold5_ts.rds",
                             pred_type,
                             dataset)
        
        ts_array <- readRDS(filename1)
        
        tmp1 <- readRDS(filename2)
        
        tmp2 <- readRDS(filename3)
        
        tmp3 <- readRDS(filename4)
        
        tmp4 <- readRDS(filename5)
        
        ts_array$test_statistics_species <- abind(ts_array$test_statistics_species,
                                                  tmp1$test_statistics_species,
                                                  tmp2$test_statistics_species,
                                                  tmp3$test_statistics_species,
                                                  tmp4$test_statistics_species,
                                                  along = 1)
        
        ts_array$test_statistics_site <- abind(ts_array$test_statistics_site,
                                               tmp1$test_statistics_site,
                                               tmp2$test_statistics_site,
                                               tmp3$test_statistics_site,
                                               tmp4$test_statistics_site,
                                               along = 1)
        
      }
      
      ## JSDMs
      
      ### Any dataset EXCEPT birds
      
      if(model %in% JSDM_models & dataset != "birds"){
        
        filename1 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold1_%3$s_ts.rds",
                             model,
                             dataset,
                             pred_type)
        
        filename2 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold2_%3$s_ts.rds",
                             model,
                             dataset,
                             pred_type)
        
        filename3 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold3_%3$s_ts.rds",
                             model,
                             dataset,
                             pred_type)
        
        filename4 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold4_%3$s_ts.rds",
                             model,
                             dataset,
                             pred_type)
        
        filename5 <- sprintf("outputs/test_statistics/%1$s_%2$s_fold5_%3$s_ts.rds",
                             model,
                             dataset,
                             pred_type)
        
        ts_array <- readRDS(filename1)
        
        tmp1 <- readRDS(filename2)
        
        tmp2 <- readRDS(filename3)
        
        tmp3 <- readRDS(filename4)
        
        tmp4 <- readRDS(filename5)
        
        if(pred_type != "condLOI"){
          
          ts_array$test_statistics_species <- abind(ts_array$test_statistics_species,
                                                    tmp1$test_statistics_species,
                                                    tmp2$test_statistics_species,
                                                    tmp3$test_statistics_species,
                                                    tmp4$test_statistics_species,
                                                    along = 1)
          
          ts_array$test_statistics_site <- abind(ts_array$test_statistics_site,
                                                 tmp1$test_statistics_site,
                                                 tmp2$test_statistics_site,
                                                 tmp3$test_statistics_site,
                                                 tmp4$test_statistics_site,
                                                 along = 1)
        }
        
        if(pred_type == "condLOI"){
          
          for(i in seq_len(length(ts_array))){
            
            ts_array[[i]]$test_statistics_species <- abind(ts_array[[i]]$test_statistics_species,
                                                           tmp1[[i]]$test_statistics_species,
                                                           tmp2[[i]]$test_statistics_species,
                                                           tmp3[[i]]$test_statistics_species,
                                                           tmp4[[i]]$test_statistics_species,
                                                           along = 1)
            
            ts_array[[i]]$test_statistics_site <- abind(ts_array[[i]]$test_statistics_site,
                                                        tmp1[[i]]$test_statistics_site,
                                                        tmp2[[i]]$test_statistics_site,
                                                        tmp3[[i]]$test_statistics_site,
                                                        tmp4[[i]]$test_statistics_site,
                                                        along = 1)
          }
          
          
        }
        
      }
      
      ### Bird dataset
      
      if(model %in% JSDM_models & dataset == "birds"){
        
        #### Find all files that match pattern
        
        file_list <- intersect(list.files(path = "outputs/test_statistics",
                                          pattern = sprintf("%1$s_%2$s_fold",
                                                            model,
                                                            dataset),
                                          full.names = TRUE),
                               list.files(path = "outputs/test_statistics",
                                          pattern = pred_type,
                                          full.names = TRUE))
        
        #### Read in first file to start array list
        
        ts_array <- readRDS(file_list[1])
        
        #### Read in other files and abind together
        
        if(pred_type != "condLOI"){
          
          for(extra_file in 2:length(file_list)){
            
            tmp_file <- readRDS(file_list[extra_file])
            
            ts_array$test_statistics_species <- abind(ts_array$test_statistics_species,
                                                      tmp_file$test_statistics_species,
                                                      along = 1)
            
            ts_array$test_statistics_site <- abind(ts_array$test_statistics_site,
                                                   tmp_file$test_statistics_site,
                                                   along = 1)
            
          }
        }
        
        if(pred_type == "condLOI"){
          
          for(extra_file in 2:length(file_list)){
            
            tmp_file <- readRDS(file_list[extra_file])
            
            for(i in seq_len(length(tmp_file))){
              
              ts_array[[i]]$test_statistics_species <- abind(ts_array[[i]]$test_statistics_species,
                                                             tmp_file[[i]]$test_statistics_species,
                                                             along = 1)
              
              ts_array[[i]]$test_statistics_site <- abind(ts_array[[i]]$test_statistics_site,
                                                          tmp_file[[i]]$test_statistics_site,
                                                          along = 1)
              
            }
          }
        }
        
      }
      
      ####################################################
      ### Summarise Test Statistics And Fill Dataframe ###
      ####################################################
      
      ## All prediction types that aren't condLOI
      
      if(pred_type != "condLOI"){
        
        ### Test statistics by species
        
        for(ts in dimnames(ts_array$test_statistics_species)[[2]]){
          
          #### Calculate summaries
          
          ts_mean <- mean(ts_array$test_statistics_species[ , ts, ],
                          na.rm = TRUE)
          
          ts_median <- median(ts_array$test_statistics_species[ , ts, ],
                              na.rm = TRUE)
          
          ts_lower <- quantile(ts_array$test_statistics_species[ , ts, ],
                               probs = 0.025,
                               na.rm = TRUE)[[1]]
          
          ts_upper <- quantile(ts_array$test_statistics_species[ , ts, ],
                               probs = 0.975,
                               na.rm = TRUE)[[1]]
          
          ts_NA <- proportionNA(ts_array$test_statistics_species[ , ts, ])
          
          #### Write to dataframe
          
          ts_df[row_index, ] <- list(model,
                                     dataset,
                                     pred_type,
                                     ifelse(pred_type %in% binary_predictions,
                                            "binary",
                                            "probability"),
                                     ts,
                                     ts_mean,
                                     ts_median,
                                     ts_lower,
                                     ts_upper,
                                     ts_NA)
          
          #### Increase row index counter
          
          row_index <- row_index + 1
          
        }
        
        ### Test statistics by site
        
        for(ts in dimnames(ts_array$test_statistics_site)[[2]]){
          
          #### Calculate summaries
          
          ts_mean <- mean(ts_array$test_statistics_site[ , ts, ],
                          na.rm = TRUE)
          
          ts_median <- median(ts_array$test_statistics_site[ , ts, ],
                              na.rm = TRUE)
          
          ts_lower <- quantile(ts_array$test_statistics_site[ , ts, ],
                               probs = 0.025,
                               na.rm = TRUE)[[1]]
          
          ts_upper <- quantile(ts_array$test_statistics_site[ , ts, ],
                               probs = 0.975,
                               na.rm = TRUE)[[1]]
          
          ts_NA <- proportionNA(ts_array$test_statistics_site[ , ts, ])
          
          #### Write to dataframe
          
          ts_df[row_index, ] <- list(model,
                                     dataset,
                                     pred_type,
                                     ifelse(pred_type %in% binary_predictions,
                                            "binary",
                                            "probability"),
                                     ts,
                                     ts_mean,
                                     ts_median,
                                     ts_lower,
                                     ts_upper,
                                     ts_NA)
          
          #### Increase row index counter
          
          row_index <- row_index + 1
          
        }
      }
      
      if(pred_type == "condLOI"){
        
        for(i in seq_len(length(ts_array))){
          
          tmp_array <- ts_array[[i]]
          
          ### Test statistics by species
          
          for(ts in dimnames(tmp_array$test_statistics_species)[[2]]){
            
            #### Calculate summaries
            
            ts_mean <- mean(tmp_array$test_statistics_species[ , ts, ],
                            na.rm = TRUE)
            
            ts_median <- median(tmp_array$test_statistics_species[ , ts, ],
                                na.rm = TRUE)
            
            ts_lower <- quantile(tmp_array$test_statistics_species[ , ts, ],
                                 probs = 0.025,
                                 na.rm = TRUE)[[1]]
            
            ts_upper <- quantile(tmp_array$test_statistics_species[ , ts, ],
                                 probs = 0.975,
                                 na.rm = TRUE)[[1]]
            
            ts_NA <- proportionNA(tmp_array$test_statistics_species[ , ts, ])
            
            #### Write to dataframe
            
            cond_pred_type <- c("condLOI_low", "condLOI_med", "condLOI_high")[i]
            
            ts_df[row_index, ] <- list(model,
                                       dataset,
                                       cond_pred_type,
                                       ifelse(cond_pred_type %in% binary_predictions,
                                              "binary",
                                              "probability"),
                                       ts,
                                       ts_mean,
                                       ts_median,
                                       ts_lower,
                                       ts_upper,
                                       ts_NA)
            
            #### Increase row index counter
            
            row_index <- row_index + 1
            
          }
          
          ### Test statistics by site
          
          for(ts in dimnames(tmp_array$test_statistics_site)[[2]]){
            
            #### Calculate summaries
            
            ts_mean <- mean(tmp_array$test_statistics_site[ , ts, ],
                            na.rm = TRUE)
            
            ts_median <- median(tmp_array$test_statistics_site[ , ts, ],
                                na.rm = TRUE)
            
            ts_lower <- quantile(tmp_array$test_statistics_site[ , ts, ],
                                 probs = 0.025,
                                 na.rm = TRUE)[[1]]
            
            ts_upper <- quantile(tmp_array$test_statistics_site[ , ts, ],
                                 probs = 0.975,
                                 na.rm = TRUE)[[1]]
            
            ts_NA <- proportionNA(tmp_array$test_statistics_site[ , ts, ])
            
            #### Write to dataframe
            
            cond_pred_type <- c("condLOI_low", "condLOI_med", "condLOI_high")[i]
            
            ts_df[row_index, ] <- list(model,
                                       dataset,
                                       cond_pred_type,
                                       ifelse(cond_pred_type %in% binary_predictions,
                                              "binary",
                                              "probability"),
                                       ts,
                                       ts_mean,
                                       ts_median,
                                       ts_lower,
                                       ts_upper,
                                       ts_NA)
            
            #### Increase row index counter
            
            row_index <- row_index + 1
            
          }
        }
      }
    }
    
  }
}

##################################################
### Final Test Statistic Summary Modifications ###
##################################################

## Remove excess rows

ts_df <- ts_df[-(row_index:nrow(ts_df)), ]

## Replace "incorrect" prediction types

# ts_df$prediction_type <-  gsub(x = ts_df$prediction_type,
#                                pattern = "SSDM_bin",
#                                replacement = "marginal_bin")
# 
# ts_df$prediction_type <-  gsub(x = ts_df$prediction_type,
#                                pattern = "SSDM_prob",
#                                replacement = "marginal_prob")
# 
# ts_df$prediction_type <-  gsub(x = ts_df$prediction_type,
#                                pattern = "SESAM",
#                                replacement = "marginal_bin")

###################################################
### Save Test Statistics Summary Output To File ###
###################################################

saveRDS(object = ts_df,
        file = "outputs/test_statistics/test_statistics_summary.rds")

########################
########################
### SPECIES RICHNESS ###
########################
########################

###############################################################
### Create empty dataframe to store species richness output ###
###############################################################

n_row <- length(model_options) * length(dataset_options) * 
  (length(prediction_options) + 2) # +2 for extra condLOI low/med/high

sr_df <- data.frame(model = factor(character(n_row),
                                   levels = model_order[model_order %in% model_options]),
                    dataset = factor(character(n_row),
                                     levels = dataset_options),
                    prediction_type = factor(character(n_row),
                                             levels = prediction_levels),
                    prediction_class = factor(character(n_row),
                                              levels = c("binary",
                                                         "probability")),
                    test_statistic = factor(character(n_row),
                                            levels = "species_richness_difference"),
                    mean = numeric(n_row),
                    median = numeric(n_row),
                    lower = numeric(n_row),
                    upper = numeric(n_row),
                    porportion_NA = numeric(n_row))

##########################################
### Summarise Species Richness Metrics ###
##########################################

## Set up a row index counter
## Increase AFTER usage

row_index <- 1

## Loop over model

for(model in model_options){
  
  ## Loop over datasets
  
  for(dataset in dataset_options){
    
    ## Loop over prediction type
    
    for(pred_type in prediction_options){
      
      ########################################
      ### Combination Compatibility Checks ###
      ########################################
      
      ## Skip combinations of loop iterators that are non-compatible
      
      ### SSDM
      
      if(model == "SSDM" & pred_type %notin% c("SSDM_bin", "SSDM_prob")){
        
        next()
        
      }
      
      ### SESAM
      
      if(model == "SESAM" & pred_type != "SESAM"){
        
        next()
        
      }
      
      ### JSDMs
      
      if(model %in% JSDM_models & pred_type %notin% c("marginal_bin",
                                                      "marginal_prob",
                                                      "condLOI",
                                                      "joint")){
        
        next()
        
      }
      
      #################
      ### Load Data ###
      #################
      
      ## SSDM
      
      if(model == "SSDM"){
        
        filename1 <- sprintf("outputs/species_richness/%1$s_%2$s_fold1_SR.rds",
                             pred_type,
                             dataset)
        
        filename2 <- sprintf("outputs/species_richness/%1$s_%2$s_fold2_SR.rds",
                             pred_type,
                             dataset)
        
        filename3 <- sprintf("outputs/species_richness/%1$s_%2$s_fold3_SR.rds",
                             pred_type,
                             dataset)
        
        filename4 <- sprintf("outputs/species_richness/%1$s_%2$s_fold4_SR.rds",
                             pred_type,
                             dataset)
        
        filename5 <- sprintf("outputs/species_richness/%1$s_%2$s_fold5_SR.rds",
                             pred_type,
                             dataset)
        
        sr_array <- readRDS(filename1)
        
        tmp1 <- readRDS(filename2)
        
        tmp2 <- readRDS(filename3)
        
        tmp3 <- readRDS(filename4)
        
        tmp4 <- readRDS(filename5)
        
        sr_array <- abind(sr_array,
                          tmp1,
                          tmp2,
                          tmp3,
                          tmp4,
                          along = 1)
        
      }
      
      ## SESAM
      
      if(model == "SESAM"){
        
        filename1 <- sprintf("outputs/species_richness/%1$s_%2$s_fold1_SR.rds",
                             model,
                             dataset)
        
        filename2 <- sprintf("outputs/species_richness/%1$s_%2$s_fold2_SR.rds",
                             model,
                             dataset)
        
        filename3 <- sprintf("outputs/species_richness/%1$s_%2$s_fold3_SR.rds",
                             model,
                             dataset)
        
        filename4 <- sprintf("outputs/species_richness/%1$s_%2$s_fold4_SR.rds",
                             model,
                             dataset)
        
        filename5 <- sprintf("outputs/species_richness/%1$s_%2$s_fold5_SR.rds",
                             model,
                             dataset)
        
        sr_array <- readRDS(filename1)
        
        tmp1 <- readRDS(filename2)
        
        tmp2 <- readRDS(filename3)
        
        tmp3 <- readRDS(filename4)
        
        tmp4 <- readRDS(filename5)
        
        sr_array <- abind(sr_array,
                          tmp1,
                          tmp2,
                          tmp3,
                          tmp4,
                          along = 1)
        
      }
      
      ## JSDMs
      
      ### Any dataset EXCEPT birds
      
      if(model %in% JSDM_models & dataset != "birds"){
        
        filename1 <- sprintf("outputs/species_richness/%1$s_%2$s_fold1_%3$s_SR.rds",
                             model,
                             dataset,
                             pred_type)

        filename2 <- sprintf("outputs/species_richness/%1$s_%2$s_fold2_%3$s_SR.rds",
                             model,
                             dataset,
                             pred_type)

        filename3 <- sprintf("outputs/species_richness/%1$s_%2$s_fold3_%3$s_SR.rds",
                             model,
                             dataset,
                             pred_type)

        filename4 <- sprintf("outputs/species_richness/%1$s_%2$s_fold4_%3$s_SR.rds",
                             model,
                             dataset,
                             pred_type)

        filename5 <- sprintf("outputs/species_richness/%1$s_%2$s_fold5_%3$s_SR.rds",
                             model,
                             dataset,
                             pred_type)
        
        sr_array <- readRDS(filename1)
        
        tmp1 <- readRDS(filename2)
        
        tmp2 <- readRDS(filename3)
        
        tmp3 <- readRDS(filename4)
        
        tmp4 <- readRDS(filename5)
        
        sr_array <- abind(sr_array,
                          tmp1,
                          tmp2,
                          tmp3,
                          tmp4,
                          along = 1)
        
      }
      
      ### Bird dataset
      
      if(model %in% JSDM_models & dataset == "birds"){
        
        #### Find all files that match pattern
        
        file_list <- intersect(list.files(path = "outputs/species_richness",
                                          pattern = sprintf("%1$s_%2$s_fold",
                                                            model,
                                                            dataset),
                                          full.names = TRUE),
                               list.files(path = "outputs/species_richness",
                                          pattern = pred_type,
                                          full.names = TRUE))
        
        #### Read in first file to start array list
        
        sr_array <- readRDS(file_list[1])
        
        #### Read in other files and abind together
        
        for(extra_file in 2:length(file_list)){
          
          tmp_file <- readRDS(file_list[extra_file])
          
          sr_array <- abind(sr_array,
                            tmp_file,
                            along = 1)
          
        }
        
      }
      
      #####################################################
      ### Summarise Species Richness And Fill Dataframe ###
      #####################################################
      
      ## All prediction types that aren't condLOI
      
      if(pred_type != "condLOI"){
        
        ### Calculate summaries
        
        sr_mean <- mean(sr_array[ , "Difference", ],
                        na.rm = TRUE)
        
        sr_median <- median(sr_array[ , "Difference", ],
                            na.rm = TRUE)
        
        sr_lower <- quantile(sr_array[ , "Difference", ],
                             probs = 0.025,
                             na.rm = TRUE)[[1]]
        
        sr_upper <- quantile(sr_array[ , "Difference", ],
                             probs = 0.975,
                             na.rm = TRUE)[[1]]
        
        sr_NA <- proportionNA(sr_array[ , "Difference", ])
        
        #### Write to dataframe
        
        sr_df[row_index, ] <- list(model,
                                   dataset,
                                   pred_type,
                                   ifelse(pred_type %in% binary_predictions,
                                          "binary",
                                          "probability"),
                                   "species_richness_difference",
                                   sr_mean,
                                   sr_median,
                                   sr_lower,
                                   sr_upper,
                                   sr_NA)
        
        #### Increase row index counter
        
        row_index <- row_index + 1
        
      }
      
      if(pred_type == "condLOI"){
        
        for(i in seq_len(dim(sr_array)[4])){
          
          tmp_array <- sr_array[ , , , i]
          
          #### Calculate summaries
          
          sr_mean <- mean(tmp_array[ , "Difference", ],
                          na.rm = TRUE)
          
          sr_median <- median(tmp_array[ , "Difference", ],
                              na.rm = TRUE)
          
          sr_lower <- quantile(tmp_array[ , "Difference", ],
                               probs = 0.025,
                               na.rm = TRUE)[[1]]
          
          sr_upper <- quantile(tmp_array[ , "Difference", ],
                               probs = 0.975,
                               na.rm = TRUE)[[1]]
          
          sr_NA <- proportionNA(tmp_array[ , "Difference", ])
          
          #### Write to dataframe
          
          cond_pred_type <- c("condLOI_low", "condLOI_med", "condLOI_high")[i]
          
          sr_df[row_index, ] <- list(model,
                                     dataset,
                                     cond_pred_type,
                                     ifelse(cond_pred_type %in% binary_predictions,
                                            "binary",
                                            "probability"),
                                     "species_richness_difference",
                                     sr_mean,
                                     sr_median,
                                     sr_lower,
                                     sr_upper,
                                     sr_NA)
          
          #### Increase row index counter
          
          row_index <- row_index + 1
          
        }
      }
    }
    
  }
}

####################################################
### Final Species Richness Summary Modifications ###
####################################################

## Remove excess rows

sr_df <- sr_df[-(row_index:nrow(sr_df)), ]

## Replace "incorrect" prediction types
# 
# sr_df$prediction_type <-  gsub(x = sr_df$prediction_type,
#                                pattern = "SSDM_bin",
#                                replacement = "marginal_bin")
# 
# sr_df$prediction_type <-  gsub(x = sr_df$prediction_type,
#                                pattern = "SSDM_prob",
#                                replacement = "marginal_prob")
# 
# sr_df$prediction_type <-  gsub(x = sr_df$prediction_type,
#                                pattern = "SESAM",
#                                replacement = "marginal_bin")

####################################################
### Save Species Richness Summary Output To File ###
####################################################

saveRDS(object = sr_df,
        file = "outputs/species_richness/species_richness_summary.rds")

################
################
### PLOTTING ###
################
################

############################
##### Set Colour Scale #####
############################

col_palette <- brewer.pal(8, "Dark2")

colour <- c()

if("SSDM" %in% model_options){
  colour <- c(colour, col_palette[8])
}
if("SESAM" %in% model_options){
  colour <- c(colour, col_palette[7])
}
if("MPR" %in% model_options){
  colour <- c(colour, col_palette[6])
}
if("HPR" %in% model_options){
  colour <- c(colour, col_palette[5])
}
if("LPR" %in% model_options){
  colour <- c(colour, col_palette[4])
}
if("DPR" %in% model_options){
  colour <- c(colour, col_palette[3])
}
if("HLR_S" %in% model_options){
  colour <- c(colour, col_palette[2])
}
if("HLR_NS" %in% model_options){
  colour <- c(colour, col_palette[1])
}

colour_inverse <- rev(colour)   # because coord_flip() 

###########################################
### Test Statistic-Specific Axis Limits ###
###########################################

graph_customisation <- read.csv("scripts/analysis/test_statistic_axis_limits.csv",
                                stringsAsFactors = FALSE)

#####################################################
### Set Prediction "Pairs" to put JSDMs and SSDMs ###
###       on same plot where applicable           ###
#####################################################

pred_pairs <- list(marginal_bin = c("marginal_bin", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Marginal"),
                   marginal_prob = c("marginal_prob", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Marginal"),
                   condLOI_low = c("condLOI_low", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (low)"),
                   condLOI_med = c("condLOI_med", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (med)"),
                   condLOI_high = c("condLOI_high", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (high)"),
                   joint = c("joint", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Joint"))

############################
######## Make JPEGs ########
############################

dodge <- position_dodge(width = 0.5)  

for(dataset in dataset_options){
  
  for(pair in seq_len(length(pred_pairs))){
    
    for(ts in c(ts_species, ts_site)){
      
      ## Open PDF for editing
      
      pdf(file = sprintf("outputs/test_statistics/plots/test_statistic_plot_%1$s_%2$s_%3$s.pdf",
                         dataset,
                         names(pred_pairs[pair]),
                         ts))
      
      ## Create plot
      
      print(ggplot(ts_df[ts_df$dataset == dataset & 
                           ts_df$test_statistic == ts &
                           ts_df$prediction_type %in% pred_pairs[[pair]][1:3], ],
                   aes(x = model,
                       y = mean,
                       colour = model)) +
              geom_point(position = dodge,
                         size = 2) + 
              geom_errorbar(aes(ymax = upper,
                                ymin = lower,),
                            position = dodge,
                            width = 0.2) +
              ylim(c(graph_customisation[graph_customisation$test_statistic == ts, "y_min"],
                     graph_customisation[graph_customisation$test_statistic == ts, "y_max"])) +
              xlab("Model") +
              ylab(eval(parse(text = graph_customisation[graph_customisation$test_statistic == ts, "graph_label"]))) +
              scale_colour_manual(values = colour,
                                  breaks = levels(ts_df$model)) +
              theme_bw() +
              theme(legend.position = "none",
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) +
              ggtitle(label = pred_pairs[[pair]][4],
                      subtitle = pred_pairs[[pair]][5]) +
              geom_vline(xintercept = 2.5))
        
      
      ## Close PDF and save
      
      dev.off()
      
    }
  }
}



