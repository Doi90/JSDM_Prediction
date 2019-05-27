#######################################################
#######################################################
#######################################################
###                                                 ###
###     JSDM TEST STATISTIC SUMMARY EXTRACTION      ###
###                                                 ###
###   This script extracts the summary values of    ###
### the JSDM prediction test statistic outputs and  ###
### puts them in the format required to do some     ###
### mixed effects modelling.                        ###
###                                                 ###
###   This script is called as part of an sbatch    ###
### submission process so that it runs once per     ###
### dataset.                                        ###
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
library(PassButter)
library(stringr)

message("Packages loaded")

##################################
### Set combination parameters ###
##################################

## Models

model_options <- c("MPR",
                     #"HPR",
                     #"LPR",
                     #"DPR",
                     #"HLR_NS",
                     #"HLR_S",
                   "SSDM")#,
                   #"SESAM")

model_order <- c("SSDM",
                #"SESAM",
                 "MPR")#,
#  "HPR",
#  "LPR",
#  "DPR",
#  "HLR_NS",
#  "HLR_S")

JSDM_models <- model_options[1]#:6]

SSDM_models <- model_options[2]#7:8]

## Datasets

dataset_options <- c("frog")#,
#"eucalypt")#,
# "bird",
# "butterfly",
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

## Folds

fold_options <- 1:5

## Prediction types

prediction_options <- c("marginal_bin",
                        "marginal_prob",
                        "condLOI",
                        "condLOI_marg",
                        "joint",
                        "SESAM",
                        "SSDM_bin",
                        "SSDM_prob")

prediction_levels <- c("marginal_bin",
                       "marginal_prob",
                       "condLOI_low",
                       "condLOI_med",
                       "condLOI_high",
                       "condLOI_marg_low",
                       "condLOI_marg_med",
                       "condLOI_marg_high",
                       "condLOI_marg",      # Joint likelihood doesn't utilise prediction so no L/M/H designation
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
                             "condLOI_marg_low",
                             "condLOI_marg_med",
                             "condLOI_marg_high",
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

## Test statistic / prediction type compatibility

binary_ts <- c("TP",
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
               "Accuracy",
               "PPV",
               "FOR",
               "FDR",
               "NPV",
               "F_1",
               "Youden_J",
               "Kappa",
               "Binomial",
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

prob_ts <- c("AUC",
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
             "Youden_J")

## Chapter

if(length(model_options) == 2){
  
  chapter <- "Ch2"
  
} else {
  
  chapter <- "Ch3"
  
}

#############################
### Load Function Scripts ###
#############################

source("scripts/analysis/proportionNA_function.R")
#source("scripts/analysis/notin_function.R")

#######################
#######################
### TEST STATISTICS ###
#######################
#######################

##############################################################
### Create empty dataframes to store test statistic output ###
##############################################################

ts_df_species <- data.frame(model = factor(character(),
                                           levels = model_order[model_order %in% model_options]),
                            dataset = factor(character(),
                                             levels = dataset_options),
                            fold = numeric(),
                            species = character(),
                            prediction_type = factor(character(),
                                                     levels = prediction_levels),
                            prediction_class = factor(character(),
                                                      levels = c("binary",
                                                                 "probability")),
                            test_statistic = factor(character(),
                                                    levels = ts_species),
                            mean = numeric(),
                            stringsAsFactors = FALSE)

ts_df_site <- data.frame(model = factor(character(),
                                        levels = model_order[model_order %in% model_options]),
                         dataset = factor(character(),
                                          levels = dataset_options),
                         fold = numeric(),
                         site = character(),
                         prediction_type = factor(character(),
                                                  levels = prediction_levels),
                         prediction_class = factor(character(),
                                                   levels = c("binary",
                                                              "probability")),
                         test_statistic = factor(character(),
                                                 levels = ts_species),
                         mean = numeric(),
                         stringsAsFactors = FALSE)

## Loop over datasets

for(dataset in dataset_options){
  
  ## Loop over folds
  
  for(fold in fold_options){
    
    ## Loop over prediction types
    
    for(pred_type in prediction_options){
      
      ## Loop over models
      
      for(model in model_options){
        
        ########################################
        ### Combination Compatibility Checks ###
        ########################################
        
        ## Skip combinations of loop iterators that are non-compatible
        
        ### SSDM
        
        if(model == "SSDM" & pred_type %nin% c("SSDM_bin", "SSDM_prob")){
          
          next()
          
        }
        
        ### SESAM
        
        if(model == "SESAM" & pred_type != "SESAM"){
          
          next()
          
        }
        
        ### JSDMs
        
        if(model %in% JSDM_models & pred_type %nin% c("marginal_bin",
                                                      "marginal_prob",
                                                      "condLOI",
                                                      "condLOI_marg",
                                                      "joint")){
          
          next()
          
        }
        
        ######################################
        ### Load test statistic posteriors ###
        ######################################
        
        ## SSDMs
        
        if(model == "SSDM"){
          
          filename <- sprintf("outputs/test_statistics/%1$s_%2$s_fold%3$s_ts.rds",
                              pred_type,
                              dataset,
                              fold)
          
          ts_array <- readRDS(filename)
          
        }
        
        ## SESAM
        
        if(model == "SESAM"){
          
          filename <- sprintf("outputs/test_statistics/%1$s_%2$s_fold%3$s_ts.rds",
                              pred_type,
                              dataset,
                              fold)
          
          ts_array <- readRDS(filename)
          
        }
        
        ## JSDMs
        
        ### For non-bird datasets
        
        if(model %in% JSDM_models & dataset != "birds"){
          
          filename <- sprintf("outputs/test_statistics/%1$s_%2$s_fold%3$s_%4$s_ts.rds",
                              model,
                              dataset,
                              fold,
                              pred_type)
          
          ts_array <- readRDS(filename)
          
        }
        
        ## For the bird dataset
        
        if(model %in% JSDM_models & dataset == "birds"){
          
          ### Find all files that match pattern
          
          file_list <- intersect(list.files(path = "outputs/test_statistics",
                                            pattern = sprintf("%1$s_%2$s_fold%3$s",
                                                              model,
                                                              dataset,
                                                              fold),
                                            full.names = TRUE),
                                 list.files(path = "outputs/test_statistics",
                                            pattern = pred_type,
                                            full.names = TRUE))
          
          ### Read in first file to start array list
          
          ts_array <- readRDS(file_list[1])
          
          ### Read in other files and abind together
          
          if(pred_type != "condLOI"){
            
            for(extra_file in 2:length(file_list)){
              
              tmp_file <- readRDS(file_list[extra_file])
              
              ts_array$test_statistics_species <- abind(ts_array$test_statistics_species,
                                                        tmp_file$test_statistics_species,
                                                        along = 3)
              
              ts_array$test_statistics_site <- abind(ts_array$test_statistics_site,
                                                     tmp_file$test_statistics_site,
                                                     along = 3)
              
            }
          }
          
          if(pred_type == "condLOI"){
            
            for(extra_file in 2:length(file_list)){
              
              tmp_file <- readRDS(file_list[extra_file])
              
              for(i in seq_len(length(tmp_file))){
                
                ts_array[[i]]$test_statistics_species <- abind(ts_array[[i]]$test_statistics_species,
                                                               tmp_file[[i]]$test_statistics_species,
                                                               along = 3)
                
                ts_array[[i]]$test_statistics_site <- abind(ts_array[[i]]$test_statistics_site,
                                                            tmp_file[[i]]$test_statistics_site,
                                                            along = 3)
                
              }
            }
          }
          
        }
        
        ###################################################
        ### Create temporary dataframes to store ouputs ###
        ###################################################
        
        if(pred_type %nin% c("condLOI", "condLOI_marg")){
          
          n_species <- dim(ts_array$test_statistics_species)[1]
          
          n_site <- dim(ts_array$test_statistics_site)[1]
          
        }
        
        if(pred_type %in% c("condLOI", "condLOI_marg")){
          
          n_species <- dim(ts_array[[1]]$test_statistics_species)[1]
          
          n_site <- dim(ts_array[[1]]$test_statistics_site)[1]
          
        }
        
        tmp_species <- data.frame(model = factor(character(n_species * length(ts_species)),
                                                 levels = model_order[model_order %in% model_options]),
                                  dataset = factor(character(n_species * length(ts_species)),
                                                   levels = dataset_options),
                                  fold = numeric(n_species * length(ts_species)),
                                  species = character(n_species * length(ts_species)),
                                  prediction_type = factor(character(n_species * length(ts_species)),
                                                           levels = prediction_levels),
                                  prediction_class = factor(character(n_species * length(ts_species)),
                                                            levels = c("binary",
                                                                       "probability")),
                                  test_statistic = factor(character(n_species * length(ts_species)),
                                                          levels = ts_species),
                                  mean = numeric(n_species * length(ts_species)),
                                  stringsAsFactors = FALSE)
        
        tmp_site <- data.frame(model = factor(character(n_site * length(ts_site)),
                                              levels = model_order[model_order %in% model_options]),
                               dataset = factor(character(n_site * length(ts_site)),
                                                levels = dataset_options),
                               fold = numeric(n_site * length(ts_site)),
                               site = character(n_site * length(ts_site)),
                               prediction_type = factor(character(n_site * length(ts_site)),
                                                        levels = prediction_levels),
                               prediction_class = factor(character(n_site * length(ts_site)),
                                                         levels = c("binary",
                                                                    "probability")),
                               test_statistic = factor(character(n_site * length(ts_site)),
                                                       levels = ts_site),
                               mean = numeric(n_site * length(ts_site)),
                               stringsAsFactors = FALSE)
        
        #########################
        ### Extract Summaries ###
        #########################
        
        if(pred_type %nin% c("condLOI", "condLOI_marg")){
          
          row_index_spp <- 1
          
          for(spp in seq_len(dim(ts_array$test_statistics_species)[1])){
            
            for(ts in ts_species){
              
              tmp_species[row_index_spp, ] <- list(model,
                                                   dataset,
                                                   fold,
                                                   sprintf("%s_%s",
                                                           dataset,
                                                           spp),
                                                   pred_type,
                                                   ifelse(pred_type %in% binary_predictions,
                                                          "binary",
                                                          "probability"),
                                                   ts,
                                                   mean(ts_array$test_statistics_species[spp, ts, ][is.finite(ts_array$test_statistics_species[spp, ts, ])],
                                                        na.rm = TRUE))
              
              row_index_spp <- row_index_spp + 1
              
            }
          }
          
          row_index_site <- 1
          
          for(site in seq_len(dim(ts_array$test_statistics_site)[1])){
            
            for(ts in ts_site){
              
              tmp_site[row_index_site, ] <- list(model,
                                                 dataset,
                                                 fold,
                                                 sprintf("%s_fold%s_%s",
                                                         dataset,
                                                         fold,
                                                         site),
                                                 pred_type,
                                                 ifelse(pred_type %in% binary_predictions,
                                                        "binary",
                                                        "probability"),
                                                 ts,
                                                 mean(ts_array$test_statistics_site[site, ts, ][is.finite(ts_array$test_statistics_site[site, ts, ])],
                                                      na.rm = TRUE))
              
              row_index_site <- row_index_site + 1
              
            }
          }
          
        }
        
        if(pred_type == "condLOI"){
          
          tmp_species <- rbind(tmp_species,
                               tmp_species,
                               tmp_species)
          
          tmp_site <- rbind(tmp_site,
                            tmp_site,
                            tmp_site)
          
          cond_pred_type <- c("condLOI_low",
                              "condLOI_med",
                              "condLOI_high")
          
          row_index_species <- 1
          
          for(cond in seq_len(length(ts_array))){
            
            for(spp in seq_len(dim(ts_array[[cond]]$test_statistics_species)[1])){
              
              for(ts in ts_species){
                
                tmp_species[row_index_species, ] <- list(model,
                                                         dataset,
                                                         fold,
                                                         sprintf("%s_%s",
                                                                 dataset,
                                                                 spp),
                                                         cond_pred_type[cond],
                                                         ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                                "binary",
                                                                "probability"),
                                                         ts,
                                                         mean(ts_array[[cond]]$test_statistics_species[spp, ts, ][is.finite(ts_array[[cond]]$test_statistics_species[spp, ts, ])],
                                                              na.rm = TRUE))
                
                row_index_species <- row_index_species + 1
                
              }
              
            }
          }
          
          row_index_site <- 1
          
          for(cond in seq_len(length(ts_array))){
            
            for(site in seq_len(dim(ts_array[[cond]]$test_statistics_site)[1])){
              
              for(ts in ts_site){
                
                tmp_site[row_index_site, ] <- list(model,
                                                   dataset,
                                                   fold,
                                                   sprintf("%s_fold%s_%s",
                                                           dataset,
                                                           fold,
                                                           site),
                                                   cond_pred_type[cond],
                                                   ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                          "binary",
                                                          "probability"),
                                                   ts,
                                                   mean(ts_array[[cond]]$test_statistics_site[site, ts, ][is.finite(ts_array[[cond]]$test_statistics_site[site, ts, ])],
                                                        na.rm = TRUE))
                
                row_index_site <- row_index_site + 1
                
              }
              
            }
          }
          
        }
        
        if(pred_type == "condLOI_marg"){
          
          tmp_species <- rbind(tmp_species,
                               tmp_species,
                               tmp_species)
          
          tmp_site <- rbind(tmp_site,
                            tmp_site,
                            tmp_site)
          
          cond_pred_type <- c("condLOI_marg_low",
                              "condLOI_marg_med",
                              "condLOI_marg_high")
          
          row_index_species <- 1
          
          for(cond in seq_len(length(ts_array))){
            
            for(spp in seq_len(dim(ts_array[[cond]]$test_statistics_species)[1])){
              
              for(ts in ts_species){
                
                tmp_species[row_index_species, ] <- list(model,
                                                         dataset,
                                                         fold,
                                                         sprintf("%s_%s",
                                                                 dataset,
                                                                 spp),
                                                         cond_pred_type[cond],
                                                         ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                                "binary",
                                                                "probability"),
                                                         ts,
                                                         mean(ts_array[[cond]]$test_statistics_species[spp, ts, ][is.finite(ts_array[[cond]]$test_statistics_species[spp, ts, ])],
                                                              na.rm = TRUE))
                
                row_index_species <- row_index_species + 1
                
              }
              
            }
          }
          
          row_index_site <- 1
          
          for(cond in seq_len(length(ts_array))){
            
            for(site in seq_len(dim(ts_array[[cond]]$test_statistics_site)[1])){
              
              for(ts in ts_site){
                
                tmp_site[row_index_site, ] <- list(model,
                                                   dataset,
                                                   fold,
                                                   sprintf("%s_fold%s_%s",
                                                           dataset,
                                                           fold,
                                                           site),
                                                   cond_pred_type[cond],
                                                   ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                          "binary",
                                                          "probability"),
                                                   ts,
                                                   mean(ts_array[[cond]]$test_statistics_site[site, ts, ][is.finite(ts_array[[cond]]$test_statistics_site[site, ts, ])],
                                                        na.rm = TRUE))
                
                row_index_site <- row_index_site + 1
                
              }
              
            }
          }
          
        }
        
        #########################################
        ### Add to complete dataframe storage ###
        #########################################
        
        ts_df_species <- rbind(ts_df_species,
                               tmp_species)
        
        ts_df_site <- rbind(ts_df_site,
                            tmp_site)
        
        
      }
    }
  }
}

##################################################
### Final Test Statistic Summary Modifications ###
##################################################

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

## Correct the "site" IDs

### Load site IDS from file

for(dataset in dataset_options){
  
  command <- sprintf("%1$s_site_ids <- readRDS('data/%1$s/site_ids.rds')",
                     dataset)
  
  eval(parse(text = command))
  
}

### Perform conversion

for(i in seq_len(nrow(ts_df_site))){
  
  old_id <- ts_df_site[i, "site"]
  
  old_id_facets <- str_split(old_id, "_")[[1]]
  
  new_id <- sprintf("%s_site_ids[[%s]][%s]", 
                    old_id_facets[1],
                    str_split(old_id_facets[2], "")[[1]][5],
                    old_id_facets[3])
  
  ts_df_site[i, "site"] <- as.character(eval(parse(text = new_id)))
  
}

## Fold as factor

ts_df_species$fold <- as.factor(ts_df_species$fold)

ts_df_site$fold <- as.factor(ts_df_site$fold)

###################################################
### Save Test Statistics Summary Output To File ###
###################################################

saveRDS(object = ts_df_species,
        file = sprintf("outputs/test_statistics/test_statistics_species_summary_MME_%s.rds",
                       chapter))

saveRDS(object = ts_df_site,
        file = sprintf("outputs/test_statistics/test_statistics_site_summary_MME_%s.rds",
                       chapter))

########################
########################
### SPECIES RICHNESS ###
########################
########################

###############################################################
### Create empty dataframe to store species richness output ###
###############################################################

sr_df <- data.frame(model = factor(character(),
                                   levels = model_order[model_order %in% model_options]),
                    dataset = factor(character(),
                                     levels = dataset_options),
                    fold = numeric(),
                    site = character(),
                    prediction_type = factor(character(),
                                             levels = prediction_levels),
                    prediction_class = factor(character(),
                                              levels = c("binary",
                                                         "probability")),
                    test_statistic = factor(character(),
                                            levels = ts_species),
                    mean = numeric(),
                    stringsAsFactors = FALSE)

##########################################
### Summarise Species Richness Metrics ###
##########################################

## Loop over model

for(model in model_options){
  
  ## Loop over fold
  
  for(fold in fold_options){
    
    ## Loop over datasets
    
    for(dataset in dataset_options){
      
      ## Loop over prediction type
      
      for(pred_type in prediction_options){
        
        ########################################
        ### Combination Compatibility Checks ###
        ########################################
        
        ## Skip combinations of loop iterators that are non-compatible
        
        ### SSDM
        
        if(model == "SSDM" & pred_type %nin% c("SSDM_bin", "SSDM_prob")){
          
          next()
          
        }
        
        ### SESAM
        
        if(model == "SESAM" & pred_type != "SESAM"){
          
          next()
          
        }
        
        ### JSDMs
        
        if(model %in% JSDM_models & pred_type %nin% c("marginal_bin",
                                                      "marginal_prob",
                                                      "condLOI",
                                                      "condLOI_marg",
                                                      "joint")){
          
          next()
          
        }
        
        #################
        ### Load Data ###
        #################
        
        ## SSDM
        
        if(model == "SSDM"){
          
          filename <- sprintf("outputs/species_richness/%1$s_%2$s_fold%3$s_SR.rds",
                              pred_type,
                              dataset,
                              fold)
          
          sr_array <- readRDS(filename)
          
        }
        
        ## SESAM
        
        if(model == "SESAM"){
          
          filename <- sprintf("outputs/species_richness/%1$s_%2$s_fold%3$s_SR.rds",
                              model,
                              dataset,
                              fold)
          
          sr_array <- readRDS(filename)
          
        }
        
        ## JSDMs
        
        ### Any dataset EXCEPT birds
        
        if(model %in% JSDM_models & dataset != "birds"){
          
          filename <- sprintf("outputs/species_richness/%1$s_%2$s_fold%3$s_%4$s_SR.rds",
                              model,
                              dataset,
                              fold,
                              pred_type)
          
          sr_array <- readRDS(filename)
          
        }
        
        ### Bird dataset
        
        if(model %in% JSDM_models & dataset == "birds"){
          
          #### Find all files that match pattern
          
          file_list <- intersect(list.files(path = "outputs/species_richness",
                                            pattern = sprintf("%1$s_%2$s_fold%3$s",
                                                              model,
                                                              dataset,
                                                              fold),
                                            full.names = TRUE),
                                 list.files(path = "outputs/species_richness",
                                            pattern = pred_type,
                                            full.names = TRUE))
          
          #### Read in first file to start array list
          
          sr_array <- readRDS(file_list[1])
          
          #### Read in other files and abind together
          
          if(pred_type != "condLOI"){
            
            for(extra_file in 2:length(file_list)){
              
              tmp_file <- readRDS(file_list[extra_file])
              
              sr_array <- abind(sr_array,
                                tmp_file,
                                along = 1)
              
            }
          }
          
          if(pred_type == "condLOI"){
            
            for(extra_file in 2:length(file_list)){
              
              tmp_file <- readRDS(file_list[extra_file])
              
              for(i in seq_len(length(tmp_file))){
                
                sr_array[[i]] <- abind(sr_array[[i]],
                                       tmp_file[[i]],
                                       along = 3)
                
              }
            }
          }
        }
        
        ###################################################
        ### Create temporary dataframes to store ouputs ###
        ###################################################
        
        n_site <- dim(sr_array)[1]
        
        tmp_sr <- data.frame(model = factor(character(n_site),
                                            levels = model_order[model_order %in% model_options]),
                             dataset = factor(character(n_site),
                                              levels = dataset_options),
                             fold = numeric(n_site),
                             site = character(n_site),
                             prediction_type = factor(character(n_site),
                                                      levels = prediction_levels),
                             prediction_class = factor(character(n_site),
                                                       levels = c("binary",
                                                                  "probability")),
                             test_statistic = factor(character(n_site),
                                                     levels = "species_richness_difference"),
                             mean = numeric(n_site),
                             stringsAsFactors = FALSE)
        
        #####################################################
        ### Summarise Species Richness And Fill Dataframe ###
        #####################################################
        
        ## All prediction types that aren't condLOI
        
        if(pred_type %nin% c("condLOI", "condLOI_marg")){
          
          row_index <- 1
          
          for(j in seq_len(n_site)){
            
            tmp_sr[row_index, ] <- list(model,
                                        dataset,
                                        fold,
                                        sprintf("%s_fold%s_%s",
                                                dataset,
                                                fold,
                                                j),
                                        pred_type,
                                        ifelse(pred_type %in% binary_predictions,
                                               "binary",
                                               "probability"),
                                        "species_richness_difference",
                                        mean(sr_array[ j, "Difference", ],
                                             na.rm = TRUE))
            
            row_index <- row_index + 1
            
          }
        }
        
        if(pred_type == "condLOI"){
          
          tmp_sr <- rbind(tmp_sr,
                          tmp_sr,
                          tmp_sr)
          
          cond_pred_type <- c("condLOI_low",
                              "condLOI_med",
                              "condLOI_high")
          
          row_index <- 1
          
          for(cond in seq_len(dim(sr_array)[4])){
            
            tmp_array <- sr_array[ , , , cond]
            
            for(j in seq_len(n_site)){
              
              tmp_sr[row_index, ] <- list(model,
                                          dataset,
                                          fold,
                                          sprintf("%s_fold%s_%s",
                                                  dataset,
                                                  fold,
                                                  j),
                                          cond_pred_type[cond],
                                          ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                 "binary",
                                                 "probability"),
                                          "species_richness_difference",
                                          mean(sr_array[ j, "Difference", , cond],
                                               na.rm = TRUE))
              
              row_index <- row_index + 1
              
            }
          }
        }
        
        if(pred_type == "condLOI_marg"){
          
          tmp_sr <- rbind(tmp_sr,
                          tmp_sr,
                          tmp_sr)
          
          cond_pred_type <- c("condLOI_marg_low",
                              "condLOI_marg_med",
                              "condLOI_marg_high")
          
          row_index <- 1
          
          for(cond in seq_len(dim(sr_array)[4])){
            
            tmp_array <- sr_array[ , , , cond]
            
            for(j in seq_len(n_site)){
              
              tmp_sr[row_index, ] <- list(model,
                                          dataset,
                                          fold,
                                          sprintf("%s_fold%s_%s",
                                                  dataset,
                                                  fold,
                                                  j),
                                          cond_pred_type[cond],
                                          ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                 "binary",
                                                 "probability"),
                                          "species_richness_difference",
                                          mean(sr_array[ j, "Difference", , cond],
                                               na.rm = TRUE))
              
              row_index <- row_index + 1
              
            }
          }
        }
        
        #########################################
        ### Add to complete dataframe storage ###
        #########################################
        
        sr_df <- rbind(sr_df,
                       tmp_sr)
        
      }
    }
  }
}

####################################################
### Final Species Richness Summary Modifications ###
####################################################

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

## Correct the "site" IDs

### Load site IDS from file

for(dataset in dataset_options){
  
  command <- sprintf("%1$s_site_ids <- readRDS('data/%1$s/site_ids.rds')",
                     dataset)
  
  eval(parse(text = command))
  
}

### Perform conversion

for(i in seq_len(nrow(sr_df))){
  
  old_id <- sr_df[i, "site"]
  
  old_id_facets <- str_split(old_id, "_")[[1]]
  
  new_id <- sprintf("%s_site_ids[[%s]][%s]", 
                    old_id_facets[1],
                    str_split(old_id_facets[2], "")[[1]][5],
                    old_id_facets[3])
  
  sr_df[i, "site"] <- as.character(eval(parse(text = new_id)))
  
}

## Fold as factor

sr_df$fold <- as.factor(sr_df$fold)

####################################################
### Save Species Richness Summary Output To File ###
####################################################

saveRDS(object = sr_df,
        file = sprintf("outputs/species_richness/species_richness_summary_%s.rds",
                       chapter))

######################
######################
### LOG-LIKELIHOOD ###
######################
######################

#############################################################
### Create empty dataframe to store log-likelihood output ###
#############################################################

ll_i_df <- data.frame(model = factor(character(),
                                     levels = model_order[model_order %in% model_options]),
                      dataset = factor(character(),
                                       levels = dataset_options),
                      fold = numeric(),
                      site = character(),
                      prediction_type = factor(character(),
                                               levels = prediction_levels),
                      prediction_class = factor(character(),
                                                levels = c("binary",
                                                           "probability")),
                      test_statistic = factor(character(),
                                              levels = ts_species),
                      mean = numeric(),
                      stringsAsFactors = FALSE)

ll_j_df <- data.frame(model = factor(character(),
                                     levels = model_order[model_order %in% model_options]),
                      dataset = factor(character(),
                                       levels = dataset_options),
                      fold = numeric(),
                      species = character(),
                      prediction_type = factor(character(),
                                               levels = prediction_levels),
                      prediction_class = factor(character(),
                                                levels = c("binary",
                                                           "probability")),
                      test_statistic = factor(character(),
                                              levels = ts_species),
                      mean = numeric(),
                      stringsAsFactors = FALSE)

##########################################
### Summarise log-likelihood Metrics ###
##########################################

## Loop over model

for(model in model_options){
  
  ## Loop over fold
  
  for(fold in fold_options){
    
    ## Loop over datasets
    
    for(dataset in dataset_options){
      
      ## Loop over prediction type
      
      for(pred_type in prediction_options){
        
        ########################################
        ### Combination Compatibility Checks ###
        ########################################
        
        ## Skip combinations of loop iterators that are non-compatible
        
        if(pred_type %nin% c("SSDM_prob", "marginal_prob", "condLOI_marg")){
          
          next()
          
        }
        
        ### SSDM
        
        if(model == "SSDM" & pred_type != "SSDM_prob"){
          
          next()
          
        }
        
        ### SESAM
        
        if(model == "SESAM" & pred_type != "SESAM"){
          
          next()
          
        }
        
        ### JSDMs
        
        if(model %in% JSDM_models & pred_type %nin% c("marginal_prob",
                                                      "condLOI_marg")){
          
          next()
          
        }
        
        ## Pred type name check
        
        pred_name <- ifelse(pred_type == "marginal_prob",
                            "marg",
                            pred_type)
        
        #################
        ### Load Data ###
        #################
        
        ## SSDM
        
        if(model == "SSDM"){
          
          filename_i <- sprintf("outputs/likelihood/%1$s_%2$s_fold%3$s_independent_likelihood.rds",
                                pred_type,
                                dataset,
                                fold)
          
          ll_i_array <- readRDS(filename_i)
          
          filename_j <- sprintf("outputs/likelihood/%1$s_%2$s_fold%3$s_joint_likelihood.rds",
                                pred_type,
                                dataset,
                                fold)
          
          ll_j_array <- readRDS(filename_j)
          
        }
        
        ## JSDMs
        
        ### Any dataset EXCEPT birds
        
        if(model %in% JSDM_models & dataset != "birds"){
          
          filename_i <- sprintf("outputs/likelihood/%1$s_%2$s_fold%3$s_%4$s_independent_likelihood.rds",
                                model,
                                dataset,
                                fold,
                                pred_name)
          
          ll_i_array <- readRDS(filename_i)
          
          filename_j <- sprintf("outputs/likelihood/%1$s_%2$s_fold%3$s_%4$s_joint_likelihood.rds",
                                model,
                                dataset,
                                fold,
                                pred_name)
          
          ll_j_array <- readRDS(filename_j)
          
        }
        
        # ### Bird dataset
        # 
        # if(model %in% JSDM_models & dataset == "birds"){
        #   
        #   #### Find all files that match pattern
        #   
        #   file_list <- intersect(list.files(path = "outputs/likelihood",
        #                                     pattern = sprintf("%1$s_%2$s_fold%3$s",
        #                                                       model,
        #                                                       dataset,
        #                                                       fold),
        #                                     full.names = TRUE),
        #                          list.files(path = "outputs/likelihood",
        #                                     pattern = pred_type,
        #                                     full.names = TRUE))
        #   
        #   #### Read in first file to start array list
        #   
        #   sr_array <- readRDS(file_list[1])
        #   
        #   #### Read in other files and abind together
        #   
        #   if(pred_type != "condLOI"){
        #     
        #     for(extra_file in 2:length(file_list)){
        #       
        #       tmp_file <- readRDS(file_list[extra_file])
        #       
        #       sr_array <- abind(sr_array,
        #                         tmp_file,
        #                         along = 1)
        #       
        #     }
        #   }
        #   
        #   if(pred_type == "condLOI"){
        #     
        #     for(extra_file in 2:length(file_list)){
        #       
        #       tmp_file <- readRDS(file_list[extra_file])
        #       
        #       for(i in seq_len(length(tmp_file))){
        #         
        #         sr_array[[i]] <- abind(sr_array[[i]],
        #                                tmp_file[[i]],
        #                                along = 3)
        #         
        #       }
        #     }
        #   }
        # }
        
        ###################################################
        ### Create temporary dataframes to store ouputs ###
        ###################################################
        
        if(model %in% JSDM_models){
          
          n_spp <- dim(ll_i_array)[1]
          
          n_site <- dim(ll_j_array)[1]
          
          tmp_ll_j <- data.frame(model = factor(character(n_site),
                                                levels = model_order[model_order %in% model_options]),
                                 dataset = factor(character(n_site),
                                                  levels = dataset_options),
                                 fold = numeric(n_site),
                                 site = character(n_site),
                                 prediction_type = factor(character(n_site),
                                                          levels = prediction_levels),
                                 prediction_class = factor(character(n_site),
                                                           levels = c("binary",
                                                                      "probability")),
                                 test_statistic = factor(character(n_site),
                                                         levels = "joint_log_likelihood"),
                                 mean = numeric(n_site),
                                 stringsAsFactors = FALSE)
          
          tmp_ll_i <- data.frame(model = factor(character(n_spp),
                                                levels = model_order[model_order %in% model_options]),
                                 dataset = factor(character(n_spp),
                                                  levels = dataset_options),
                                 fold = numeric(n_spp),
                                 species = character(n_spp),
                                 prediction_type = factor(character(n_spp),
                                                          levels = prediction_levels),
                                 prediction_class = factor(character(n_spp),
                                                           levels = c("binary",
                                                                      "probability")),
                                 test_statistic = factor(character(n_spp),
                                                         levels = "independent_log_likelihood"),
                                 mean = numeric(n_spp),
                                 stringsAsFactors = FALSE)
          
          
          ###############################################
          ### Summarise Likelihood And Fill Dataframe ###
          ###############################################
          
          ## Joint log-likelihood
          
          row_index_j <- 1
          
          for(j in seq_len(n_site)){
            
            tmp_ll_j[row_index_j, ] <- list(model,
                                            dataset,
                                            fold,
                                            sprintf("%s_fold%s_%s",
                                                    dataset,
                                                    fold,
                                                    j),
                                            pred_type,
                                            ifelse(pred_type %in% binary_predictions,
                                                   "binary",
                                                   "probability"),
                                            "joint_log_likelihood",
                                            mean(ll_j_array[j, ],
                                                 na.rm = TRUE))
            
            row_index_j <- row_index_j + 1
            
          }
          
          ## Indepedent log-likelihood
          
          ### All prediction types that aren't condLOI
          
          if(pred_type != "condLOI_marg"){
            
            row_index_i <- 1
            
            for(i in seq_len(n_spp)){
              
              tmp_ll_i[row_index_i, ] <- list(model,
                                              dataset,
                                              fold,
                                              sprintf("%s_fold%s_%s",
                                                      dataset,
                                                      fold,
                                                      i),
                                              pred_type,
                                              ifelse(pred_type %in% binary_predictions,
                                                     "binary",
                                                     "probability"),
                                              "independent_log_likelihood",
                                              mean(ll_i_array[i, ],
                                                   na.rm = TRUE))
              
              row_index_i <- row_index_i + 1
              
            }
          }
          
          if(pred_type == "condLOI_marg"){
            
            tmp_ll_i <- rbind(tmp_ll_i,
                              tmp_ll_i,
                              tmp_ll_i)
            
            cond_pred_type <- c("condLOI_marg_low",
                                "condLOI_marg_med",
                                "condLOI_marg_high")
            
            row_index_i <- 1
            
            for(cond in seq_len(dim(ll_i_array)[3])){
              
              tmp_array <- ll_i_array[ , , cond]
              
              for(i in seq_len(n_spp)){
                
                tmp_ll_i[row_index_i, ] <- list(model,
                                                dataset,
                                                fold,
                                                sprintf("%s_fold%s_%s",
                                                        dataset,
                                                        fold,
                                                        i),
                                                cond_pred_type[cond],
                                                ifelse(cond_pred_type[cond] %in% binary_predictions,
                                                       "binary",
                                                       "probability"),
                                                "independent_log_likelihood",
                                                mean(tmp_array[i, ],
                                                     na.rm = TRUE))
                
                row_index_i <- row_index_i + 1
                
              }
            }
          }
        }
        
        if(model %in% SSDM_models){
          
          n_spp <- length(ll_i_array)
          
          n_site <- length(ll_j_array)
          
          tmp_ll_j <- data.frame(model = factor(character(n_site),
                                                levels = model_order[model_order %in% model_options]),
                                 dataset = factor(character(n_site),
                                                  levels = dataset_options),
                                 fold = numeric(n_site),
                                 site = character(n_site),
                                 prediction_type = factor(character(n_site),
                                                          levels = prediction_levels),
                                 prediction_class = factor(character(n_site),
                                                           levels = c("binary",
                                                                      "probability")),
                                 test_statistic = factor(character(n_site),
                                                         levels = "joint_log_likelihood"),
                                 mean = numeric(n_site),
                                 stringsAsFactors = FALSE)
          
          tmp_ll_i <- data.frame(model = factor(character(n_spp),
                                                levels = model_order[model_order %in% model_options]),
                                 dataset = factor(character(n_spp),
                                                  levels = dataset_options),
                                 fold = numeric(n_spp),
                                 species = character(n_spp),
                                 prediction_type = factor(character(n_spp),
                                                          levels = prediction_levels),
                                 prediction_class = factor(character(n_spp),
                                                           levels = c("binary",
                                                                      "probability")),
                                 test_statistic = factor(character(n_spp),
                                                         levels = "independent_log_likelihood"),
                                 mean = numeric(n_spp),
                                 stringsAsFactors = FALSE)
          
          
          ###############################################
          ### Summarise Likelihood And Fill Dataframe ###
          ###############################################
          
          ## Joint log-likelihood
          
          row_index_j <- 1
          
          for(j in seq_len(n_site)){
            
            tmp_ll_j[row_index_j, ] <- list(model,
                                            dataset,
                                            fold,
                                            sprintf("%s_fold%s_%s",
                                                    dataset,
                                                    fold,
                                                    j),
                                            pred_type,
                                            ifelse(pred_type %in% binary_predictions,
                                                   "binary",
                                                   "probability"),
                                            "joint_log_likelihood",
                                            ll_j_array[j])
            
            row_index_j <- row_index_j + 1
            
          }
          
          ## Indepedent log-likelihood
          
          row_index_i <- 1
          
          for(i in seq_len(n_spp)){
            
            tmp_ll_i[row_index_i, ] <- list(model,
                                            dataset,
                                            fold,
                                            sprintf("%s_fold%s_%s",
                                                    dataset,
                                                    fold,
                                                    i),
                                            pred_type,
                                            ifelse(pred_type %in% binary_predictions,
                                                   "binary",
                                                   "probability"),
                                            "independent_log_likelihood",
                                            ll_i_array[i])
            
            row_index_i <- row_index_i + 1
            
          }
        }
        
        #########################################
        ### Add to complete dataframe storage ###
        #########################################
        
        ll_i_df <- rbind(ll_i_df,
                         tmp_ll_i)
        
        ll_j_df <- rbind(ll_j_df,
                         tmp_ll_j)
        
      }
    }
  }
}

##################################################
### Final Log likelihood Summary Modifications ###
##################################################

## Correct the "site" IDs

### Load site IDS from file

for(dataset in dataset_options){
  
  command <- sprintf("%1$s_site_ids <- readRDS('data/%1$s/site_ids.rds')",
                     dataset)
  
  eval(parse(text = command))
  
}

### Perform conversion

for(i in seq_len(nrow(ll_j_df))){
  
  old_id <- ll_j_df[i, "site"]
  
  old_id_facets <- str_split(old_id, "_")[[1]]
  
  new_id <- sprintf("%s_site_ids[[%s]][%s]", 
                    old_id_facets[1],
                    str_split(old_id_facets[2], "")[[1]][5],
                    old_id_facets[3])
  
  ll_j_df[i, "site"] <- as.character(eval(parse(text = new_id)))
  
}

## Fold as factor

ll_i_df$fold <- as.factor(ll_i_df$fold)

ll_j_df$fold <- as.factor(ll_j_df$fold)

###################################################
### Save Log Likelihood Summary Outputs To File ###
###################################################

saveRDS(object = ll_i_df,
        file = sprintf("outputs/likelihood/independent_likelihood_summary_%s.rds",
                       chapter))

saveRDS(object = ll_j_df,
        file = sprintf("outputs/likelihood/joint_likelihood_summary_%s.rds",
                       chapter))

# ################
# ################
# ### PLOTTING ###
# ################
# ################
# 
# ############################
# ##### Set Colour Scale #####
# ############################
# 
# col_palette <- brewer.pal(8, "Dark2")
# 
# colour <- c()
# 
# if("SSDM" %in% model_options){
#   colour <- c(colour, col_palette[8])
# }
# if("SESAM" %in% model_options){
#   colour <- c(colour, col_palette[7])
# }
# if("MPR" %in% model_options){
#   colour <- c(colour, col_palette[6])
# }
# if("HPR" %in% model_options){
#   colour <- c(colour, col_palette[5])
# }
# if("LPR" %in% model_options){
#   colour <- c(colour, col_palette[4])
# }
# if("DPR" %in% model_options){
#   colour <- c(colour, col_palette[3])
# }
# if("HLR_S" %in% model_options){
#   colour <- c(colour, col_palette[2])
# }
# if("HLR_NS" %in% model_options){
#   colour <- c(colour, col_palette[1])
# }
# 
# colour_inverse <- rev(colour)   # because coord_flip() 
# 
# ###########################################
# ### Test Statistic-Specific Axis Limits ###
# ###########################################
# 
# graph_customisation <- read.csv("scripts/analysis/test_statistic_axis_limits.csv",
#                                 stringsAsFactors = FALSE)
# 
# #####################################################
# ### Set Prediction "Pairs" to put JSDMs and SSDMs ###
# ###       on same plot where applicable           ###
# #####################################################
# 
# pred_pairs <- list(marginal_bin = c("marginal_bin", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Marginal"),
#                    marginal_prob = c("marginal_prob", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Marginal"),
#                    condLOI_low = c("condLOI_low", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (low)"),
#                    condLOI_med = c("condLOI_med", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (med)"),
#                    condLOI_high = c("condLOI_high", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (high)"),
#                    joint = c("joint", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Joint"))
# 
# ############################
# ######## Make JPEGs ########
# ############################
# 
# dodge <- position_dodge(width = 0.5)  
# 
# for(dataset in dataset_options){
#   
#   for(pair in seq_len(length(pred_pairs))){
#     
#     for(ts in c(ts_species, ts_site)){
#       
#       ## Open PDF for editing
#       
#       pdf(file = sprintf("outputs/test_statistics/plots/test_statistic_plot_%1$s_%2$s_%3$s.pdf",
#                          dataset,
#                          names(pred_pairs[pair]),
#                          ts))
#       
#       ## Create plot
#       
#       print(ggplot(ts_df[ts_df$dataset == dataset & 
#                            ts_df$test_statistic == ts &
#                            ts_df$prediction_type %in% pred_pairs[[pair]][1:3], ],
#                    aes(x = model,
#                        y = mean,
#                        colour = model)) +
#               geom_point(position = dodge,
#                          size = 2) + 
#               geom_errorbar(aes(ymax = upper,
#                                 ymin = lower,),
#                             position = dodge,
#                             width = 0.2) +
#               ylim(c(graph_customisation[graph_customisation$test_statistic == ts, "y_min"],
#                      graph_customisation[graph_customisation$test_statistic == ts, "y_max"])) +
#               xlab("Model") +
#               ylab(eval(parse(text = graph_customisation[graph_customisation$test_statistic == ts, "graph_label"]))) +
#               scale_colour_manual(values = colour,
#                                   breaks = levels(ts_df$model)) +
#               theme_bw() +
#               theme(legend.position = "none",
#                     panel.grid.minor = element_blank(),
#                     panel.grid.major = element_blank()) +
#               ggtitle(label = pred_pairs[[pair]][4],
#                       subtitle = pred_pairs[[pair]][5]) +
#               geom_vline(xintercept = 2.5))
#       
#       
#       ## Close PDF and save
#       
#       dev.off()
#       
#     }
#   }
# }
# 
# 
# 
