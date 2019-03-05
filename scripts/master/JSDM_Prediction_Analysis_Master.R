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

message("Packages loaded")

##################################
### Set combination parameters ###
##################################

## Models

model_options <- c("MPR",
                   # "HPR",
                   # "LPR",
                   # "DPR",
                   # "HLR_NS",
                   # "HLR_S",
                   "SSDM")#,
                   #"SESAM")

JSDM_models <- model_options[1]

SSDM_models <- model_options[2]

## Datasets

dataset_options <- c("frog")#,
                     # "eucalypt",
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

## Folds

fold_options <- 1:5

## Prediction types

prediction_options <- c("marginal_bin",
                        "marginal_prob",
                        "condLOI",
                        "joint",
                        "SESAM",
                        "SSDM_bin",
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
  length(fold_options) * (length(prediction_options) + 2) * # +2 for extra condLOI low/med/high
  (length(ts_site) + length(ts_species))                    # add number of ts together

ts_df <- data.frame(model = character(n_row),
                    dataset = character(n_row),
                    fold = numeric(n_row),
                    prediction_type = character(n_row),
                    test_statistic = character(n_row),
                    mean = numeric(n_row),
                    median = numeric(n_row),
                    lower = numeric(n_row),
                    upper = numeric(n_row),
                    porportion_NA = numeric(n_row),
                    stringsAsFactors = FALSE)

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
    
    ## Loop over fold
    
    for(fold in fold_options){
      
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
          
          filename <- sprintf("outputs/test_statistics/%1$s_%2$s_fold%3$s_ts.rds",
                              pred_type,
                              dataset,
                              fold)
          
          ts_array <- readRDS(filename)
          
        }
        
        ## SESAM
        
        if(model == "SESAM"){
          
          filename <- sprintf("outputs/test_statistics/%1$s_%2$s_fold%3$s_ts.rds",
                              model,
                              dataset,
                              fold)
          
          ts_array <- readRDS(filename)
          
        }
        
        ## JSDMs
        
        ### Any dataset EXCEPT birds
        
        if(model %in% JSDM_models & dataset != "birds"){
          
          filename <- sprintf("outputs/test_statistics/%1$s_%2$s_fold%3$s_%4$s_ts.rds",
                              model,
                              dataset,
                              fold,
                              pred_type)
          
          ts_array <- readRDS(filename)
          
        }
        
        ### Bird dataset
        
        if(model %in% JSDM_models & dataset == "birds"){
          
          #### Find all files that match pattern
          
          partial_filename <- sprintf("%1$s_%2$s_fold%3$s_%4$s",
                                      model,
                                      dataset,
                                      fold,
                                      pred_type)
          
          file_list <- list.files(path = "outputs/test_statistics",
                                  pattern = partial_filename,
                                  full.names = TRUE)
          
          #### Read in first file to start array list
          
          ts_array <- readRDS(file_list[1])
          
          #### Read in other files and abind together
          
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
                                       fold,
                                       pred_type,
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
                                       fold,
                                       pred_type,
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
                                         fold,
                                         cond_pred_type,
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
                                         fold,
                                         cond_pred_type,
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
  length(fold_options) * (length(prediction_options) + 2) # +2 for extra condLOI low/med/high

sr_df <- data.frame(model = character(n_row),
                    dataset = character(n_row),
                    fold = numeric(n_row),
                    prediction_type = character(n_row),
                    test_statistic = character(n_row),
                    mean = numeric(n_row),
                    median = numeric(n_row),
                    lower = numeric(n_row),
                    upper = numeric(n_row),
                    porportion_NA = numeric(n_row),
                    stringsAsFactors = FALSE)

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
    
    ## Loop over fold
    
    for(fold in fold_options){
      
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
          
          partial_filename <- sprintf("%1$s_%2$s_fold%3$s_%4$s",
                                      model,
                                      dataset,
                                      fold,
                                      pred_type)
          
          file_list <- list.files(path = "outputs/species_richness",
                                  pattern = partial_filename,
                                  full.names = TRUE)
          
          #### Read in first file to start array list
          
          sr_array <- readRDS(file_list[1])
          
          #### Read in other files and abind together
          
          for(extra_file in 2:length(file_list)){
            
            tmp_file <- readRDS(file_list[extra_file])
            
            sr_array <- abind(sr_array,
                              tmp_file,
                              along = 3)
            
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
                                     fold,
                                     pred_type,
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
                                       fold,
                                       cond_pred_type,
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

model <- unique(df$Model)
col_palette <- brewer.pal(6, "Dark2")
colour <- c()

if("MPR" %in% model){
  colour <- c(colour, col_palette[6])
}
if("HPR" %in% model){
  colour <- c(colour, col_palette[5])
}
if("LPR" %in% model){
  colour <- c(colour, col_palette[4])
}
if("DPR" %in% model){
  colour <- c(colour, col_palette[3])
}
if("HLR-S" %in% model){
  colour <- c(colour, col_palette[2])
}
if("HLR-NS" %in% model){
  colour <- c(colour, col_palette[1])
}

colour_inverse <- rev(colour)   # because coord_flip() 


############################
######## Make JPEGs ########
############################

dodge <- position_dodge(width=0.5)  

for(i in unique(df$Species)){
  
  # Create blank PNG file
  
  file.name <- paste("Beta_", i, ".pdf", sep = "")
  
  # Create plot
  
  ggplot(df[df$Species == i,],aes(x = Coefficient, y = Posterior.Mean,
                                  colour = Model)) + 
    geom_point(position=dodge) +
    geom_errorbar(aes(ymax=Upper,ymin=Lower),position = dodge) +
    theme(legend.position = "right") +
    ylab("Posterior Mean") + 
    xlab("Variable") +
    coord_flip() +
    ggtitle(i) +
    scale_colour_manual(values = colour_inverse, breaks = rev(levels(df$Model))) +
    scale_x_discrete(name = "", limits = rev(levels(df$Coefficient))) +
    theme_bw()
  
  ggsave(file.name, units = "in", width = 7, height = 7)
}
