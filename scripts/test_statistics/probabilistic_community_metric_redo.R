#############################################################
#############################################################
###                                                       ###
###  COMMUNITY DISSIMILARITY ON PROBABLISTIC PREDICTIONS  ###
###                                                       ###
###       This script calculates the various test         ###
###    statistics used in this analysis for the JSDMs     ###
###                                                       ###
#############################################################
#############################################################

#####################
### Load Packages ###
#####################

library(TruncatedNormal)
library(Metrics)
library(caret)
library(vegan)
library(psych)
library(abind)

#######################
### Loop parameters ###
#######################

model_options <- c("MPR")#,
                   # "HPR",
                   # "LPR",
                   # "DPR",
                   # "HLR_NS",
                   # "HLR_S")

dataset_options <- c("frog")#,
                     # "eucalypt",
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

fold_options <- 1:5

###################################
### New test statistic function ###
###################################

## Normal function for everything except conditional prediction

test_statistic_replace <- function(observed = NULL,
                                   predictions = NULL){
  
  if(is.null(observed)){
    stop("observed not supplied")
  }
  
  if(is.null(predictions)){
    stop("observed not supplied")
  }
  
  ## Create empty arrays to store results ----
  
  ### Metrics calculated by site
  
  metric_names_site <- c(
    # Community dissimilarity metrics
    "Binomial", "Bray", "Canberra", "Euclidean", "Gower", "Gower_alt",
    "Horn", "Jaccard", "Kulczynski", "Mahalanobis", "Manhattan", "Mountford", "Raup"
  )
  
  test_statistics_site <- array(NA,
                                dim = c(nrow(observed),
                                        length(metric_names_site),
                                        dim(predictions)[3]),
                                dimnames = list(NULL,
                                                metric_names_site,
                                                NULL))
  
  ###----
  
  ## Calculate test statistics by site ----
  
  for(s in seq_len(dim(predictions)[3])){
    
    for(i in seq_len(dim(predictions)[1])){
      
      ## Extract observed and predicted data
      
      obs_data <- unlist(observed[i, ])
      
      pred_data <- predictions[i, , s]
      
      ## Implement a catch for handling the NA predictions
      
      if(anyNA(pred_data)){
        
        test_statistics_site[ i, , s] <- rep(NA, length(metric_names_site))
        
        next()
        
      }
      
      ## Calculate test statistics
      
      ### Community Dissimilarity ----
      
      #### Binary classification of predictions
      
      bin_pred <- rbinom(n = length(pred_data),
                            size = 1,
                            prob = pred_data)
      
      #### Calulate community dissimilarity indices
      
      Binomial <- tryCatch(expr = vegdist(rbind(obs_data,
                                                bin_pred),
                                          method = "binomial", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Bray <- tryCatch(expr = vegdist(rbind(obs_data,
                                            bin_pred),
                                      method = "bray", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Canberra <- tryCatch(expr = vegdist(rbind(obs_data,
                                                bin_pred),
                                          method = "canberra", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Euclidean <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "euclidean", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Gower <- tryCatch(expr = vegdist(rbind(obs_data,
                                             bin_pred),
                                       method = "gower", 
                                       binary = TRUE),
                        error = function(err){ return(NA) })
      
      Gower_alt <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "altGower", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Horn <- tryCatch(expr = vegdist(rbind(obs_data,
                                            bin_pred),
                                      method = "horn", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Jaccard <- tryCatch(expr = vegdist(rbind(obs_data,
                                               bin_pred),
                                         method = "jaccard", 
                                         binary = TRUE),
                          error = function(err){ return(NA) })
      
      Kulczynski <- tryCatch(expr = vegdist(rbind(obs_data,
                                                  bin_pred),
                                            method = "kulczynski", 
                                            binary = TRUE),
                             error = function(err){ return(NA) })
      
      Mahalanobis <- tryCatch(expr = vegdist(rbind(obs_data,
                                                   bin_pred),
                                             method = "mahalanobis", 
                                             binary = TRUE),
                              error = function(err){ return(NA) })
      
      Manhattan <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "manhattan", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Mountford <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "mountford", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Raup <- tryCatch(expr = vegdist(rbind(obs_data,
                                            bin_pred),
                                      method = "raup", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      ###----
      
      test_statistics_site[ i, , s] <- c(
        # Community dissimilarity metrics
        Binomial, Bray, Canberra, Euclidean, Gower, Gower_alt,
        Horn, Jaccard, Kulczynski, Mahalanobis, Manhattan, Mountford, Raup
      )
      
    }
  }
  
  ###----
  
  ts_out <- test_statistics_site
  
  return(ts_out)
  
}

## Function for conditional predictions only

test_statistic_conditional_replace <- function(observed = NULL,
                                               predictions = NULL,
                                               dataset_id = NULL,
                                               species_id = NULL){
  
  if(is.null(observed)){
    stop("observed not supplied")
  }
  
  if(is.null(predictions)){
    stop("observed not supplied")
  }
  
  if(is.null(dataset_id)){
    stop("dataset_id not supplied")
  }
  
  if(is.null(species_id)){
    stop("species_id not supplied")
  }
  
  ## Identify the species being left in
  
  if(dataset_id == "frog"){
    
    species_left_in_IDs <- c(9, 4, 6) # Lit_rani, Lim_tas, Lit_ewing
    
  }
  
  if(dataset_id == "eucalypt"){
    
    species_left_in_IDs <- c(8, 9, 3) # OVA, WIL, BAX
    
  }
  
  if(dataset_id == "bird"){
    
    species_left_in_IDs <- c(16, 19, 289) # Common_Goldeneye, Canada_Goose, American_Robin
    
  }
  
  if(dataset_id == "butterfly"){
    
    species_left_in_IDs <- c(5, 6, 43) # high_brown_fritillary, dark.green_fritillary, common_blue
    
  }
  
  if(dataset_id %in% c("sim1random",
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
                       "sim10spatial")){
    
    spp_LOI <- readRDS("data/simulated_datasets_species_left_in.rds")
    
    species_left_in_IDs <- as.numeric(spp_LOI[dataset_id, ])
    
  }
  
  species_left_in <- species_left_in_IDs[species_id]
  
  ## Create empty arrays to store results ----
  
  ### Metrics calculated by site
  
  metric_names_site <- c(
    # Community dissimilarity metrics
    "Binomial", "Bray", "Canberra", "Euclidean", "Gower", "Gower_alt",
    "Horn", "Jaccard", "Kulczynski", "Mahalanobis", "Manhattan", "Mountford", "Raup"
  )
  
  test_statistics_site <- array(NA,
                                dim = c(nrow(observed),
                                        length(metric_names_site),
                                        dim(predictions)[3]),
                                dimnames = list(NULL,
                                                metric_names_site,
                                                NULL))
  
  
  ###----
  
  ## Calculate test statistics by site ----
  
  for(s in seq_len(dim(predictions)[3])){
    
    for(i in seq_len(dim(predictions)[1])){
      
      ## Extract observed and predicted data
      
      obs_data <- unlist(observed[i, -species_left_in])
      
      pred_data <- predictions[i, -species_left_in, s]
      
      ## Calculate test statistics
      
      ### Community Dissimilarity ----
      
      #### Binary classification of predictions
      
      bin_pred <- rbinom(length(pred_data),
                         1,
                         pred_data)
      
      #### Calulate community dissimilarity indices
      
      Binomial <- tryCatch(expr = vegdist(rbind(obs_data,
                                                bin_pred),
                                          method = "binomial", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Bray <- tryCatch(expr = vegdist(rbind(obs_data,
                                            bin_pred),
                                      method = "bray", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Canberra <- tryCatch(expr = vegdist(rbind(obs_data,
                                                bin_pred),
                                          method = "canberra", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Euclidean <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "euclidean", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Gower <- tryCatch(expr = vegdist(rbind(obs_data,
                                             bin_pred),
                                       method = "gower", 
                                       binary = TRUE),
                        error = function(err){ return(NA) })
      
      Gower_alt <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "altGower", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Horn <- tryCatch(expr = vegdist(rbind(obs_data,
                                            bin_pred),
                                      method = "horn", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Jaccard <- tryCatch(expr = vegdist(rbind(obs_data,
                                               bin_pred),
                                         method = "jaccard", 
                                         binary = TRUE),
                          error = function(err){ return(NA) })
      
      Kulczynski <- tryCatch(expr = vegdist(rbind(obs_data,
                                                  bin_pred),
                                            method = "kulczynski", 
                                            binary = TRUE),
                             error = function(err){ return(NA) })
      
      Mahalanobis <- tryCatch(expr = vegdist(rbind(obs_data,
                                                   bin_pred),
                                             method = "mahalanobis", 
                                             binary = TRUE),
                              error = function(err){ return(NA) })
      
      Manhattan <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "manhattan", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Mountford <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 bin_pred),
                                           method = "mountford", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Raup <- tryCatch(expr = vegdist(rbind(obs_data,
                                            bin_pred),
                                      method = "raup", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      ###----
      
      test_statistics_site[ i, , s] <- c(
        # Community dissimilarity metrics
        Binomial, Bray, Canberra, Euclidean, Gower, Gower_alt,
        Horn, Jaccard, Kulczynski, Mahalanobis, Manhattan, Mountford, Raup
      )
      
    }
  }
  
  ###----
  
  ts_out <- test_statistics_site
  
  return(ts_out)
  
}

############################
### Loop over iterations ###
############################

for(model_id in model_options){
  
  for(dataset_id in dataset_options){
    
    for(fold_id in fold_options){
      
      #################
      ### Load Data ###
      #################
      
      ## Observed / Testing Data
      
      command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test.csv')", 
                         dataset_id,                          # Need to build command to read in
                         fold_id)                             # specific files for this CV fold
      
      y_test <- eval(parse(text = command))                   # Evaluate command to read in data
      
      y_test <- y_test[ , -1]                                 # Remove rownames
      
      ## Marginal - probabilities ----
      
      ### Load Data
      
      filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_prob.rds",
                          model_id,
                          dataset_id,
                          fold_id)
      
      marg_pred_prob <- readRDS(filename)
      
      filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_marginal_prob_ts.rds",
                          model_id,
                          dataset_id,
                          fold_id)
      
      marg_ts_prob <- readRDS(filename)
      
      ### Calculate Test Statistics
      
      marg_ts_prob$test_statistics_site <- test_statistic_replace(observed = y_test,
                                                                  predictions = marg_pred_prob)
      
      ### Save To File
      
      filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_marginal_prob_ts.rds",
                          model_id,
                          dataset_id,
                          fold_id)
      
      saveRDS(marg_ts_prob,
              filename)
      
      ### Memory purge
      
      rm(marg_pred_prob,
         marg_ts_prob)
      
      ## Conditional Marginal - Leave One In ----
      
      ### Load Data
      
      filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI_marg.rds",
                          model_id,
                          dataset_id,
                          fold_id)
      
      cond_marg_LOI_pred <- readRDS(filename)
      
      filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_condLOI_marg_ts.rds",
                          model_id,
                          dataset_id,
                          fold_id)
      
      cond_marg_LOI_ts <- readRDS(filename)
      
      ### Calculate Test Statistics
      
      # Loop over array 4th dimension, calculate test statistic separately
      # Save as list, turn back into 4D array
      
      ts_list <- vector(mode = "list",
                        length = dim(cond_marg_LOI_pred)[4])
      
      for(i in seq_len(dim(cond_marg_LOI_pred)[4])){
        
        ts_list[[i]] <- test_statistic_conditional_replace(observed = y_test,
                                                           predictions = cond_marg_LOI_pred[ , , , i],
                                                           dataset_id = dataset_id,
                                                           species_id = i)
        
      }
      
      for(i in seq_len(dim(cond_marg_LOI_pred)[4])){
        
        cond_marg_LOI_ts[[i]]$test_statistics_site <- ts_list[[i]]
        
      }
      
      ### Save To File
      
      filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_condLOI_marg_ts.rds",
                          model_id,
                          dataset_id,
                          fold_id)
      
      saveRDS(cond_marg_LOI_ts,
              filename)
      
      ### Memory purge
      
      rm(cond_marg_LOI_pred,
         cond_marg_LOI_ts)
      
    }
  }
}
