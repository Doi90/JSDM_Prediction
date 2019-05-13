######################################################
######################################################
###                                                ###
###           Test Statistics Function             ###
###                                                ###
###    This script defines the various test        ###
###       statistic function.                      ###
###                                                ###
######################################################
######################################################

## Normal function for everything except conditional prediction

test_statistic <- function(observed = NULL,
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
  
  ### Metrics calculated by species
  
  metric_names_species <- c(
    # Threshold-independent
    "AUC", "bias", "MSE", "R2", "RMSE", "SSE", "Pearson", "Spearman", "Kendall",
    # Threshold-dependant
    "TP", "FP", "TN", "FN",
    "TPR", "FPR", "TNR", "FNR", "PLR", "NLR", "DOR", "Prevalence", "Accuracy",
    "PPV", "FOR", "FDR", "NPV", "F_1", "Youden_J", "Kappa"
  )
  
  test_statistics_species <- array(NA,
                                   dim = c(ncol(observed),
                                           length(metric_names_species),
                                           dim(predictions)[3]),
                                   dimnames = list(NULL,
                                                   metric_names_species,
                                                   NULL))
  
  ###----
  
  ## Calculate test statistics by species ----
  
  for(s in seq_len(dim(predictions)[3])){
    
    for(j in seq_len(dim(predictions)[2])){
      
      ## Extract observed and predicted data
      
      iteration_data <- cbind(unlist(observed[ , j]),
                              predictions[ , j, s])
      
      iteration_data <- iteration_data[!apply(iteration_data, 1, anyNA), , drop = FALSE]
      
      obs_data <- iteration_data[ , 1]
      
      pred_data <- iteration_data[ , 2]
      
      # ## Implement a catch for handling the NA predictions
      # 
      # if(anyNA(pred_data)){
      #   
      #   test_statistics[ , j, s] <- rep(NA, 42)
      #   
      #   NA_count <- NA_count + 1
      #   
      #   next()
      #   
      # }
      
      ## Calculate test statistics
      
      ### Threshold independent ----
      
      AUC <- tryCatch(expr = auc(actual = obs_data,
                                 predicted = pred_data),
                      error = function(err){ return(NA) })
      
      bias <- tryCatch(expr = bias(actual = obs_data,
                                   predicted = pred_data),
                       error = function(err){ return(NA) })
      
      MSE <- tryCatch(expr = mse(actual = obs_data,
                                 predicted = pred_data),
                      error = function(err){ return(NA) })
      
      R2 <- tryCatch(expr = R2(obs = obs_data,
                               pred = pred_data),
                     error = function(err){ return(NA) })
      
      RMSE <- tryCatch(expr = rmse(actual = obs_data,
                                   predicted = pred_data),
                       error = function(err){ return(NA) })
      
      SSE <- tryCatch(expr = sse(actual = obs_data,
                                 predicted = pred_data),
                      error = function(err){ return(NA) })
      
      Pearson <- tryCatch(expr = cor(obs_data,
                                     pred_data,
                                     method = "pearson"),
                          error = function(err){ return(NA) })
      
      Spearman <- tryCatch(expr = cor(obs_data,
                                      pred_data,
                                      method = "spearman"),
                           error = function(err){ return(NA) })
      
      Kendall <- tryCatch(expr = cor(obs_data,
                                     pred_data,
                                     method = "kendall"),
                          error = function(err){ return(NA) })
      
      ###----
      
      ## Set a threshold for binary classification for remaining metrics
      
      threshold <- 0.5
      
      ### Threshold dependent ----
      
      #### Create a binary confusion matrix and extract values
      
      # confusion <- table(pred_data > threshold, obs_data)
      # 
      # TP <- confusion[2,2]  # True positives
      # FP <- confusion[2,1]  # False positives
      # TN <- confusion[1,1]  # True negatives
      # FN <- confusion[1,2]  # False negatives
      
      #### Create a probabilistic confusion matrix and extract values
      
      TP <- sum(pred_data * obs_data)
      FP <- sum(pred_data * (1 - obs_data))
      TN <- sum((1 - pred_data) * (1 - obs_data))
      FN <- sum((1 - pred_data) * obs_data)
      
      #### Calculate metrics
      
      TPR <- tryCatch(expr = TP / (TP + FN),
                      error = function(err){ return(NA) })                          # True Positive Rate / Sensitivity
      
      FPR <- tryCatch(expr = FP / (FP + TN),
                      error = function(err){ return(NA) })                         # False Positive Rate
      
      TNR <- tryCatch(expr = TN / (FP + TN),
                      error = function(err){ return(NA) })                         # True Negative Rate / Specificity
      
      FNR <- tryCatch(expr = FN / (TP + FN),
                      error = function(err){ return(NA) })                        # False Negative Rate
      
      PLR <- tryCatch(expr = TPR / FPR,
                      error = function(err){ return(NA) })                               # Positive Likelihood Ratio
      
      NLR <- tryCatch(expr = FNR / TNR,
                      error = function(err){ return(NA) })                               # Negative Likelihood Ratio
      
      DOR <- tryCatch(expr = PLR / NLR,
                      error = function(err){ return(NA) })                               # Diagnositic Odds Ratio
      
      Prevalence <- tryCatch(expr = (TP + FN) / (TP + FP + TN + FN),
                             error = function(err){ return(NA) })
      
      Accuracy <- tryCatch(expr = (TP + TN) / (TP + FP + TN + FN),
                           error = function(err){ return(NA) })
      
      PPV <- tryCatch(expr = TP / (TP + FP),
                      error = function(err){ return(NA) })                          # Positive Predictive Value
      
      FOR <- tryCatch(expr = FN / (FN + TN),
                      error = function(err){ return(NA) })                          # False Omission Rate
      
      FDR <- tryCatch(expr = FP / (TP + FP),
                      error = function(err){ return(NA) })                          # False Discovery Rate  
      
      NPV <- tryCatch(expr = TN / (FN + TN),
                      error = function(err){ return(NA) })                          # Negative Predictive Value
      
      F_1 <- tryCatch(expr = 2 / ((1 / TPR) + (1 / PPV)),
                      error = function(err){ return(NA) })             # F1 Score. Equivalent to Sorenson-Dice?
      
      Youden_J <- tryCatch(expr = TPR + TNR - 1,
                           error = function(err){ return(NA) })                      # Youden's J statistic
      
      Kappa <- tryCatch(expr = cohen.kappa(cbind(obs_data,           # Cohen's Kappa
                                                 pred_data))$kappa,
                        error = function(err){ return(NA) })
      
      # Kappa <- ((TP + TN) - 
      #             (((TP + FN) * (TP + FP) + (FP + TN) * (FN + TN)) /
      #                (TP + FP + TN + FN))) /
      #          ((TP + FP + TN + FN) -
      #             (((TP + FN) * (TP + FP) + (FP + TN) * (FN + TN)) /
      #                (TP + FP + TN + FN))) 
      
      ###----
      
      test_statistics_species[ j, , s] <- c(
        # Threshold-independent
        AUC, bias, MSE, R2, RMSE, SSE, Pearson, Spearman, Kendall,
        # Threshold-dependant
        TP, FP, TN, FN, TPR, FPR, TNR, FNR, PLR, NLR, DOR, Prevalence,
        Accuracy, PPV, FOR, FDR, NPV, F_1, Youden_J, Kappa
      )
      
    }
  }
  
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
      
      ### Set a threshold for binary classification
      
      threshold <- 0.5
      
      ### Community Dissimilarity ----
      
      #### Binary classification of predictions
      
      pred_thresh <- as.numeric(pred_data > threshold)
      
      #### Calulate community dissimilarity indices
      
      Binomial <- tryCatch(expr = vegdist(rbind(obs_data,
                                                pred_thresh),
                                          method = "binomial", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Bray <- tryCatch(expr = vegdist(rbind(obs_data,
                                            pred_thresh),
                                      method = "bray", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Canberra <- tryCatch(expr = vegdist(rbind(obs_data,
                                                pred_thresh),
                                          method = "canberra", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Euclidean <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "euclidean", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Gower <- tryCatch(expr = vegdist(rbind(obs_data,
                                             pred_thresh),
                                       method = "gower", 
                                       binary = TRUE),
                        error = function(err){ return(NA) })
      
      Gower_alt <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "altGower", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Horn <- tryCatch(expr = vegdist(rbind(obs_data,
                                            pred_thresh),
                                      method = "horn", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Jaccard <- tryCatch(expr = vegdist(rbind(obs_data,
                                               pred_thresh),
                                         method = "jaccard", 
                                         binary = TRUE),
                          error = function(err){ return(NA) })
      
      Kulczynski <- tryCatch(expr = vegdist(rbind(obs_data,
                                                  pred_thresh),
                                            method = "kulczynski", 
                                            binary = TRUE),
                             error = function(err){ return(NA) })
      
      Mahalanobis <- tryCatch(expr = vegdist(rbind(obs_data,
                                                   pred_thresh),
                                             method = "mahalanobis", 
                                             binary = TRUE),
                              error = function(err){ return(NA) })
      
      Manhattan <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "manhattan", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Mountford <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "mountford", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Raup <- tryCatch(expr = vegdist(rbind(obs_data,
                                            pred_thresh),
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
  
  ts_out <- list(test_statistics_species = test_statistics_species,
                 test_statistics_site = test_statistics_site)
  
  return(ts_out)
  
}

## Function for conditional predictions only

test_statistic_conditional <- function(observed = NULL,
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
  
  ### Metrics calculated by species
  
  metric_names_species <- c(
    # Threshold-independent
    "AUC", "bias", "MSE", "R2", "RMSE", "SSE", "Pearson", "Spearman", "Kendall",
    # Threshold-dependant
    "TP", "FP", "TN", "FN",
    "TPR", "FPR", "TNR", "FNR", "PLR", "NLR", "DOR", "Prevalence", "Accuracy",
    "PPV", "FOR", "FDR", "NPV", "F_1", "Youden_J", "Kappa"
  )
  
  test_statistics_species <- array(NA,
                                   dim = c(ncol(observed),
                                           length(metric_names_species),
                                           dim(predictions)[3]),
                                   dimnames = list(NULL,
                                                   metric_names_species,
                                                   NULL))
  
  ###----
  
  ## Calculate test statistics by species ----
  
  for(s in seq_len(dim(predictions)[3])){
    
    for(j in seq_len(dim(predictions)[2])){
      
      ## Don't calculate test statistics for the known species (not predictions)
      
      if(j == species_left_in){
        
        test_statistics_species[ j, , s] <- NA
        
        next()
        
      }
      
      ## Extract observed and predicted data
      
      iteration_data <- cbind(unlist(observed[ , j]),
                              predictions[ , j, s])
      
      iteration_data <- iteration_data[!apply(iteration_data, 1, anyNA), , drop = FALSE]
      
      obs_data <- iteration_data[ , 1]
      
      pred_data <- iteration_data[ , 2]
      
      # ## Implement a catch for handling the NA predictions
      # 
      # if(anyNA(pred_data)){
      #   
      #   test_statistics[ , j, s] <- rep(NA, 42)
      #   
      #   NA_count <- NA_count + 1
      #   
      #   next()
      #   
      # }
      
      ## Calculate test statistics
      
      ### Threshold independent ----
      
      AUC <- tryCatch(expr = auc(actual = obs_data,
                                 predicted = pred_data),
                      error = function(err){ return(NA) })
      
      bias <- tryCatch(expr = bias(actual = obs_data,
                                   predicted = pred_data),
                       error = function(err){ return(NA) })
      
      MSE <- tryCatch(expr = mse(actual = obs_data,
                                 predicted = pred_data),
                      error = function(err){ return(NA) })
      
      R2 <- tryCatch(expr = R2(obs = obs_data,
                               pred = pred_data),
                     error = function(err){ return(NA) })
      
      RMSE <- tryCatch(expr = rmse(actual = obs_data,
                                   predicted = pred_data),
                       error = function(err){ return(NA) })
      
      SSE <- tryCatch(expr = sse(actual = obs_data,
                                 predicted = pred_data),
                      error = function(err){ return(NA) })
      
      Pearson <- tryCatch(expr = cor(obs_data,
                                     pred_data,
                                     method = "pearson"),
                          error = function(err){ return(NA) })
      
      Spearman <- tryCatch(expr = cor(obs_data,
                                      pred_data,
                                      method = "spearman"),
                           error = function(err){ return(NA) })
      
      Kendall <- tryCatch(expr = cor(obs_data,
                                     pred_data,
                                     method = "kendall"),
                          error = function(err){ return(NA) })
      
      ###----
      
      ## Set a threshold for binary classification for remaining metrics
      
      threshold <- 0.5
      
      ### Threshold dependent ----
      
      #### Create a binary confusion matrix and extract values
      
      # confusion <- table(pred_data > threshold, obs_data)
      # 
      # TP <- confusion[2,2]  # True positives
      # FP <- confusion[2,1]  # False positives
      # TN <- confusion[1,1]  # True negatives
      # FN <- confusion[1,2]  # False negatives
      
      #### Create a probabilistic confusion matrix and extract values
      
      TP <- sum(pred_data * obs_data)
      FP <- sum(pred_data * (1 - obs_data))
      TN <- sum((1 - pred_data) * (1 - obs_data))
      FN <- sum((1 - pred_data) * obs_data)
      
      #### Calculate metrics
      
      TPR <- tryCatch(expr = TP / (TP + FN),
                      error = function(err){ return(NA) })                          # True Positive Rate / Sensitivity
      
      FPR <- tryCatch(expr = FP / (FP + TN),
                      error = function(err){ return(NA) })                         # False Positive Rate
      
      TNR <- tryCatch(expr = TN / (FP + TN),
                      error = function(err){ return(NA) })                         # True Negative Rate / Specificity
      
      FNR <- tryCatch(expr = FN / (TP + FN),
                      error = function(err){ return(NA) })                        # False Negative Rate
      
      PLR <- tryCatch(expr = TPR / FPR,
                      error = function(err){ return(NA) })                               # Positive Likelihood Ratio
      
      NLR <- tryCatch(expr = FNR / TNR,
                      error = function(err){ return(NA) })                               # Negative Likelihood Ratio
      
      DOR <- tryCatch(expr = PLR / NLR,
                      error = function(err){ return(NA) })                               # Diagnositic Odds Ratio
      
      Prevalence <- tryCatch(expr = (TP + FN) / (TP + FP + TN + FN),
                             error = function(err){ return(NA) })
      
      Accuracy <- tryCatch(expr = (TP + TN) / (TP + FP + TN + FN),
                           error = function(err){ return(NA) })
      
      PPV <- tryCatch(expr = TP / (TP + FP),
                      error = function(err){ return(NA) })                          # Positive Predictive Value
      
      FOR <- tryCatch(expr = FN / (FN + TN),
                      error = function(err){ return(NA) })                          # False Omission Rate
      
      FDR <- tryCatch(expr = FP / (TP + FP),
                      error = function(err){ return(NA) })                          # False Discovery Rate  
      
      NPV <- tryCatch(expr = TN / (FN + TN),
                      error = function(err){ return(NA) })                          # Negative Predictive Value
      
      F_1 <- tryCatch(expr = 2 / ((1 / TPR) + (1 / PPV)),
                      error = function(err){ return(NA) })             # F1 Score. Equivalent to Sorenson-Dice?
      
      Youden_J <- tryCatch(expr = TPR + TNR - 1,
                           error = function(err){ return(NA) })                      # Youden's J statistic
      
      Kappa <- tryCatch(expr = cohen.kappa(cbind(obs_data,           # Cohen's Kappa
                                                 pred_data))$kappa,
                        error = function(err){ return(NA) })
      
      # Kappa <- ((TP + TN) - 
      #             (((TP + FN) * (TP + FP) + (FP + TN) * (FN + TN)) /
      #                (TP + FP + TN + FN))) /
      #          ((TP + FP + TN + FN) -
      #             (((TP + FN) * (TP + FP) + (FP + TN) * (FN + TN)) /
      #                (TP + FP + TN + FN))) 
      
      ###----
      
      test_statistics_species[ j, , s] <- c(
        # Threshold-independent
        AUC, bias, MSE, R2, RMSE, SSE, Pearson, Spearman, Kendall,
        # Threshold-dependant
        TP, FP, TN, FN, TPR, FPR, TNR, FNR, PLR, NLR, DOR, Prevalence,
        Accuracy, PPV, FOR, FDR, NPV, F_1, Youden_J, Kappa
      )
      
    }
  }
  
  ###----
  
  ## Calculate test statistics by site ----
  
  for(s in seq_len(dim(predictions)[3])){
    
    for(i in seq_len(dim(predictions)[1])){
      
      ## Extract observed and predicted data
      
      obs_data <- unlist(observed[i, -species_left_in])
      
      pred_data <- predictions[i, -species_left_in, s]
      
      ## Implement a catch for handling the NA predictions
      
      if(anyNA(pred_data)){
        
        test_statistics_site[ i, , s] <- rep(NA, length(metric_names_site))
        
        next()
        
      }
      
      ## Calculate test statistics
      
      ### Set a threshold for binary classification
      
      threshold <- 0.5
      
      ### Community Dissimilarity ----
      
      #### Binary classification of predictions
      
      pred_thresh <- as.numeric(pred_data > threshold)
      
      #### Calulate community dissimilarity indices
      
      Binomial <- tryCatch(expr = vegdist(rbind(obs_data,
                                                pred_thresh),
                                          method = "binomial", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Bray <- tryCatch(expr = vegdist(rbind(obs_data,
                                            pred_thresh),
                                      method = "bray", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Canberra <- tryCatch(expr = vegdist(rbind(obs_data,
                                                pred_thresh),
                                          method = "canberra", 
                                          binary = TRUE),
                           error = function(err){ return(NA) })
      
      Euclidean <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "euclidean", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Gower <- tryCatch(expr = vegdist(rbind(obs_data,
                                             pred_thresh),
                                       method = "gower", 
                                       binary = TRUE),
                        error = function(err){ return(NA) })
      
      Gower_alt <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "altGower", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Horn <- tryCatch(expr = vegdist(rbind(obs_data,
                                            pred_thresh),
                                      method = "horn", 
                                      binary = TRUE),
                       error = function(err){ return(NA) })
      
      Jaccard <- tryCatch(expr = vegdist(rbind(obs_data,
                                               pred_thresh),
                                         method = "jaccard", 
                                         binary = TRUE),
                          error = function(err){ return(NA) })
      
      Kulczynski <- tryCatch(expr = vegdist(rbind(obs_data,
                                                  pred_thresh),
                                            method = "kulczynski", 
                                            binary = TRUE),
                             error = function(err){ return(NA) })
      
      Mahalanobis <- tryCatch(expr = vegdist(rbind(obs_data,
                                                   pred_thresh),
                                             method = "mahalanobis", 
                                             binary = TRUE),
                              error = function(err){ return(NA) })
      
      Manhattan <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "manhattan", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Mountford <- tryCatch(expr = vegdist(rbind(obs_data,
                                                 pred_thresh),
                                           method = "mountford", 
                                           binary = TRUE),
                            error = function(err){ return(NA) })
      
      Raup <- tryCatch(expr = vegdist(rbind(obs_data,
                                            pred_thresh),
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
  
  ts_out <- list(test_statistics_species = test_statistics_species,
                 test_statistics_site = test_statistics_site)
  
  return(ts_out)
  
}
