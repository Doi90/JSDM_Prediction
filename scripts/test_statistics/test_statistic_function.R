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

test_statistic <- function(observed,
                           predictions){

  ## Create empty array to store results
  
  metric_names <- c(
    # Threshold-independent
    "AUC", "bias", "MSE", "R2", "RMSE", "SSE", "Pearson", "Spearman", "Kendall",
    # Threshold-dependant
    "TP", "FP", "TN", "FN",
    "TPR", "FPR", "TNR", "FNR", "PLR", "NLR", "DOR", "Prevalence", "Accuracy",
    "PPV", "FOR", "FDR", "NPV", "F_1", "Youden_J", "Kappa",
    # Community dissimilarity metrics
    "Binomial", "Bray", "Canberra", "Cao", "Chao", "Euclidean", "Gower", "Gower_alt",
    "Horn", "Jaccard", "Kulczynski", "Mahalanobis", "Manhattan", "Morisita", "Mountford", "Raup"
  )
  
  test_statistics <- array(NA,
                           dim = c(nrow(observed),
                                   length(metric_names),
                                   dim(predictions)[3]),
                           dimnames = list(NULL,
                                           metric_names,
                                           NULL))
  
  ## Calculate test statistics
  
  for(a in seq_len(dim(predictions)[3])){
    
    for(i in seq_len(dim(predictions)[1])){
      
      ## Extract observed and predicted data
      
      obs_data <- observed[i, ]
      
      pred_data <- predictions[i, , a]
      
      ## Calculate test statistics
      
      ### Threshold independent ----
      
      AUC <- auc(actual = obs_data,
                 predicted = pred_data)
      
      bias <- bias(actual = obs_data,
                   predicted = pred_data)
      
      MSE <- mse(actual = obs_data,
                 predicted = pred_data)
      
      R2 <- R2(obs = obs_data,
               pred = pred_data)
      
      RMSE <- rmse(actual = obs_data,
                   predicted = pred_data)
      
      SSE <- sse(actual = obs_data,
                 predicted = pred_data)
      
      Pearson <- cor(obs_data,
                     pred_data,
                     method = "pearson")
      
      Spearman <- cor(obs_data,
                      pred_data,
                      method = "spearman")
      
      Kendall <- cor(obs_data,
                     pred_data,
                     method = "kendall")
      
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
      TN <- sum(1 - pred_data) * (1 - obs_data)
      FN <- sum((1 - pred_data) * obs_data)
      
      #### Calculate metrics
      
      TPR <- TP / (TP + FN)                          # True Positive Rate / Sensitivity
      FPR <- FP / (FP + TN)                          # False Positive Rate
      TNR <- TN / (FP + TN)                          # True Negative Rate / Specificity
      FNR <- FN / (TP + FN)                          # False Negative Rate
      PLR <- TPR / FPR                               # Positive Likelihood Ratio
      NLR <- FNR / TNR                               # Negative Likelihood Ratio
      DOR <- PLR / NLR                               # Diagnositic Odds Ratio
      Prevalence <- (TP + FN) / (TP + FN + TN + FN)
      Accuracy <- (TP + TN) / (TP + FN + TN + FN)
      PPV <- TP / (TP + FP)                          # Positive Predictive Value
      FOR <- FN / (FN + TN)                          # False Omission Rate
      FDR <- FP / (TP + FP)                          # False Discovery Rate  
      NPV <- TN / (FN + TN)                          # Negative Predictive Value
      F_1 <- 2 / ((1 / TPR) + (1 / PPV))             # F1 Score. Equivalent to Sorenson-Dice?
      Youden_J <- TPR + TNR -1                       # Youden's J statistic
      Kappa <- cohen.kappa(cbind(obs_data,           # Cohen's Kappa
                                 pred_data))$kappa
      # Kappa <- ((TP + TN) - 
      #             (((TP + FN) * (TP + FP) + (FP + TN) + (FN + TN)) /
      #                (TP + FP + TN + FN))) /
      #          ((TP + FP + TN + FN) -
      #             (((TP + FN) * (TP + FP) + (FP + TN) * (FN + TN)) /
      #                (TP + FP + TN + FN))) 
      
      ###----
      
      ### Community Dissimilarity ----
      
      #### Binary classification of predictions
      
      pred_thresh <- as.numeric(pred_data > threshold)
      
      #### Calulate community dissimilarity indices
      
      Binomial <- vegdist(rbind(obs_data,
                                pred_thresh),
                          method = "binomial", 
                          binary = TRUE)
      
      Bray <- vegdist(rbind(obs_data,
                            pred_thresh),
                      method = "bray",
                      binary = TRUE)
      
      Canberra <- vegdist(rbind(obs_data,
                                pred_thresh),
                          method = "canberra",
                          binary = TRUE)
      
      Cao <- vegdist(rbind(obs_data,
                           pred_thresh),
                     method = "cao", 
                     binary = TRUE)
      
      Chao <- vegdist(rbind(obs_data,
                            pred_thresh),
                      method = "chao", 
                      binary = TRUE)
      
      Euclidean <- vegdist(rbind(obs_data,
                                 pred_thresh),
                           method = "euclidean",
                           binary = TRUE)
      
      Gower <- vegdist(rbind(obs_data,
                             pred_thresh),
                       method = "gower", 
                       binary = TRUE)
      
      Gower_alt <- vegdist(rbind(obs_data,
                                 pred_thresh),
                           method = "altGower", 
                           binary = TRUE)
      
      Horn <- vegdist(rbind(obs_data,
                            pred_thresh),
                      method = "horn", 
                      binary = TRUE)
      
      Jaccard <- vegdist(rbind(obs_data,
                               pred_thresh),
                         method = "jaccard",
                         binary = TRUE)
      
      Kulczynski <- vegdist(rbind(obs_data,
                                  pred_thresh),
                            method = "kulczynski",
                            binary = TRUE)
      
      Mahalanobis <- vegdist(rbind(obs_data,
                                   pred_thresh),
                             method = "mahalanobis", 
                             binary = TRUE)
      
      Manhattan <- vegdist(rbind(obs_data,
                                 pred_thresh),
                           method = "manhattan",
                           binary = TRUE)
      
      Morisita <- vegdist(rbind(obs_data,
                                pred_thresh),
                          method = "morisita", 
                          binary = TRUE)
      
      Mountford <- vegdist(rbind(obs_data,
                                 pred_thresh),
                           method = "mountford",
                           binary = TRUE)
      
      Raup <- vegdist(rbind(obs_data,
                            pred_thresh),
                      method = "raup", 
                      binary = TRUE)
      
      ###----
      
      test_statistics[i, , a] <- c(
        # Threshold-independent
        AUC, bias, MSE, R2, RMSE, SSE, Pearson, Spearman, Kendall,
        # Threshold-dependant
        TP, FP, TN, FN, TPR, FPR, TNR, FNR, PLR, NLR, DOR, Prevalence,
        Accuracy, PPV, FOR, FDR, NPV, F_1, Youden_J, Kappa,
        # Community dissimilarity metrics
        Binomial, Bray, Canberra, Cao, Chao, Euclidean, Gower, Gower_alt,
        Horn, Jaccard, Kulczynski, Mahalanobis, Manhattan, Morisita, Mountford, Raup
      )
      
    }
  }
  return(test_statistics)
}
