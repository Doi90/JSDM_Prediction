######################################################
######################################################
###                                                ###
###           Calculate Test Statistics            ###
###                                                ###
###    This script calculates the various test     ###
###       statistics used in this analysis.        ###
###                                                ###
######################################################
######################################################

#####################
### Load Packages ###
#####################

library(Metrics)
library(caret)
library(vegan)

#################
### Load Data ###
#################

## Observed / Testing Data ----

# Matrix with <sites> rows and <species> columns

## Predicted Data ----

# Array with <sites> rows, <species> columns, and <samples> slices

#----

#################################
### Calculate Test Statistics ###
#################################

## Create empty array to store results

test_statistics <- array(NA,
                         dim = c(ncol(y),
                                 7,
                                 dim(predictions)[3]))

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
    
    RMSE <- rmse(actual = obs_data,
                 predicted = pred_data)
    
    SSE <- sse(actual = obs_data,
               predicted = pred_data)
    
    ###----
    
    ## Set a threshold for binary classification for remaining metrics
    
    threshold <- 0.5
    
    ### Threshold dependent ----
    
    #### Create a confusion matrix and extract values
    
    confusion <- table(pred_data > threshold, obs_data)
    
    TP <- confusion[2,2]  # True positives
    FP <- confusion[2,1]  # False positives
    TN <- confusion[1,1]  # True negatives
    FN <- confusion[1,2]  # False negatives
    
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
    
    ###----
    
    ### Community Dissimilarity ----
    
    ###----
    
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
 
  }
}

