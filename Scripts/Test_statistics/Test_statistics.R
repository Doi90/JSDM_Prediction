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
    
    apk <- apk(actual = obs_data,
               predicted = pred_data,
               k = length(pred_data))
    
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
    
    confusion <- table(pred_data, obs_data)
    
    TP <- confusion[2,2]  # True positives
    FP <- confusion[2,1]  # False positives
    TN <- confusion[1,1]  # True negatives
    FN <- confusion[1,2]  # False negatives
    
    TPR <- TP / (TP + FN) # True Positive Rate / Sensitivity
    
    FPR <- FP / (FP + TN) # False Positive Rate
    
    TNR <- TN / (FP + TN) # True Negative Rate / Specificity
    
    FNR <- FN / (TP + FN) # False Negative Rate
    
    PLR <- TPR / FPR # Positive Likelihood Ratio
    
    NLR <- FNR / TNR # Negative Likelihood Ratio
    
    DOR <- PLR / NLR # Diagnositic Odds Ratio
    
    Prevalence <- (TP + FN) / (TP + FN + TN + FN)
    
    Accuracy <- (TP + TN) / (TP + FN + TN + FN)
    
    PPV <- TP / (TP + FP) # Positive Predictive Value
    
    FOR <- FN / (FN + TN) # False Omission Rate
    
    FDR <- FP / (TP + FP) # False Discovery Rate  
    
    NPV <- TN / (FN + TN) # Negative Predictive Value
    
    F_1 <- 2 / ((1 / TPR) + (1 / PPV)) # F1 Score. Equivalent to Sorenson-Dice?
    
    Manhattan <- vegdist(rbind(obs_data,
                               pred_data),
                         method = "manhattan",
                         binary = TRUE)
    
    Euclidean <- vegdist(rbind(obs_data,
                               pred_data),
                         method = "euclidean",
                         binary = TRUE)
    
    Canberra <- vegdist(rbind(obs_data,
                              pred_data),
                        method = "canberra",
                        binary = TRUE)
    
    Bray <- vegdist(rbind(obs_data,
                          pred_data),
                    method = "bray",
                    binary = TRUE)
    
    Kulczynski <- vegdist(rbind(obs_data,
                                pred_data),
                          method = "kulczynski",
                          binary = TRUE)
    
    Jaccard <- vegdist(rbind(obs_data,
                             pred_data),
                       method = "jaccard",
                       binary = TRUE)
    
    Gower <- vegdist(rbind(obs_data,
                           pred_data),
                     method = "gower", 
                     binary = TRUE)
    
    altGower <- vegdist(rbind(obs_data,
                              pred_data),
                        method = "altGower", 
                        binary = TRUE)
    
    Morisita <- vegdist(rbind(obs_data,
                              pred_data),
                        method = "morisita", 
                        binary = TRUE)
    
    Horn <- vegdist(rbind(obs_data,
                          pred_data),
                    method = "horn", 
                    binary = TRUE)
    
    Mountford <- vegdist(rbind(obs_data,
                               pred_data),
                         method = "mountford",
                         binary = TRUE)
    
    Raup <- vegdist(rbind(obs_data,
                          pred_data),
                    method = "raup", 
                    binary = TRUE)
    
    Binomial <- vegdist(rbind(obs_data,
                              pred_data),
                        method = "binomial", 
                        binary = TRUE)
    
    Chao <- vegdist(rbind(obs_data,
                          pred_data),
                    method = "chao", 
                    binary = TRUE)
    
    Cao <- vegdist(rbind(obs_data,
                         pred_data),
                   method = "cao", 
                   binary = TRUE)
    
    Mahalanobis <- vegdist(rbind(obs_data,
                                 pred_data),
                           method = "mahalanobis", 
                           binary = TRUE)
    
  }
}

