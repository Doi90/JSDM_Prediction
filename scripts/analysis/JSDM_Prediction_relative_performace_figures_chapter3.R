#####################
### Load Packages ###
#####################

library(nlme)
library(ggplot2)
library(RColorBrewer)
library(PassButter)
library(psych)
library(logitnorm)
library(MASS)

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

dataset_options <- c("frog",
                     "eucalypt",
                     #"bird",
                     #"butterfly",
                     "sim1random",
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
                             "condLOI_marg",
                             "SSDM_prob")


## Define prediction combinations for subsetting

pred_sets <- list(marginal_bin = c("marginal_bin", "SSDM_bin", "SESAM"),
                  marginal_prob = c("marginal_prob", "SSDM_prob"),
                  condLOI_low = c("condLOI_low", "SSDM_bin", "SESAM"),
                  condLOI_med = c("condLOI_med", "SSDM_bin", "SESAM"),
                  condLOI_high = c("condLOI_high", "SSDM_bin", "SESAM"),
                  condLOI_marg_low = c("condLOI_marg_low", "condLOI_marg", "SSDM_prob"),
                  condLOI_marg_med = c("condLOI_marg_med", "condLOI_marg", "SSDM_prob"),
                  condLOI_marg_high = c("condLOI_marg_high", "condLOI_marg", "SSDM_prob"),
                  joint = c("joint", "SSDM_bin", "SESAM"))

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
                #"TP",
                #"FP",
                #"TN",
                #"FN",
                "TPR",
                "FPR",
                "TNR",
                "FNR",
                "PLR",
                "NLR",
                "DOR",
                #"Prevalence",
                "Accuracy",
                "PPV",
                "FOR",
                "FDR",
                "NPV",
                "F_1",
                "Youden_J",
                "Kappa")

### By Site

ts_site <- c(#"Binomial",
  "Bray",
  "Canberra",
  #"Euclidean",
  "Gower",
  "Gower_alt",
  #"Horn",
  "Jaccard",
  "Kulczynski",
  #"Mahalanobis",
  #"Manhattan",
  "Mountford",
  "Raup")

## Test statistic / prediction type compatibility

binary_ts <- c(#"TP",
  #"FP",
  #"TN",
  #"FN",
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
  #"Binomial",
  "Bray",
  "Canberra",
  #"Euclidean",
  "Gower",
  "Gower_alt",
  #"Horn",
  "Jaccard",
  "Kulczynski",
  #"Mahalanobis",
  #"Manhattan",
  "Mountford",
  "Raup",
  "species_richness_difference")

prob_ts <- c("AUC",
             "bias",
             "MSE",
             "R2",
             "RMSE",
             "SSE",
             "Pearson",
             "Spearman",
             "Kendall",
             #"TP",
             #"FP",
             #"TN",
             #"FN",
             "TPR",
             "FPR",
             "TNR",
             "FNR",
             "PLR",
             "NLR",
             "DOR",
             #"Prevalence",
             "Accuracy",
             "PPV",
             "FOR",
             "FDR",
             "NPV",
             "F_1",
             "Youden_J",
             "Kappa",
             #"Binomial",
             "Bray",
             "Canberra",
             #"Euclidean",
             "Gower",
             "Gower_alt",
             #"Horn",
             "Jaccard",
             "Kulczynski",
             #"Mahalanobis",
             #"Manhattan",
             "Mountford",
             "Raup",
             "independent_log_likelihood",
             "joint_log_likelihood",
             "species_richness_difference")

## Test Statistic range

ts_Inf_Inf <- c("bias",
                "Pearson",
                "Spearman",
                "Kendall",
                "species_richness_difference",
                "independent_log_likelihood",
                "joint_log_likelihood")

ts_0_Inf <- c("MSE",
              "RMSE",
              "SSE",
              "TP",
              "FP",
              "TN",
              "FN",
              "PLR",
              "NLR",
              "DOR",
              "Binomial",
              "Euclidean",
              "Mahalanobis",
              "Manhattan")

ts_0_1 <- c("AUC",
            "R2",
            "TPR",
            "FPR",
            "TNR",
            "FNR",
            "Prevalence",
            "Accuracy",
            "PPV",
            "FOR",
            "FDR",
            "NPV",
            "F_1",
            "Youden_J",
            "Kappa",
            "Bray",
            "Canberra",
            "Gower",
            "Gower_alt",
            "Horn",
            "Jaccard",
            "Kulczynski",
            "Mountford",
            "Raup")

## Test statistic subsets

subsets <- list(ts_0_1_H = c("Accuracy",
                             "AUC",
                             "Kappa",
                             "F_1",
                             "PPV",
                             "NPV",
                             "R2",
                             "TPR",
                             "TNR",
                             "Youden_J"),
                ts_0_1_L_a = c("FDR",
                               "FNR",
                               "FOR",
                               "FPR"),
                ts_0_1_L_b = c("Bray",
                               "Canberra",
                               "Gower",
                               "Gower_alt",
                               "Jaccard",
                               "Kulczynski",
                               "Mountford",
                               "Raup"),
                ts_1_1 = c("Kendall",
                           "Pearson",
                           "Spearman"),
                ts_Inf_Inf_H = c("independent_log_likelihood",
                                 "joint_log_likelihood"),
                ts_Inf_Inf_0 = c("bias",
                                 "species_richness_difference"),
                ts_0_Inf_L = c("MSE",
                               "RMSE",
                               "SSE",
                               "PLR",
                               "NLR",
                               "DOR"))

## Transformation functions

trans_log <- function(x){
  
  x[x == 0] <- .Machine$double.eps
  
  tmp <- log(x)
  
  return(tmp)
  
}

trans_logit <- function(x){
  
  x[x == 0] <- .Machine$double.eps
  
  x[x == 1] <- 1 - .Machine$double.eps
  
  tmp <- logit(x)
  
  return(tmp)
  
  
}

trans_exp <- function(x){
  
  tmp <- exp(x)
  
  return(tmp)
  
}

trans_inv_logit <- function(x){
  
  tmp <- plogis(x)
  
  return(tmp)
  
}

## Chapter

if(length(model_options) == 2){
  
  chapter <- "Ch2"
  
} else {
  
  chapter <- "Ch3"
  
}

#################################################
### Create empty dataframe to store plot info ###
#################################################

n_row <- length(model_options) *
  length(pred_sets) * length(c(ts_site, ts_species))

plot_df <- data.frame(model = character(n_row),
                      ts = character(n_row),
                      mean_diff = numeric(n_row),
                      upper = numeric(n_row),
                      lower = numeric(n_row),
                      pred_type = character(n_row),
                      stringsAsFactors = FALSE)

###################################################
### Extract values to create plotting dataframe ###
###################################################

row_index <- 1

for(prediction in seq_len(length(pred_sets))){
  
  for(ts in c(ts_species, 
              ts_site, 
              "species_richness_difference",
              "independent_log_likelihood",
              "joint_log_likelihood")){
    
    ## Check ts/prediction compatibility
    
    ### Binary test statistics for binary predictions
    
    if(names(pred_sets)[prediction] %in% binary_predictions & ts %nin% binary_ts){
      
      next()
      
    }
    
    ### Probabilistic test statistics for probabilistic predictions
    
    if(names(pred_sets)[prediction] %in% probability_predictions & ts %nin% prob_ts){
      
      next()
      
    }
    
    ## Read in file
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        names(pred_sets[prediction]),
                        ts,
                        chapter)
    
    if(!file.exists(filename)){
      
      message(sprintf("No model for: %s - %s",
                      names(pred_sets[prediction]),
                      ts))
      
      next()
      
    }
    
    mem_model <- readRDS(filename)
    
    if(names(pred_sets)[prediction] %in% probability_predictions){
      
      ## Extract coefficients
      
      coef_matrix <- coef(summary(mem_model))
      
      ## Extract variance/covariance matrix
      
      vc_matrix <- vcov(mem_model)
      
      ## Difference between two means
      
      #### Incorrect Method
      # ### Difference in mean of fixed effects mew(A-B)
      # ### Variance of difference in mean of fixed effects
      # ### var(A-B) = var(A) + var(B) + 2 * covar(A,B)
      # 
      # mean_diff <- coef_matrix["modelSSDM", "Value"] - coef_matrix["modelMPR", "Value"]
      # 
      # var_diff <- vc_matrix["modelSSDM", "modelSSDM"] +
      #   vc_matrix["modelMPR", "modelMPR"] +
      #   2 * vc_matrix["modelSSDM", "modelMPR"]
      # 
      # se_diff <- sqrt(var_diff)
      # 
      # if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
      #   
      #   quantile_fun <- qnorm
      #   
      #   mean <- mean_diff
      #   
      # }
      # 
      # if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
      #   
      #   quantile_fun <- qlnorm
      #   
      #   mean <- exp(mean_diff + (se_diff ^ 2 / 2))
      #   
      # }
      # 
      # if(mem_model$transformed == TRUE & ts %in% ts_0_1){
      #   
      #   quantile_fun <- qlogitnorm
      #   
      #   mean <- mean(plogis(rnorm(1000000, mean_diff, se_diff)))
      #   
      # }
      
      ## New Correct Method. Credit: Martin Ingram
      
      means <- c(coef_matrix["modelSSDM", "Value"],
                 coef_matrix["modelMPR", "Value"],
                 coef_matrix["modelHPR", "Value"],
                 coef_matrix["modelLPR", "Value"],
                 coef_matrix["modelDPR", "Value"],
                 coef_matrix["modelHLR_NS", "Value"],
                 coef_matrix["modelHLR_S", "Value"])
      
      SSDM_mean <- coef_matrix["modelSSDM", "Value"]
      
      SSDM_se <- coef_matrix["modelSSDM", "Std.Error"]
      
      covar_matrix <- matrix(c(vc_matrix["modelSSDM", "modelSSDM"], vc_matrix["modelSSDM", "modelMPR"], vc_matrix["modelSSDM", "modelHPR"], vc_matrix["modelSSDM", "modelLPR"], vc_matrix["modelSSDM", "modelDPR"], vc_matrix["modelSSDM", "modelHLR_NS"], vc_matrix["modelSSDM", "modelHLR_S"],
                               vc_matrix["modelMPR", "modelSSDM"], vc_matrix["modelMPR", "modelMPR"], vc_matrix["modelMPR", "modelHPR"], vc_matrix["modelMPR", "modelLPR"], vc_matrix["modelMPR", "modelDPR"], vc_matrix["modelMPR", "modelHLR_NS"], vc_matrix["modelMPR", "modelHLR_S"],
                               vc_matrix["modelHPR", "modelSSDM"], vc_matrix["modelHPR", "modelMPR"], vc_matrix["modelHPR", "modelHPR"], vc_matrix["modelHPR", "modelLPR"], vc_matrix["modelHPR", "modelDPR"], vc_matrix["modelHPR", "modelHLR_NS"], vc_matrix["modelHPR", "modelHLR_S"],
                               vc_matrix["modelLPR", "modelSSDM"], vc_matrix["modelLPR", "modelMPR"], vc_matrix["modelLPR", "modelHPR"], vc_matrix["modelLPR", "modelLPR"], vc_matrix["modelLPR", "modelDPR"], vc_matrix["modelLPR", "modelHLR_NS"], vc_matrix["modelLPR", "modelHLR_S"],
                               vc_matrix["modelDPR", "modelSSDM"], vc_matrix["modelDPR", "modelMPR"], vc_matrix["modelDPR", "modelHPR"], vc_matrix["modelDPR", "modelLPR"], vc_matrix["modelDPR", "modelDPR"], vc_matrix["modelDPR", "modelHLR_NS"], vc_matrix["modelDPR", "modelHLR_S"],
                               vc_matrix["modelHLR_NS", "modelSSDM"], vc_matrix["modelHLR_NS", "modelMPR"], vc_matrix["modelHLR_NS", "modelHPR"], vc_matrix["modelHLR_NS", "modelLPR"], vc_matrix["modelHLR_NS", "modelDPR"], vc_matrix["modelHLR_NS", "modelHLR_NS"], vc_matrix["modelHLR_NS", "modelHLR_S"],
                               vc_matrix["modelHLR_S", "modelSSDM"], vc_matrix["modelHLR_S", "modelMPR"], vc_matrix["modelHLR_S", "modelHPR"], vc_matrix["modelHLR_S", "modelLPR"], vc_matrix["modelHLR_S", "modelDPR"], vc_matrix["modelHLR_S", "modelHLR_NS"], vc_matrix["modelHLR_S", "modelHLR_S"]),
                             nrow = 7, 
                             ncol = 7,
                             byrow = TRUE)
      
      draws <- mvrnorm(n = 1E6, 
                       means, 
                       covar_matrix)
      
      SSDM_orig_scale <- draws[ , 1]
      
      MPR_orig_scale <- draws[, 2]
      
      HPR_orig_scale <- draws[, 3]
      
      LPR_orig_scale <- draws[, 4]
      
      DPR_orig_scale <- draws[, 5]
      
      HLR_NS_orig_scale <- draws[, 6]
      
      HLR_S_orig_scale <- draws[, 7]
      
      if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
        
        SSDM_coef <- SSDM_mean
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
        
        SSDM_orig_scale <- exp(SSDM_orig_scale)
        
        MPR_orig_scale <- exp(MPR_orig_scale)
        
        HPR_orig_scale <- exp(HPR_orig_scale)
        
        LPR_orig_scale <- exp(LPR_orig_scale)
        
        DPR_orig_scale <- exp(DPR_orig_scale)
        
        HLR_NS_orig_scale <- exp(HLR_NS_orig_scale)
        
        HLR_S_orig_scale <- exp(HLR_S_orig_scale)
        
        SSDM_coef <- exp(SSDM_mean + (SSDM_se ^ 2 / 2)) 
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_1){
        
        SSDM_orig_scale <- plogis(SSDM_orig_scale)
        
        MPR_orig_scale <- plogis(MPR_orig_scale)
        
        HPR_orig_scale <- plogis(HPR_orig_scale)
        
        LPR_orig_scale <- plogis(LPR_orig_scale)
        
        DPR_orig_scale <- plogis(DPR_orig_scale)
        
        HLR_NS_orig_scale <- plogis(HLR_NS_orig_scale)
        
        HLR_S_orig_scale <- plogis(HLR_S_orig_scale)
        
        SSDM_coef <- mean(plogis(rnorm(1000000, SSDM_mean, SSDM_se)))
        
      }
      
      ## Divide by backtransformed regression coefficient here
      
      MPR_diff_samples <- MPR_orig_scale - SSDM_orig_scale
      
      MPR_diff_samples <- (MPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      HPR_diff_samples <- HPR_orig_scale - SSDM_orig_scale
      
      HPR_diff_samples <- (HPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      LPR_diff_samples <- LPR_orig_scale - SSDM_orig_scale
      
      LPR_diff_samples <- (LPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      DPR_diff_samples <- DPR_orig_scale - SSDM_orig_scale
      
      DPR_diff_samples <- (DPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      HLR_NS_diff_samples <- HLR_NS_orig_scale - SSDM_orig_scale
      
      HLR_NS_diff_samples <- (HLR_NS_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      HLR_S_diff_samples <- HLR_S_orig_scale - SSDM_orig_scale
      
      HLR_S_diff_samples <- (HLR_S_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      ## Save to data frame
      
      plot_df[row_index, ] <- list("MPR",
                                   ts,
                                   mean(MPR_diff_samples,
                                        na.rm = TRUE),
                                   quantile(MPR_diff_samples, 0.975,
                                            na.rm = TRUE),
                                   quantile(MPR_diff_samples, 0.025,
                                            na.rm = TRUE),
                                   names(pred_sets)[prediction])
      
      plot_df[row_index + 1, ] <- list("HPR",
                                       ts,
                                       mean(HPR_diff_samples,
                                            na.rm = TRUE),
                                       quantile(HPR_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(HPR_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 2, ] <- list("LPR",
                                       ts,
                                       mean(LPR_diff_samples,
                                            na.rm = TRUE),
                                       quantile(LPR_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(LPR_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 3, ] <- list("DPR",
                                       ts,
                                       mean(DPR_diff_samples,
                                            na.rm = TRUE),
                                       quantile(DPR_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(DPR_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 4, ] <- list("HLR_NS",
                                       ts,
                                       mean(HLR_NS_diff_samples,
                                            na.rm = TRUE),
                                       quantile(HLR_NS_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(HLR_NS_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 5, ] <- list("HLR_S",
                                       ts,
                                       mean(HLR_S_diff_samples,
                                            na.rm = TRUE),
                                       quantile(HLR_S_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(HLR_S_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      row_index <- row_index + 6
      
    }
    
    if(names(pred_sets)[prediction] %in% binary_predictions){
      
      ## Extract coefficients
      
      coef_matrix <- coef(summary(mem_model))
      
      ## Extract variance/covariance matrix
      
      vc_matrix <- vcov(mem_model)
      
      ## Difference between two means
      
      #### Incorrect Method
      # ### Difference in mean of fixed effects mew(A-B)
      # ### Variance of difference in mean of fixed effects
      # ### var(A-B) = var(A) + var(B) + 2 * covar(A,B)
      # 
      # mean_diff <- coef_matrix["modelSSDM", "Value"] - coef_matrix["modelMPR", "Value"]
      # 
      # var_diff <- vc_matrix["modelSSDM", "modelSSDM"] +
      #   vc_matrix["modelMPR", "modelMPR"] +
      #   2 * vc_matrix["modelSSDM", "modelMPR"]
      # 
      # se_diff <- sqrt(var_diff)
      # 
      # if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
      #   
      #   quantile_fun <- qnorm
      #   
      #   mean <- mean_diff
      #   
      # }
      # 
      # if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
      #   
      #   quantile_fun <- qlnorm
      #   
      #   mean <- exp(mean_diff + (se_diff ^ 2 / 2))
      #   
      # }
      # 
      # if(mem_model$transformed == TRUE & ts %in% ts_0_1){
      #   
      #   quantile_fun <- qlogitnorm
      #   
      #   mean <- mean(plogis(rnorm(1000000, mean_diff, se_diff)))
      #   
      # }
      
      ## New Correct Method. Credit: Martin Ingram
      
      means <- c(coef_matrix["modelSSDM", "Value"],
                 coef_matrix["modelMPR", "Value"],
                 coef_matrix["modelHPR", "Value"],
                 coef_matrix["modelLPR", "Value"],
                 coef_matrix["modelDPR", "Value"],
                 coef_matrix["modelHLR_NS", "Value"],
                 coef_matrix["modelHLR_S", "Value"],
                 coef_matrix["modelSESAM", "Value"])
      
      SSDM_mean <- coef_matrix["modelSSDM", "Value"]
      
      SSDM_se <- coef_matrix["modelSSDM", "Std.Error"]
      
      covar_matrix <- matrix(c(vc_matrix["modelSSDM", "modelSSDM"], vc_matrix["modelSSDM", "modelMPR"], vc_matrix["modelSSDM", "modelHPR"], vc_matrix["modelSSDM", "modelLPR"], vc_matrix["modelSSDM", "modelDPR"], vc_matrix["modelSSDM", "modelHLR_NS"], vc_matrix["modelSSDM", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelMPR", "modelSSDM"], vc_matrix["modelMPR", "modelMPR"], vc_matrix["modelMPR", "modelHPR"], vc_matrix["modelMPR", "modelLPR"], vc_matrix["modelMPR", "modelDPR"], vc_matrix["modelMPR", "modelHLR_NS"], vc_matrix["modelMPR", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelHPR", "modelSSDM"], vc_matrix["modelHPR", "modelMPR"], vc_matrix["modelHPR", "modelHPR"], vc_matrix["modelHPR", "modelLPR"], vc_matrix["modelHPR", "modelDPR"], vc_matrix["modelHPR", "modelHLR_NS"], vc_matrix["modelHPR", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelLPR", "modelSSDM"], vc_matrix["modelLPR", "modelMPR"], vc_matrix["modelLPR", "modelHPR"], vc_matrix["modelLPR", "modelLPR"], vc_matrix["modelLPR", "modelDPR"], vc_matrix["modelLPR", "modelHLR_NS"], vc_matrix["modelLPR", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelDPR", "modelSSDM"], vc_matrix["modelDPR", "modelMPR"], vc_matrix["modelDPR", "modelHPR"], vc_matrix["modelDPR", "modelLPR"], vc_matrix["modelDPR", "modelDPR"], vc_matrix["modelDPR", "modelHLR_NS"], vc_matrix["modelDPR", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelHLR_NS", "modelSSDM"], vc_matrix["modelHLR_NS", "modelMPR"], vc_matrix["modelHLR_NS", "modelHPR"], vc_matrix["modelHLR_NS", "modelLPR"], vc_matrix["modelHLR_NS", "modelDPR"], vc_matrix["modelHLR_NS", "modelHLR_NS"], vc_matrix["modelHLR_NS", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelHLR_S", "modelSSDM"], vc_matrix["modelHLR_S", "modelMPR"], vc_matrix["modelHLR_S", "modelHPR"], vc_matrix["modelHLR_S", "modelLPR"], vc_matrix["modelHLR_S", "modelDPR"], vc_matrix["modelHLR_S", "modelHLR_NS"], vc_matrix["modelHLR_S", "modelHLR_S"], vc_matrix["modelSSDM", "modelSESAM"],
                               vc_matrix["modelSESAM", "modelSSDM"], vc_matrix["modelSESAM", "modelMPR"], vc_matrix["modelSESAM", "modelHPR"], vc_matrix["modelSESAM", "modelLPR"], vc_matrix["modelSESAM", "modelDPR"], vc_matrix["modelSESAM", "modelHLR_NS"], vc_matrix["modelSESAM", "modelHLR_S"], vc_matrix["modelSESAM", "modelSESAM"]),
                             nrow = 8, 
                             ncol = 8,
                             byrow = TRUE)
      
      draws <- mvrnorm(n = 1E6, 
                       means, 
                       covar_matrix)
      
      SSDM_orig_scale <- draws[ , 1]
      
      MPR_orig_scale <- draws[, 2]
      
      HPR_orig_scale <- draws[, 3]
      
      LPR_orig_scale <- draws[, 4]
      
      DPR_orig_scale <- draws[, 5]
      
      HLR_NS_orig_scale <- draws[, 6]
      
      HLR_S_orig_scale <- draws[, 7]
      
      SESAM_orig_scale <- draws[, 8]
      
      if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
        
        SSDM_coef <- SSDM_mean
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
        
        SSDM_orig_scale <- exp(SSDM_orig_scale)
        
        MPR_orig_scale <- exp(MPR_orig_scale)
        
        HPR_orig_scale <- exp(HPR_orig_scale)
        
        LPR_orig_scale <- exp(LPR_orig_scale)
        
        DPR_orig_scale <- exp(DPR_orig_scale)
        
        HLR_NS_orig_scale <- exp(HLR_NS_orig_scale)
        
        HLR_S_orig_scale <- exp(HLR_S_orig_scale)
        
        SESAM_orig_scale <- exp(SESAM_orig_scale)
        
        SSDM_coef <- exp(SSDM_mean + (SSDM_se ^ 2 / 2)) 
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_1){
        
        SSDM_orig_scale <- plogis(SSDM_orig_scale)
        
        MPR_orig_scale <- plogis(MPR_orig_scale)
        
        HPR_orig_scale <- plogis(HPR_orig_scale)
        
        LPR_orig_scale <- plogis(LPR_orig_scale)
        
        DPR_orig_scale <- plogis(DPR_orig_scale)
        
        HLR_NS_orig_scale <- plogis(HLR_NS_orig_scale)
        
        HLR_S_orig_scale <- plogis(HLR_S_orig_scale)
        
        SESAM_orig_scale <- plogis(SESAM_orig_scale)
        
        SSDM_coef <- mean(plogis(rnorm(1000000, SSDM_mean, SSDM_se)))
        
      }
      
      ## Divide by backtransformed regression coefficient here
      
      MPR_diff_samples <- MPR_orig_scale - SSDM_orig_scale
      
      MPR_diff_samples <- (MPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      HPR_diff_samples <- HPR_orig_scale - SSDM_orig_scale
      
      HPR_diff_samples <- (HPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      LPR_diff_samples <- LPR_orig_scale - SSDM_orig_scale
      
      LPR_diff_samples <- (LPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      DPR_diff_samples <- DPR_orig_scale - SSDM_orig_scale
      
      DPR_diff_samples <- (DPR_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      HLR_NS_diff_samples <- HLR_NS_orig_scale - SSDM_orig_scale
      
      HLR_NS_diff_samples <- (HLR_NS_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      HLR_S_diff_samples <- HLR_S_orig_scale - SSDM_orig_scale
      
      HLR_S_diff_samples <- (HLR_S_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      SESAM_diff_samples <- SESAM_orig_scale - SSDM_orig_scale
      
      SESAM_diff_samples <- (SESAM_diff_samples / abs(SSDM_coef)) * 100 #- 100
      
      ## Save to data frame
      
      plot_df[row_index, ] <- list("MPR",
                                   ts,
                                   mean(MPR_diff_samples,
                                        na.rm = TRUE),
                                   quantile(MPR_diff_samples, 0.975,
                                            na.rm = TRUE),
                                   quantile(MPR_diff_samples, 0.025,
                                            na.rm = TRUE),
                                   names(pred_sets)[prediction])
      
      plot_df[row_index + 1, ] <- list("HPR",
                                       ts,
                                       mean(HPR_diff_samples,
                                            na.rm = TRUE),
                                       quantile(HPR_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(HPR_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 2, ] <- list("LPR",
                                       ts,
                                       mean(LPR_diff_samples,
                                            na.rm = TRUE),
                                       quantile(LPR_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(LPR_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 3, ] <- list("DPR",
                                       ts,
                                       mean(DPR_diff_samples,
                                            na.rm = TRUE),
                                       quantile(DPR_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(DPR_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 4, ] <- list("HLR_NS",
                                       ts,
                                       mean(HLR_NS_diff_samples,
                                            na.rm = TRUE),
                                       quantile(HLR_NS_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(HLR_NS_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 5, ] <- list("HLR_S",
                                       ts,
                                       mean(HLR_S_diff_samples,
                                            na.rm = TRUE),
                                       quantile(HLR_S_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(HLR_S_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      plot_df[row_index + 6, ] <- list("SESAM",
                                       ts,
                                       mean(SESAM_diff_samples,
                                            na.rm = TRUE),
                                       quantile(SESAM_diff_samples, 0.975,
                                                na.rm = TRUE),
                                       quantile(SESAM_diff_samples, 0.025,
                                                na.rm = TRUE),
                                       names(pred_sets)[prediction])
      
      row_index <- row_index + 7
      
    }
  }
}

plot_df <- plot_df[plot_df$ts != "", ]

plot_df$ts_class <- NA
  
for(i in seq_len(nrow(plot_df))){
  
  plot_df$ts_class[i] <- names(subsets)[unlist(lapply(subsets, function(x) plot_df$ts[i] %in% x))]

}

################
### Plotting ###
################

dodge <- position_dodge(width = 0.5)

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


## Plot

for(pred_type in unique(plot_df$pred_type)){
    
    for(subset in names(subsets)){
      
      tmp_df <- plot_df[plot_df$pred_type == pred_type &
                          plot_df$ts %in% subsets[[subset]], ]
      
      if(nrow(tmp_df) == 0){
        
        next()
        
      }
      
      tmp_df$ts <- factor(tmp_df$ts,
                          levels = sort(subsets[[subset]]))
      
      ## Make plot and save to file
      
      filename <- sprintf("outputs/test_statistics/plots/relative_performance_%1$s_%2$s_%3$s_MME.pdf",
                          pred_type,
                          subset,
                          chapter)
      
      pdf(filename)
      
      ### Start plot
      
      tmp_plot <- ggplot(tmp_df,
                         aes(x = ts,
                             y = mean_diff,
                             colour = model)) +
        geom_hline(yintercept = 0) +
        geom_point(position = dodge,
                   size = 2) + 
        geom_errorbar(aes(ymax = upper,
                          ymin = lower,),
                      position = dodge,
                      width = 0.1) +
        xlab("") +
        ylab("Relative Performance Difference (%)") +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text.x = element_text(angle = 45,
                                         hjust = 1))
      
      ### Axis labels
      
      if(subset == "ts_0_1_H"){
        
        tmp_plot <- tmp_plot +
          scale_x_discrete(labels = c("Accuracy" = "Accuracy",
                                      "AUC" = "AUC",
                                      "Kappa" = "Cohen's Kappa",
                                      "F_1" = eval(parse(text = "expression(F[1])")),
                                      "PPV" = "Positive predictive performance",
                                      "NPV" = "Negative predictive performance",
                                      "R2" = eval(parse(text = "expression(R^2)")),
                                      "TNR" = "True negative rate",
                                      "TPR" = "True positive rate",
                                      "Youden_J" = "Youden's J"),
                           drop = FALSE)
        
      }
      
      if(subset == "ts_0_1_L_a"){
        
        tmp_plot <- tmp_plot + 
          scale_x_discrete(labels = c("FDR" = "False discovery rate",
                                      "FNR" = "False negative rate",
                                      "FOR" = "False omission rate",
                                      "FPR" = "False positive rate"),
                           drop = FALSE)
        
      }
      
      if(subset == "ts_0_1_L_b"){
        
        tmp_plot <- tmp_plot + 
          scale_x_discrete(labels = c("Bray" = "Bray-Curtis dissimilarity",
                                      "Canberra" = "Canberra index",
                                      "Gower" = "Gower index",
                                      "Gower_alt" = "Gower index (alternative)",
                                      "Jaccard" = "Jaccard index",
                                      "Kulczynski" = "Kulczynski index",
                                      "Mountford" = "Mountford index",
                                      "Raup" = "Raup-Crick dissimilarity"),
                           drop = FALSE)
        
      }
      
      if(subset == "ts_1_1"){
        
        tmp_plot <- tmp_plot +
          scale_x_discrete(labels = c("Kendall" = "Kendall rank correlation coefficient",
                                      "Spearman" = "Spearman's rank correlation coefficient",
                                      "Pearson" = "Pearson correlation coefficient"),
                           drop = FALSE)
        
      }
      
      if(subset == "ts_Inf_Inf_H"){
        
        tmp_plot <- tmp_plot +
          scale_x_discrete(labels = c("independent_log_likelihood" = "Independent log likelihood",
                                      "joint_log_likelihood" = "Joint log likelihood"),
                           drop = FALSE)
        
      }
      
      if(subset == "ts_Inf_Inf_0"){
        
        tmp_plot <- tmp_plot +
          scale_x_discrete(labels = c("bias" = "Mean error",
                                      "species_richness_difference" = "Species richness difference"),
                           drop = FALSE)
        
      }
      
      if(subset == "ts_0_Inf_L"){
        
        tmp_plot <- tmp_plot +
          scale_x_discrete(labels= c("MSE" = "Mean square error",
                                     "RMSE" = "Root mean square error",
                                     "SSE" = "Sum of squared errors"),
                           drop = FALSE)
        
      }
      
      ### Main title
      
      if(pred_type == "marginal_bin"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Marginal (binary)")
        
      }
      
      if(pred_type == "marginal_prob"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Marginal (probabilistic)")
        
      }
      
      if(pred_type == "condLOI_low"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Conditional (Low)")
        
      }
      
      if(pred_type == "condLOI_med"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Conditional (Medium)")
        
      }
      
      if(pred_type == "condLOI_high"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Conditional (High)")
        
      }
      
      if(pred_type == "condLOI_marg_low"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Conditional marginal (Low)")
        
      }
      
      if(pred_type == "condLOI_marg_med"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Conditional marginal (Medium)")
        
      }
      
      if(pred_type == "condLOI_marg_high"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Conditional marginal (High)")
        
      }
      
      if(pred_type == "joint"){
        
        tmp_plot <- tmp_plot +
          ggtitle(label = "Relative test statistic performance (SSDM - JSDM)",
                  subtitle = "Prediction type: Joint")
        
      }
      
      y_range <- ggplot_build(tmp_plot)$layout$panel_scales_y[[1]]$range$range
      
      if(subset %in% names(subsets)[c(2,3,7)]){
        
        # tmp_plot <- tmp_plot +
        #   annotate("segment", 
        #            x = 0.8, 
        #            xend = 0.8, 
        #            y = 0 - diff_range / 10, 
        #            yend = 0 + diff_range / 10, 
        #            colour = "black", 
        #            size = 0.5,
        #            arrow = arrow(ends = "last",
        #                          length = unit(0.1, "inches")))
        # 
        tmp_plot <- tmp_plot +
          annotate("text",
                   label = "JSDM better",
                   x = 0.7,
                   y = y_range[2] - ((y_range[2] - y_range[1]) / 10) * 9)
        
      }
      
      if(subset %in% names(subsets)[c(1,4,5)]){
        
        # tmp_plot <- tmp_plot +
        #   annotate("segment", 
        #            x = 0.5, 
        #            xend = 0.5, 
        #            y = 0 - diff_range / 10, 
        #            yend = 0 + diff_range / 10, 
        #            colour = "black", 
        #            size = 0.5,
        #            arrow = arrow(ends = "first",
        #                          length = unit(0.1, "inches")))
        
        tmp_plot <- tmp_plot +
          annotate("text",
                   label = "JSDM better",
                   x = 0.7,
                   y = y_range[1] + ((y_range[2] - y_range[1]) / 10) * 9)
        
      }
      
      if(subset %in% names(subsets)[c(6)]){
        
      }
      
      print(tmp_plot)
      
      dev.off()
      
    }
  }