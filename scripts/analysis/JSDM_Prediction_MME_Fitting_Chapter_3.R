


#####################
### Load Packages ###
#####################

library(nlme)
library(ggplot2)
library(RColorBrewer)
library(PassButter)
library(psych)
library(logitnorm)
library(AICcmodavg)

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
                     # "bird",
                     # "butterfly",
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

################################################
################################################
### LOAD ALL DATA AND CALCULATE BONFERRONI P ###
################################################
################################################

#######################
### Test Statistics ###
#######################

## Load Data

ts_df_species <- readRDS(sprintf("outputs/test_statistics/test_statistics_species_summary_MME_%s.rds",
                                 chapter))

ts_df_site <- readRDS(sprintf("outputs/test_statistics/test_statistics_site_summary_MME_%s.rds",
                              chapter))

ts_df_species$fold <- as.factor(ts_df_species$fold)

ts_df_site$fold <- as.factor(ts_df_site$fold)

ts_df_species$model <- as.factor(ts_df_species$model)

ts_df_species$model <- relevel(ts_df_species$model, "SSDM")

ts_df_site$model <- as.factor(ts_df_site$model)

ts_df_site$model <- relevel(ts_df_site$model, "SSDM")

ts_df_species <- ts_df_species[order(ts_df_species$model), ]

ts_df_site <- ts_df_site[order(ts_df_site$model), ]

ts_df_species$SSDM <- as.numeric(ts_df_species$model == "SSDM")

ts_df_species$SESAM <- as.numeric(ts_df_species$model == "SESAM")

ts_df_species$MPR <- as.numeric(ts_df_species$model == "MPR")

ts_df_species$HPR <- as.numeric(ts_df_species$model == "HPR")

ts_df_species$LPR <- as.numeric(ts_df_species$model == "LPR")

ts_df_species$DPR <- as.numeric(ts_df_species$model == "DPR")

ts_df_species$HLR_NS <- as.numeric(ts_df_species$model == "HLR_NS")

ts_df_species$HLR_S <- as.numeric(ts_df_species$model == "HLR_S")

ts_df_site$SSDM <- as.numeric(ts_df_site$model == "SSDM")

ts_df_site$SESAM <- as.numeric(ts_df_site$model == "SESAM")

ts_df_site$MPR <- as.numeric(ts_df_site$model == "MPR")

ts_df_site$HPR <- as.numeric(ts_df_site$model == "HPR")

ts_df_site$LPR <- as.numeric(ts_df_site$model == "LPR")

ts_df_site$DPR <- as.numeric(ts_df_site$model == "DPR")

ts_df_site$HLR_NS <- as.numeric(ts_df_site$model == "HLR_NS")

ts_df_site$HLR_S <- as.numeric(ts_df_site$model == "HLR_S")

## Empty dataframe to store Kolmogorov-Smirnov test outputs

n_row_ts <- (length(unique(ts_df_species$test_statistic)) + 
               length(unique(ts_df_site$test_statistic))) *
  length(pred_sets)

ks_df <- data.frame(pred_type = character(n_row_ts),
                    test_statistic = character(n_row_ts),
                    p_value = numeric(n_row_ts),
                    ks_D = numeric(n_row_ts),
                    model_fit_success = numeric(n_row_ts),
                    stringsAsFactors = FALSE)

########################
### Species Richness ###
########################

## Load Data

sr_df <- readRDS(sprintf("outputs/species_richness/species_richness_summary_%s.rds",
                         chapter))

sr_df$fold <- as.factor(sr_df$fold)

sr_df$model <- as.factor(sr_df$model)

sr_df$model <- relevel(sr_df$model, "SSDM")

sr_df$SSDM <- as.numeric(sr_df$model == "SSDM")

sr_df$SESAM <- as.numeric(sr_df$model == "SESAM")

sr_df$MPR <- as.numeric(sr_df$model == "MPR")

sr_df$HPR <- as.numeric(sr_df$model == "HPR")

sr_df$LPR <- as.numeric(sr_df$model == "LPR")

sr_df$DPR <- as.numeric(sr_df$model == "DPR")

sr_df$HLR_NS <- as.numeric(sr_df$model == "HLR_NS")

sr_df$HLR_S <- as.numeric(sr_df$model == "HLR_S")

## Empty dataframe to store Kolmogorov-Smirnov test outputs

n_row_sr <- length(unique(sr_df$test_statistic)) *
  length(pred_sets)

ks_sr <- data.frame(pred_type = character(n_row_sr),
                    test_statistic = character(n_row_sr),
                    p_value = numeric(n_row_sr),
                    ks_D = numeric(n_row_sr),
                    model_fit_success = numeric(n_row_sr),
                    stringsAsFactors = FALSE)

##################
### Likelihood ###
##################

## Load Data

ll_i_df <- readRDS(sprintf("outputs/likelihood/independent_likelihood_summary_%s.rds",
                           chapter))

ll_j_df <- readRDS(sprintf("outputs/likelihood/joint_likelihood_summary_%s.rds",
                           chapter))

ll_i_df$fold <- as.factor(ll_i_df$fold)

ll_j_df$fold <- as.factor(ll_j_df$fold)

ll_i_df$model <- as.factor(ll_i_df$model)

ll_i_df$model <- relevel(ll_i_df$model, "SSDM")

ll_j_df$model <- as.factor(ll_j_df$model)

ll_j_df$model <- relevel(ll_j_df$model, "SSDM")

ll_i_df$SSDM <- as.numeric(ll_i_df$model == "SSDM")

ll_i_df$SESAM <- as.numeric(ll_i_df$model == "SESAM")

ll_i_df$MPR <- as.numeric(ll_i_df$model == "MPR")

ll_i_df$HPR <- as.numeric(ll_i_df$model == "HPR")

ll_i_df$LPR <- as.numeric(ll_i_df$model == "LPR")

ll_i_df$DPR <- as.numeric(ll_i_df$model == "DPR")

ll_i_df$HLR_NS <- as.numeric(ll_i_df$model == "HLR_NS")

ll_i_df$HLR_S <- as.numeric(ll_i_df$model == "HLR_S")

ll_j_df$SSDM <- as.numeric(ll_j_df$model == "SSDM")

ll_j_df$SESAM <- as.numeric(ll_j_df$model == "SESAM")

ll_j_df$MPR <- as.numeric(ll_j_df$model == "MPR")

ll_j_df$HPR <- as.numeric(ll_j_df$model == "HPR")

ll_j_df$LPR <- as.numeric(ll_j_df$model == "LPR")

ll_j_df$DPR <- as.numeric(ll_j_df$model == "DPR")

ll_j_df$HLR_NS <- as.numeric(ll_j_df$model == "HLR_NS")

ll_j_df$HLR_S <- as.numeric(ll_j_df$model == "HLR_S")

## Empty dataframe to store Kolmogorov-Smirnov test outputs

n_row_ll <- (length(unique(ll_i_df$test_statistic)) *
               length(pred_sets)) +
  (length(unique(ll_j_df$test_statistic)) *
     length(pred_sets))

ks_ll <- data.frame(pred_type = character(n_row_ll),
                    test_statistic = character(n_row_ll),
                    p_value = numeric(n_row_ll),
                    ks_D = numeric(n_row_ll),
                    model_fit_success = numeric(n_row_ll),
                    stringsAsFactors = FALSE)

###############################
### Bonferroni P correction ###
###############################

bonferroni_p <- 0.05 / (n_row_ts + n_row_sr + n_row_ll)

#######################
#######################
### TEST STATISTICS ###
#######################
#######################

############################
### Mixed Effects Models ###
############################

## Loop over different model iterations

row_index <- 1

## Species-based statistics

for(prediction in seq_len(length(pred_sets))){

  for(ts in ts_species){
    
    ## Check ts/prediction compatibility
    
    ### Binary test statistics for binary predictions
    
    if(names(pred_sets)[prediction] %in% binary_predictions & ts %nin% binary_ts){
      
      next()
      
    }
    
    ### Probabilistic test statistics for probabilistic predictions
    
    if(names(pred_sets)[prediction] %in% probability_predictions & ts %nin% prob_ts){
      
      next()
      
    }
    
    ## Subset dataset
    
    tmp_df <- ts_df_species[ts_df_species$prediction_type %in% pred_sets[[prediction]] &
                              ts_df_species$test_statistic == ts, ]
    
    tmp_df$species <- as.factor(tmp_df$species)
    
    ## Fit mixed effects model
    
    tmp_df <- na.omit(tmp_df)
    
    tmp_df <- tmp_df[is.finite(tmp_df$mean), ]
    
    tmp_df$mean_trans <- tmp_df$mean
    
    if(ts %in% prob_ts & names(pred_sets)[prediction] %in% probability_predictions){
      
      formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset + DPR * dataset")
    
    }
    
    if(ts %in% binary_ts & names(pred_sets)[prediction] %in% binary_predictions){
      
      formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + SESAM + HPR * dataset + DPR * dataset")
    
    }
    
    mem_model <- tryCatch(expr = nlme::lme(formula, 
                                           random = ~ 1|species,
                                           weights = varIdent(form = ~1|model),
                                           data = tmp_df,
                                           control = lmeControl(msMaxIter = 1000,
                                                                opt = "optim")),
                          error = function(err){
                            
                            message(sprintf("Model fit failed for: %s - %s",
                                            names(pred_sets[prediction]),
                                            ts))
                            
                            return(NA)
                            
                          }
    )
    
    ## If no model fit, skip remainder
    
    if(class(mem_model) != "lme"){
      
      ks_df[row_index, ] <- list(names(pred_sets[prediction]),
                                 ts,
                                 NA,
                                 NA,
                                 0)
      
      row_index <- row_index + 1
      
      next()
      
    }
    
    mem_model$transformed <- FALSE
    
    ## Test for normality to see if we need to transform
    
    ### Extract residuals
    
    residuals <- resid(mem_model, type = "pearson")
    
    ### Correct residuals for variance structure
    
    cf <- coef(mem_model$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ### Check normality with Kolmogorov-Smirnov test
    
    extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    ks <- ks.test(test_resid, 
                  pnorm, 
                  mean(test_resid), 
                  sd(test_resid))
    
    ## If fail KS-test, re-fit model with transformation
    
    if(ks$p.value < bonferroni_p){
      
      if(ts %in% ts_Inf_Inf){
        
        tmp_df$mean_trans <- tmp_df$mean
        
      }
      
      if(ts %in% ts_0_Inf){
        
        tmp_df$mean_trans <- trans_log(tmp_df$mean)
        
      }
      
      if(ts %in% ts_0_1){
        
        tmp_df$mean_trans <- trans_logit(tmp_df$mean)
        
      }
      
      tmp_df <- na.omit(tmp_df)
      
      tmp_df <- tmp_df[is.finite(tmp_df$mean), ]
      
      tmp_df$mean_trans <- tmp_df$mean
      
      if(ts %in% prob_ts & names(pred_sets)[prediction] %in% probability_predictions){
        
        formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset + DPR * dataset")
        
      }
      
      if(ts %in% binary_ts & names(pred_sets)[prediction] %in% binary_predictions){
        
        formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + SESAM + HPR * dataset + DPR * dataset")
        
      }
      
      mem_model2 <- tryCatch(expr = nlme::lme(formula,
                                             random = ~ 1|species,
                                             weights = varIdent(form = ~1|model),
                                             data = tmp_df,
                                             control = lmeControl(msMaxIter = 1000,
                                                                  opt = "optim")),
                            error = function(err){
                              
                              message(sprintf("Model fit failed for: %s - %s",
                                              names(pred_sets[prediction]),
                                              ts))
                              
                              return(NA)
                              
                            }
      )
      
      mem_model2$transformed <- TRUE
      
    
    ## Check if transformed is better or not
    
    ### Extract residuals
    
    residuals <- resid(mem_model2, type = "pearson")
    
    ### Correct residuals for variance structure
    
    cf <- coef(mem_model2$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ### Residuals qq plot
    
    extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    ks2 <- ks.test(test_resid, 
                   pnorm, 
                   mean(test_resid), 
                   sd(test_resid))
    
    if(ks2$p.value > ks$p.value){
      
      mem_model <- mem_model2
      
    }
    
    }
    
    ## Extract residuals
    
    residuals <- resid(mem_model, type = "pearson")
    
    ## Correct residuals for variance structure
    
    cf <- coef(mem_model$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ## Plot residuals and save to file
    
    ### Set filename
    
    if(length(model_options) == 2){
      
      chapter <- "Ch2"
      
    } else {
      
      chapter <- "Ch3"
      
    }
    
    filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s.pdf",
                        names(pred_sets[prediction]),
                        ts,
                        chapter)
    
    ### Plot to pdf
    
    pdf(filename)
    
    ### Residuals boxplot
    
    plot(residuals ~ factor(tmp_df$model),
         ylab = "Pearson residuals",
         xlab = "Model")
    
    ### Residuals qq plot
    
    extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    qqplot(test_resid, 
           rnorm(n = 1000, 
                 mean(test_resid), 
                 sd(test_resid)),
           ylab = "N(mean(residuals), sd(residuals))",
           xlab = "Model Pearson residuals",
           main = "QQ-plot")
    abline(0,1)
    
    dev.off()
    
    ### Check normality with Kolmogorov-Smirnov test
    
    ks <- ks.test(test_resid, 
                  pnorm, 
                  mean(test_resid), 
                  sd(test_resid))
    
    ks_df[row_index, ] <- list(names(pred_sets[prediction]),
                               ts,
                               ks$p.value,
                               ks$statistic,
                               1)
    ### Save model to file
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        names(pred_sets[prediction]),
                        ts,
                        chapter)
    
    saveRDS(mem_model,
            filename)
    
    ### Statusbar
    
    statusbar(run = row_index,
              max.run = (n_row_ts + n_row_sr + n_row_ll),
              info = sprintf("%s: %s - %-15s",
                             row_index,
                             names(pred_sets[prediction]),
                             ts))
    
    row_index <- row_index + 1
    
    suppressWarnings(rm(ks))
    suppressWarnings(rm(ks2))
    
  }
}

## Site-based statistics

for(prediction in seq_len(length(pred_sets))){

  for(ts in ts_site){
    
    ## Check ts/prediction compatibility
    
    ### Binary test statistics for binary predictions
    
    if(names(pred_sets)[prediction] %in% binary_predictions & ts %nin% binary_ts){
      
      next()
      
    }
    
    ### Probabilistic test statistics for probabilistic predictions
    
    if(names(pred_sets)[prediction] %in% probability_predictions & ts %nin% prob_ts){
      
      next()
      
    }
    
    ## Subset dataset
    
    tmp_df <- ts_df_site[ts_df_site$prediction_type %in% pred_sets[[prediction]] &
                           ts_df_site$test_statistic == ts, ]
    
    tmp_df$site <- as.factor(tmp_df$site)
    
    ## Fit mixed effects model
    
    tmp_df <- na.omit(tmp_df)
    
    tmp_df <- tmp_df[is.finite(tmp_df$mean), ]
    
    tmp_df$mean_trans <- tmp_df$mean
    
    if(ts %in% prob_ts & names(pred_sets)[prediction] %in% probability_predictions){
      
      formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset + DPR * dataset")
      
    }
    
    if(ts %in% binary_ts & names(pred_sets)[prediction] %in% binary_predictions){
      
      formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + SESAM + HPR * dataset + DPR * dataset")
      
    }
    
    mem_model <- tryCatch(expr = nlme::lme(formula,
                                           random = ~ 1|site,
                                           weights = varIdent(form = ~1|model),
                                           data = tmp_df,
                                           control = lmeControl(msMaxIter = 1000,
                                                                opt = "optim")),
                          error = function(err){
                            
                            message(sprintf("Model fit failed for: %s - %s",
                                            names(pred_sets[prediction]),
                                            ts))
                            
                            return(NA)
                            
                          }
    )
    
    ## If no model fit, skip remainder
    
    if(class(mem_model) != "lme"){
      
      ks_df[row_index, ] <- list(names(pred_sets[prediction]),
                                 ts,
                                 NA,
                                 NA,
                                 0)
      
      row_index <- row_index + 1
      
      next()
      
    }
    
    mem_model$transformed <- FALSE
    
    ## Test for normality to see if we need to transform
    
    ### Extract residuals
    
    residuals <- resid(mem_model, type = "pearson")
    
    ### Correct residuals for variance structure
    
    cf <- coef(mem_model$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ### Check normality with Kolmogorov-Smirnov test
    
    extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    ks <- ks.test(test_resid, 
                  pnorm, 
                  mean(test_resid), 
                  sd(test_resid))
    
    ## If fail KS-test, re-fit model with transformation
    
    if(ks$p.value < bonferroni_p){
      
      if(ts %in% ts_Inf_Inf){
        
        tmp_df$mean_trans <- tmp_df$mean
        
      }
      
      if(ts %in% ts_0_Inf){
        
        tmp_df$mean_trans <- trans_log(tmp_df$mean)
        
      }
      
      if(ts %in% ts_0_1){
        
        tmp_df$mean_trans <- trans_logit(tmp_df$mean)
        
      }
      
      tmp_df <- na.omit(tmp_df)
      
      tmp_df <- tmp_df[is.finite(tmp_df$mean), ]
      
      tmp_df$mean_trans <- tmp_df$mean
      
      if(ts %in% prob_ts & names(pred_sets)[prediction] %in% probability_predictions){
        
        formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset + DPR * dataset")
        
      }
      
      if(ts %in% binary_ts & names(pred_sets)[prediction] %in% binary_predictions){
        
        formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + SESAM + HPR * dataset + DPR * dataset")
        
      }
      
      mem_model2 <- tryCatch(expr = nlme::lme(formula,
                                              random = ~ 1|site,
                                              weights = varIdent(form = ~1|model),
                                              data = tmp_df,
                                              control = lmeControl(msMaxIter = 1000,
                                                                   opt = "optim")),
                             error = function(err){
                               
                               message(sprintf("Model fit failed for: %s - %s",
                                               names(pred_sets[prediction]),
                                               ts))
                               
                               return(NA)
                               
                             }
      )
      
      mem_model2$transformed <- TRUE
      
    
    ## Check if transformed is better or not
    
    ### Extract residuals
    
    residuals <- resid(mem_model2, type = "pearson")
    
    ### Correct residuals for variance structure
    
    cf <- coef(mem_model2$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ### Residuals qq plot
    
    extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    ks2 <- ks.test(test_resid, 
                   pnorm, 
                   mean(test_resid), 
                   sd(test_resid))
    
    if(ks2$p.value > ks$p.value){
      
      mem_model <- mem_model2
      
    }
    
    }
    ## Extract residuals
    
    residuals <- resid(mem_model, type = "pearson")
    
    ## Correct residuals for variance structure
    
    cf <- coef(mem_model$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ## Plot residuals and save to file
    
    ### Set filename
    
    if(length(model_options) == 2){
      
      chapter <- "Ch2"
      
    } else {
      
      chapter <- "Ch3"
      
    }
    
    filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s.pdf",
                        names(pred_sets[prediction]),
                        ts,
                        chapter)
    
    ### Plot to pdf
    
    pdf(filename)
    
    ### Residuals boxplot
    
    plot(residuals ~ factor(tmp_df$model),
         ylab = "Pearson residuals",
         xlab = "Model")
    
    ### Residuals qq plot
    
    extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    qqplot(test_resid, 
           rnorm(n = 1000, 
                 mean(test_resid), 
                 sd(test_resid)),
           ylab = "N(mean(residuals), sd(residuals))",
           xlab = "Model Pearson residuals",
           main = "QQ-plot")
    abline(0,1)
    
    dev.off()
    
    ### Check normality with Kolmogorov-Smirnov test
    
    ks <- ks.test(test_resid, 
                  pnorm, 
                  mean(test_resid), 
                  sd(test_resid))
    
    ks_df[row_index, ] <- list(names(pred_sets[prediction]),
                               ts,
                               ks$p.value,
                               ks$statistic,
                               1)
    ### Save model to file
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        names(pred_sets[prediction]),
                        ts,
                        chapter)
    
    saveRDS(mem_model,
            filename)
    
    ### Statusbar
    
    statusbar(run = row_index,
              max.run = (n_row_ts + n_row_sr + n_row_ll),
              info = sprintf("%s: %s - %-15s",
                             row_index,
                             names(pred_sets[prediction]),
                             ts))
    
    row_index <- row_index + 1
    
    suppressWarnings(rm(ks))
    suppressWarnings(rm(ks2))
  }
}

########################
########################
### SPECIES RICHNESS ###
########################
########################

## Loop over different different model iterations

row_index <- 1

for(prediction in seq_len(length(pred_sets))){

  ts <- "species_richness_difference"
  
  ## Subset dataset
  
  tmp_df <- sr_df[sr_df$prediction_type %in% pred_sets[[prediction]], ]
  
  tmp_df <- na.omit(tmp_df)
  
  tmp_df <- tmp_df[is.finite(tmp_df$mean), ]
  
  ## Fit mixed effects model
  
  tmp_df$mean_trans <- tmp_df$mean
  
  if(ts %in% prob_ts & names(pred_sets)[prediction] %in% probability_predictions){
    
    formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset + DPR * dataset")
    
  }
  
  if(ts %in% binary_ts & names(pred_sets)[prediction] %in% binary_predictions){
    
    formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + SESAM + HPR * dataset + DPR * dataset")
    
  }
  
  mem_model <- tryCatch(expr = nlme::lme(formula,
                                         random = ~ 1|site,
                                         weights = varIdent(form = ~1|model),
                                         data = tmp_df,
                                         control = lmeControl(msMaxIter = 1000,
                                                              opt = "optim")),
                        error = function(err){
                          
                          message(sprintf("Model fit failed for: %s - %s",
                                          names(pred_sets[prediction]),
                                          "species_richness_difference"))
                          
                          return(NA)
                          
                        }
  )
  
  ## If no model fit, skip remainder
  
  if(class(mem_model) != "lme"){
    
    ks_sr[row_index, ] <- list(names(pred_sets[prediction]),
                               ts,
                               NA,
                               NA,
                               0)
    
    row_index <- row_index + 1
    
    next()
    
  }
  
  mem_model$transformed <- FALSE
  
  ## Test for normality to see if we need to transform
  
  ### Extract residuals
  
  residuals <- resid(mem_model, type = "pearson")
  
  ### Correct residuals for variance structure
  
  cf <- coef(mem_model$modelStruct$varStruct, 
             unconstrained = FALSE)
  
  for(i in names(cf)){
    
    idx <- tmp_df$model == i
    
    residuals[idx] <- residuals[idx] / cf[i]
    
  }
  
  ### Check normality with Kolmogorov-Smirnov test
  
  extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
  
  if(length(extreme_id) > 0){
    
    test_resid <- residuals[-extreme_id]
    
  } else {
    
    test_resid <- residuals
    
  }
  
  ks <- ks.test(test_resid, 
                pnorm, 
                mean(test_resid), 
                sd(test_resid))
  
  ## If fail KS-test, re-fit model with transformation
  
  if(ks$p.value < bonferroni_p){
    
    if(ts %in% ts_Inf_Inf){
      
      tmp_df$mean_trans <- tmp_df$mean
      
    }
    
    if(ts %in% ts_0_Inf){
      
      tmp_df$mean_trans <- trans_log(tmp_df$mean)
      
    }
    
    if(ts %in% ts_0_1){
      
      tmp_df$mean_trans <- trans_logit(tmp_df$mean)
      
    }
    
    if(ts %in% prob_ts & names(pred_sets)[prediction] %in% probability_predictions){
      
      formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset + DPR * dataset")
      
    }
    
    if(ts %in% binary_ts & names(pred_sets)[prediction] %in% binary_predictions){
      
      formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + SESAM + HPR * dataset + DPR * dataset")
      
    }
    
    mem_model2 <- tryCatch(expr = nlme::lme(formula,
                                           random = ~ 1|site,
                                           weights = varIdent(form = ~1|model),
                                           data = tmp_df,
                                           control = lmeControl(msMaxIter = 1000,
                                                                opt = "optim")),
                          error = function(err){
                            
                            message(sprintf("Model fit failed for: %s - %s",
                                            names(pred_sets[prediction]),
                                            "species_richness_difference"))
                            
                            return(NA)
                            
                          }
    )
    
    mem_model2$transformed <- TRUE
    
  
  ## Check if transformed is better or not
  
  ### Extract residuals
  
  residuals <- resid(mem_model2, type = "pearson")
  
  ### Correct residuals for variance structure
  
  cf <- coef(mem_model2$modelStruct$varStruct, 
             unconstrained = FALSE)
  
  for(i in names(cf)){
    
    idx <- tmp_df$model == i
    
    residuals[idx] <- residuals[idx] / cf[i]
    
  }
  
  ### Residuals qq plot
  
  extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
  
  if(length(extreme_id) > 0){
    
    test_resid <- residuals[-extreme_id]
    
  } else {
    
    test_resid <- residuals
    
  }
  
  ks2 <- ks.test(test_resid, 
                 pnorm, 
                 mean(test_resid), 
                 sd(test_resid))
  
  if(ks2$p.value > ks$p.value){
    
    mem_model <- mem_model2
    
  }
  
  }
  ## Extract residuals
  
  residuals <- resid(mem_model, type = "pearson")
  
  ## Correct residuals for variance structure
  
  cf <- coef(mem_model$modelStruct$varStruct, 
             unconstrained = FALSE)
  
  for(i in names(cf)){
    
    idx <- tmp_df$model == i
    
    residuals[idx] <- residuals[idx] / cf[i]
    
  }
  
  ## Plot residuals and save to file
  
  ### Set filename
  
  if(length(model_options) == 2){
    
    chapter <- "Ch2"
    
  } else {
    
    chapter <- "Ch3"
    
  }
  
  filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s.pdf",
                      names(pred_sets[prediction]),
                      "species_richness_difference",
                      chapter)
  
  ### Plot to pdf
  
  pdf(filename)
  
  ### Residuals boxplot
  
  plot(residuals ~ tmp_df$model,
       ylab = "Pearson residuals",
       xlab = "Model")
  
  ### Residuals qq plot
  
  extreme_id <- which(tmp_df$mean == 0 | tmp_df$mean == 1)
  
  if(length(extreme_id) > 0){
    
    test_resid <- residuals[-extreme_id]
    
  } else {
    
    test_resid <- residuals
    
  }
  
  qqplot(test_resid, 
         rnorm(n = 1000, 
               mean(test_resid), 
               sd(test_resid)),
         ylab = "N(mean(residuals), sd(residuals))",
         xlab = "Model Pearson residuals",
         main = "QQ-plot")
  abline(0,1)
  
  dev.off()
  
  ### Check normality with Kolmogorov-Smirnov test
  
  ks <- ks.test(test_resid, 
                pnorm, 
                mean(test_resid), 
                sd(test_resid))
  
  ks_sr[row_index, ] <- list(names(pred_sets[prediction]),
                             ts,
                             ks$p.value,
                             ks$statistic,
                             1)
  
  ### Save model to file
  
  filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                      names(pred_sets[prediction]),
                      "species_richness_difference",
                      chapter)
  
  saveRDS(mem_model,
          filename)
  
  ### Statusbar
  
  statusbar(run = row_index,
            max.run = (n_row_ts + n_row_sr + n_row_ll),
            info = sprintf("%s: %s - SR",
                           row_index,
                           names(pred_sets[prediction])))
  
  row_index <- row_index + 1
  
  suppressWarnings(rm(ks))
  suppressWarnings(rm(ks2))
}

ks_df <- rbind(ks_df,
               ks_sr)

##################
##################
### LIKELIHOOD ###
##################
##################

############################
### Mixed Effects Models ###
############################

## Loop over different different model iterations

row_index <- 1

for(prediction in seq_len(length(pred_sets))){

  ## Only run for probabilistic predictions
  
  if(prediction %nin% c(2, 6:8)){
    
    next()
    
  }
  
  ## Subset dataset
  
  tmp_df_i <- ll_i_df[ll_i_df$prediction_type %in% pred_sets[[prediction]] &
                        ll_i_df$test_statistic == "independent_log_likelihood", ]
  
  tmp_df_i <- na.omit(tmp_df_i[is.finite(tmp_df_i$mean), ])
  
  tmp_df_i <- tmp_df_i[is.finite(tmp_df_i$mean), ]
  
  tmp_df_j <- ll_j_df[ll_j_df$prediction_type %in% pred_sets[[prediction]] &
                        ll_j_df$test_statistic == "joint_log_likelihood", ]
  
  tmp_df_j <- na.omit(tmp_df_j[is.finite(tmp_df_j$mean), ])
  
  tmp_df_j <- tmp_df_j[is.finite(tmp_df_j$mean), ]
  
  ## Fit mixed effects model
  
  tmp_df_i$mean_trans <- tmp_df_i$mean
  
  tmp_df_j$mean_trans <- tmp_df_j$mean
  
  formula <- as.formula("mean_trans ~ dataset + fold + MPR + LPR + HLR_NS + HLR_S + HPR * dataset")# + DPR * dataset")
  
  mem_model_i <- tryCatch(expr = nlme::lme(formula,
                                           random = ~ 1|species,
                                           weights = varIdent(form = ~1|model),
                                           data = tmp_df_i[is.finite(tmp_df_i$mean), ],
                                           control = lmeControl(msMaxIter = 1000,
                                                                opt = "optim")),
                          error = function(err){
                            
                            message(sprintf("Model fit failed for: %s - %s",
                                            names(pred_sets[prediction]),
                                            "independent_log_likelihood"))
                            
                            return(NA)
                            
                          }
  )
  
  mem_model_j <- tryCatch(expr = nlme::lme(formula, 
                                           random = ~ 1|site,
                                           weights = varIdent(form = ~1|model),
                                           data = tmp_df_j,
                                           control = lmeControl(msMaxIter = 1000,
                                                                opt = "optim")),
                          error = function(err){
                            
                            message(sprintf("Model fit failed for: %s - %s",
                                            names(pred_sets[prediction]),
                                            "joint_log_likelihood"))
                            
                            return(NA)
                            
                          }
  )
  
  
  if(class(mem_model_i) != "lme"){
    
    ks_ll[row_index, ] <- list(names(pred_sets[prediction]),
                               "independent_log_likelihood",
                               NA,
                               NA,
                               0)
    
    row_index <- row_index + 1
    
  } else {
    
    mem_model_i$transformed <- FALSE
    
    mem_model_j$transformed <- FALSE
    
    ## Extract residuals
    
    residuals <- resid(mem_model_i, type = "pearson")
    
    ## Correct residuals for variance structure
    
    cf <- coef(mem_model_i$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df_i$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ## Plot residuals and save to file
    
    ### Set filename
    
    if(length(model_options) == 2){
      
      chapter <- "Ch2"
      
    } else {
      
      chapter <- "Ch3"
      
    }
    
    filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s.pdf",
                        names(pred_sets[prediction]),
                        "independent_log_likelihood",
                        chapter)
    
    ### Plot to pdf
    
    pdf(filename)
    
    ### Residuals boxplot
    
    plot(residuals ~ tmp_df_i$model[is.finite(tmp_df_i$mean)], 
         ylab = "Pearson residuals",
         xlab = "Model")
    
    ### Residuals qq plot
    
    extreme_id <- which(tmp_df_i$mean == 0 | tmp_df_i$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    qqplot(test_resid, 
           rnorm(n = 1000, 
                 mean(test_resid), 
                 sd(test_resid)),
           ylab = "N(mean(residuals), sd(residuals))",
           xlab = "Model Pearson residuals",
           main = "QQ-plot")
    abline(0,1)
    
    dev.off()
    
    ### Check normality with Kolmogorov-Smirnov test
    
    ks <- ks.test(test_resid, 
                  pnorm, 
                  mean(test_resid), 
                  sd(test_resid))
    
    ks_ll[row_index, ] <- list(names(pred_sets[prediction]),
                               "independent_log_likelihood",
                               ks$p.value,
                               ks$statistic,
                               1)
    
    ### Save model to file
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        names(pred_sets[prediction]),
                        "independent_log_likelihood",
                        chapter)
    
    saveRDS(mem_model_i,
            filename)
    
    ### Statusbar
    
    statusbar(run = row_index,
              max.run = (n_row_ts + n_row_sr + n_row_ll),
              info = sprintf("%s: %s - LL_I",
                             row_index,
                             names(pred_sets[prediction])))
    
    row_index <- row_index + 1
    
  }
  
  if(class(mem_model_j) != "lme"){
    
    ks_ll[row_index, ] <- list(names(pred_sets[prediction]),
                               "joint_log_likelihood",
                               NA,
                               NA,
                               0)
    
    row_index <- row_index + 1
    
  } else {
    
    ## Extract residuals
    
    residuals <- resid(mem_model_j, type = "pearson")
    
    ## Correct residuals for variance structure
    
    cf <- coef(mem_model_j$modelStruct$varStruct, 
               unconstrained = FALSE)
    
    for(i in names(cf)){
      
      idx <- tmp_df_j$model == i
      
      residuals[idx] <- residuals[idx] / cf[i]
      
    }
    
    ## Plot residuals and save to file
    
    ### Set filename
    
    if(length(model_options) == 2){
      
      chapter <- "Ch2"
      
    } else {
      
      chapter <- "Ch3"
      
    }
    
    filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s.pdf",
                        names(pred_sets[prediction]),
                        "joint_log_likelihood",
                        chapter)
    
    ### Plot to pdf
    
    pdf(filename)
    
    ### Residuals boxplot
    
    plot(residuals ~ tmp_df_j$model, 
         ylab = "Pearson residuals",
         xlab = "Model")
    
    ### Residuals qq plot
    
    extreme_id <- which(tmp_df_j$mean == 0 | tmp_df_j$mean == 1)
    
    if(length(extreme_id) > 0){
      
      test_resid <- residuals[-extreme_id]
      
    } else {
      
      test_resid <- residuals
      
    }
    
    qqplot(test_resid, 
           rnorm(n = 1000, 
                 mean(test_resid), 
                 sd(test_resid)),
           ylab = "N(mean(residuals), sd(residuals))",
           xlab = "Model Pearson residuals",
           main = "QQ-plot")
    abline(0,1)
    
    dev.off()
    
    ### Check normality with Kolmogorov-Smirnov test
    
    ks <- ks.test(test_resid, 
                  pnorm, 
                  mean(test_resid), 
                  sd(test_resid))
    
    ks_ll[row_index, ] <- list(names(pred_sets[prediction]),
                               "joint_log_likelihood",
                               ks$p.value,
                               ks$statistic,
                               1)
    
    ### Save model to file
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        names(pred_sets[prediction]),
                        "joint_log_likelihood",
                        chapter)
    
    saveRDS(mem_model_j,
            filename)
    
    ### Statusbar
    
    statusbar(run = row_index,
              max.run = (n_row_ts + n_row_sr + n_row_ll),
              info = sprintf("%s: %s - LL_J",
                             row_index,
                             names(pred_sets[prediction])))
    
    row_index <- row_index + 1
    
  }
}

ks_df <- rbind(ks_df,
               ks_ll)

#######################################
#######################################
### KOLMOGOROV-SMIRNOV TEST SUMMARY ###
#######################################
#######################################

## Drop empty rows from incompatible ts/pred combos

ks_df <- ks_df[ks_df$test_statistic != "", ]

## Test to see if any are now still > 0.05

ks_fail <- ks_df[ks_df$p_value < bonferroni_p, ]

# ##################
# ##################
# ### MAKE PLOTS ###
# ##################
# ##################
# 
# dodge <- position_dodge(width = 0.5)
# 
# #########################
# ### Set colour scheme ###
# #########################
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
#                    condLOI_marg_low = c("condLOI_marg_low", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Conditional Marginal (low)"),
#                    condLOI_marg_med = c("condLOI_marg_med", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Conditional Marginal (med)"),
#                    condLOI_marg_high = c("condLOI_marg_high", "SSDM_bin", "SESAM", "Probabilistic predictions", "JSDM prediction: Conditional Marginal (high)"),
#                    joint = c("joint", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Joint"))
# 
# ################################################
# ### Generate plot for different combinations ###
# ################################################
# 
# for(prediction in seq_len(length(pred_sets))){
#   
#   for(ts in c(ts_species, 
#               ts_site, 
#               "species_richness_difference",
#               "independent_log_likelihood",
#               "joint_log_likelihood")){
#     
#     ## Check ts/prediction compatibility
#     
#     ### Binary test statistics for binary predictions
#     
#     if(names(pred_sets)[prediction] %in% binary_predictions & ts %nin% binary_ts){
#       
#       next()
#       
#     }
#     
#     ### Probabilistic test statistics for probabilistic predictions
#     
#     if(names(pred_sets)[prediction] %in% probability_predictions & ts %nin% prob_ts){
#       
#       next()
#       
#     }
#     
#     ## Read in file
#     
#     filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
#                         names(pred_sets[prediction]),
#                         ts,
#                         chapter)
#     
#     if(!file.exists(filename)){
#       
#       message(sprintf("No model for: %s - %s",
#                       names(pred_sets[prediction]),
#                       ts))
#       
#       next()
#       
#     }
#     
#     mem_model <- readRDS(filename)
#     
#     ## Extract values
#     
#     model_summary <- summary(mem_model)
#     
#     coefs <- model_summary$tTable
#     
#     ## Dataframe to store values
#     
#     n_row <- length(model_order)
#     
#     plot_df <- data.frame(model = factor(model_order,
#                                          levels = model_order),
#                           mean = numeric(n_row),
#                           upper = numeric(n_row),
#                           lower = numeric(n_row))
#     
#     for(i in seq_len(length(model_order))){
#       
#       model <- model_order[i]
#       
#       if(paste0("model", model) %nin% rownames(coefs)){
#         
#         plot_df[i, ] <- list(model,
#                              NA,
#                              NA,
#                              NA)
#       }
#       
#       if(paste0("model", model) %in% rownames(coefs)){
#         
#         # plot_df[i, ] <- list(model,                                         # Model
#         #                      coefs[paste0("model", model), "Value"],        # Mean
#         #                      coefs[paste0("model", model), "Value"] +       # Upper
#         #                        qnorm(0.975) * coefs[paste0("model", model), "Std.Error"],
#         #                      coefs[paste0("model", model), "Value"] +       # Lower
#         #                        qnorm(0.025) * coefs[paste0("model", model), "Std.Error"])
#         
#         pred_coef <- MEM_predict(mem = mem_model,
#                                  model = model)
#         
#         mn <- pred_coef$mean
#         se <- pred_coef$se
#         
#         if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
#           
#           quantile_fun <- qnorm
#           
#           mean <- mn
#           
#         }
#         
#         if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
#           
#           quantile_fun <- qlnorm
#           
#           mean <- exp(mn + (se ^ 2 / 2))
#           
#         }
#         
#         if(mem_model$transformed == TRUE & ts %in% ts_0_1){
#           
#           quantile_fun <- qlogitnorm
#           
#           mean <- mean(plogis(rnorm(1000000, mn, se)))
#           
#         }
#         
#         plot_df[i, ] <- list(model,
#                              mean,
#                              quantile_fun(0.975, mn, se),
#                              quantile_fun(0.025, mn, se))
#         
#       }
#     }
#     
#     ## Make plot and save to file
#     
#     filename <- sprintf("outputs/test_statistics/plots/%1$s_%2$s_%3$s_MME.pdf",
#                         names(pred_sets[prediction]),
#                         ts,
#                         chapter)
#     pdf(filename)
#     
#     tmp_plot <- tryCatch(expr = ggplot(plot_df,
#                                        aes(x = model,
#                                            y = mean,
#                                            colour = model)) +
#                            geom_hline(yintercept = plot_df[plot_df$model == "SSDM", "mean"],
#                                       linetype = 11) +
#                            geom_hline(yintercept = plot_df[plot_df$model == "SSDM", "upper"],
#                                       linetype = 11,
#                                       colour = "lightgrey") +
#                            geom_hline(yintercept = plot_df[plot_df$model == "SSDM", "lower"],
#                                       linetype = 11,
#                                       colour = "lightgrey") + 
#                            geom_point(position = dodge,
#                                       size = 2) + 
#                            geom_errorbar(aes(ymax = upper,
#                                              ymin = lower,),
#                                          position = dodge,
#                                          width = 0.1) +
#                            ylim(c(graph_customisation[graph_customisation$test_statistic == ts, "y_min"],
#                                   graph_customisation[graph_customisation$test_statistic == ts, "y_max"])) +
#                            xlab("Model") +
#                            ylab(eval(parse(text = graph_customisation[graph_customisation$test_statistic == ts, "graph_label"]))) +
#                            scale_colour_manual(values = colour,
#                                                breaks = levels(plot_df$model)) +
#                            theme_bw() +
#                            theme(legend.position = "none",
#                                  panel.grid.minor = element_blank(),
#                                  panel.grid.major = element_blank()) +
#                            ggtitle(label = pred_pairs[[prediction]][4],
#                                    subtitle = pred_pairs[[prediction]][5]),
#                          error = function(err){
#                            
#                            return(NA)
#                            
#                          })
#     
#     if(chapter == "Ch3"){
#       
#       tmp_plot <- tmp_plot + geom_vline(xintercept = 2.5,
#                                         lwd = 1.1)
#       
#     }
#     
#     tryCatch(expr = print(tmp_plot),
#              error = function(err){
#                
#                message(sprintf("Plot failed for: %s - %s",
#                                names(pred_sets)[prediction],
#                                ts))
#                
#              }
#     )
#     
#     dev.off()
#   }
# }
