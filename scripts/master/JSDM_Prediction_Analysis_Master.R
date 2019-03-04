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

##################################
### Command line arguments and ###
###       defining indices     ###
##################################

## Read in the command line arguments

command_args <- commandArgs(trailingOnly = TRUE)

dataset_index <- as.numeric(command_args[1])

## Possible ID options

dataset_options <- c("frog",
                     "eucalypt",
                     "bird")

## Extract IDs corresponding to command line argument indices

dataset_id <- dataset_options[dataset_index]

message("Command line arguments finished")

#####################
### Load packages ###
#####################

.libPaths("/home/davidpw/R/lib/3.5")

library(abind)

message("Packages loaded")

#############################
### Load Function Scripts ###
#############################

source("scripts/analysis/proportionNA_function.R")

########################################
### Summarise Test Statistic Metrics ###
########################################

## Model names

JSDM_model_names <- c("MPR", "HPR", "LPR",
                      "DPR", "HLR_NS", "HLR_S")

SSDM_model_names <- c("SSDM", "SESAM")

model_names <- c(JSDM_model_names, SSDM_model_names)

## Test statistic names

ts_spp_names <- c("AUC", "bias", "MSE", "R2", "RMSE", "SSE", "Pearson",
                  "Spearman", "Kendall", "TP", "FP", "TN", "FN", "TPR",
                  "FPR", "TNR", "FNR", "PLR", "NLR", "DOR", "Prevalence",
                  "Accuracy", "PPV", "FOR", "FDR", "NPV", "F_1",
                  "Youden_J", "Kappa")

ts_site_names <- c("Binomial", "Bray", "Canberra",
                   "Euclidean", "Gower", "Gower_alt", "Horn", "Jaccard",    
                   "Kulczynski", "Mahalanobis", "Manhattan", "Mountford",
                   "Raup")

ts_all_names <- c(ts_spp_names, ts_site_names)

###
### BINARY PREDICTIONS
###

## Marginal Prediction ----

### Create empty storage

ts_summary_marg_bin <- list(NA_proportion = matrix(data = NA,
                                                   nrow = length(ts_all_names),
                                                   ncol = length(JSDM_model_names)),
                            mean = matrix(data = NA,
                                          nrow = length(ts_all_names),
                                          ncol = length(JSDM_model_names)),
                            median = matrix(data = NA,
                                            nrow = length(ts_all_names),
                                            ncol = length(JSDM_model_names)),
                            CI_95_lower = matrix(data = NA,
                                                 nrow = length(ts_all_names),
                                                 ncol = length(JSDM_model_names)),
                            CI_95_upper = matrix(data = NA,
                                                 nrow = length(ts_all_names),
                                                 ncol = length(JSDM_model_names)))

for(i in seq_len(length(ts_summary_marg_bin))){
  
  rownames(ts_summary_marg_bin[[i]]) <- ts_all_names
  
  colnames(ts_summary_marg_bin[[i]]) <- JSDM_model_names

}

### Fill values

for(model in JSDM_model_names){
  
  ## Load test statistics for all folds
  
  ts_1 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold1_marginal_bin_ts.rds",
                          model,
                          dataset_id))
  
  ts_2 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold2_marginal_bin_ts.rds",
                          model,
                          dataset_id))
  
  ts_3 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold3_marginal_bin_ts.rds",
                          model,
                          dataset_id))
  
  ts_4 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold4_marginal_bin_ts.rds",
                          model,
                          dataset_id))
  
  ts_5 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold5_marginal_bin_ts.rds",
                          model,
                          dataset_id))
  
  ## Merge folds
  
  ts_all_spp <- abind(ts_1$test_statistics_species,
                      ts_2$test_statistics_species,
                      ts_3$test_statistics_species,
                      ts_4$test_statistics_species,
                      ts_5$test_statistics_species,
                      along = 1)
  
  ts_all_site <- abind(ts_1$test_statistics_site,
                       ts_2$test_statistics_site,
                       ts_3$test_statistics_site,
                       ts_4$test_statistics_site,
                       ts_5$test_statistics_site,
                       along = 1)
  
  ## NA_percentage
  
  for(ts in ts_spp_names){
    
    ts_summary_marg_bin$NA_proportion[ts, model] <- proportionNA(ts_all_spp[ , ts, ])
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_marg_bin$NA_proportion[ts, model] <- proportionNA(ts_all_site[ , ts, ])
    
  }
  
  ## Mean
  
  for(ts in ts_spp_names){
    
    ts_summary_marg_bin$mean[ts, model] <- mean(ts_all_spp[ , ts, ],
                                            na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_marg_bin$mean[ts, model] <- mean(ts_all_site[ , ts, ],
                                            na.rm = TRUE)
    
  }
  
  ## Median
  
  for(ts in ts_spp_names){
    
    ts_summary_marg_bin$median[ts, model] <- median(ts_all_spp[ , ts, ],
                                                na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_marg_bin$median[ts, model] <- median(ts_all_site[ , ts, ],
                                                na.rm = TRUE)
    
  }
  
  ## CI_95_lower
  
  for(ts in ts_spp_names){
    
    ts_summary_marg_bin$CI_95_lower[ts, model] <- quantile(ts_all_spp[ , ts, ],
                                                       probs = 0.025,
                                                       na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_marg_bin$CI_95_lower[ts, model] <- quantile(ts_all_site[ , ts, ],
                                                       probs = 0.025,
                                                       na.rm = TRUE)[[1]]
    
  }
  
  ## CI_95_upper
  
  for(ts in ts_spp_names){
    
    ts_summary_marg_bin$CI_95_upper[ts, model] <- quantile(ts_all_spp[ , ts, ],
                                                       probs = 0.975,
                                                       na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_marg_bin$CI_95_upper[ts, model] <- quantile(ts_all_site[ , ts, ],
                                                       probs = 0.975,
                                                       na.rm = TRUE)[[1]]
    
  }
  
}

## Save to file

filename <- sprintf("outputs/test_statistics/ts_summary_JSDM_%s_marginal_bin.rds",
                    dataset_id)

saveRDS(ts_summary_marg_bin,
        filename)

#----

## Conditional LOI Prediction ----

### Create empty storage

ts_summary_condLOI_low <- list(NA_proportion = matrix(data = NA,
                                                      nrow = length(ts_all_names),
                                                      ncol = length(JSDM_model_names)),
                               mean = matrix(data = NA,
                                             nrow = length(ts_all_names),
                                             ncol = length(JSDM_model_names)),
                               median = matrix(data = NA,
                                               nrow = length(ts_all_names),
                                               ncol = length(JSDM_model_names)),
                               CI_95_lower = matrix(data = NA,
                                                    nrow = length(ts_all_names),
                                                    ncol = length(JSDM_model_names)),
                               CI_95_upper = matrix(data = NA,
                                                    nrow = length(ts_all_names),
                                                    ncol = length(JSDM_model_names)))

ts_summary_condLOI_med <- list(NA_proportion = matrix(data = NA,
                                                      nrow = length(ts_all_names),
                                                      ncol = length(JSDM_model_names)),
                               mean = matrix(data = NA,
                                             nrow = length(ts_all_names),
                                             ncol = length(JSDM_model_names)),
                               median = matrix(data = NA,
                                               nrow = length(ts_all_names),
                                               ncol = length(JSDM_model_names)),
                               CI_95_lower = matrix(data = NA,
                                                    nrow = length(ts_all_names),
                                                    ncol = length(JSDM_model_names)),
                               CI_95_upper = matrix(data = NA,
                                                    nrow = length(ts_all_names),
                                                    ncol = length(JSDM_model_names)))

ts_summary_condLOI_high <- list(NA_proportion = matrix(data = NA,
                                                       nrow = length(ts_all_names),
                                                       ncol = length(JSDM_model_names)),
                                mean = matrix(data = NA,
                                              nrow = length(ts_all_names),
                                              ncol = length(JSDM_model_names)),
                                median = matrix(data = NA,
                                                nrow = length(ts_all_names),
                                                ncol = length(JSDM_model_names)),
                                CI_95_lower = matrix(data = NA,
                                                     nrow = length(ts_all_names),
                                                     ncol = length(JSDM_model_names)),
                                CI_95_upper = matrix(data = NA,
                                                     nrow = length(ts_all_names),
                                                     ncol = length(JSDM_model_names)))

for(i in seq_len(length(ts_summary_condLOI_low))){
  
  rownames(ts_summary_condLOI_low[[i]]) <- ts_all_names
  
  colnames(ts_summary_condLOI_low[[i]]) <- JSDM_model_names
  
}

for(i in seq_len(length(ts_summary_condLOI_med))){
  
  rownames(ts_summary_condLOI_med[[i]]) <- ts_all_names
  
  colnames(ts_summary_condLOI_med[[i]]) <- JSDM_model_names
  
}

for(i in seq_len(length(ts_summary_condLOI_high))){
  
  rownames(ts_summary_condLOI_high[[i]]) <- ts_all_names
  
  colnames(ts_summary_condLOI_high[[i]]) <- JSDM_model_names
  
}

### Fill values

for(model in JSDM_model_names){
  
  ## Load test statistics for all folds
  
  ts_1 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold1_condLOI_ts.rds",
                          model,
                          dataset_id))
  
  ts_2 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold2_condLOI_ts.rds",
                          model,
                          dataset_id))
  
  ts_3 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold3_condLOI_ts.rds",
                          model,
                          dataset_id))
  
  ts_4 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold4_condLOI_ts.rds",
                          model,
                          dataset_id))
  
  ts_5 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold5_condLOI_ts.rds",
                          model,
                          dataset_id))
  
  ## Merge folds
  
  ts_all_spp_low <- abind(ts_1[[1]]$test_statistics_species,
                          ts_2[[1]]$test_statistics_species,
                          ts_3[[1]]$test_statistics_species,
                          ts_4[[1]]$test_statistics_species,
                          ts_5[[1]]$test_statistics_species,
                          along = 1)
  
  ts_all_site_low <- abind(ts_1[[1]]$test_statistics_site,
                           ts_2[[1]]$test_statistics_site,
                           ts_3[[1]]$test_statistics_site,
                           ts_4[[1]]$test_statistics_site,
                           ts_5[[1]]$test_statistics_site,
                           along = 1)
  
  ts_all_spp_med <- abind(ts_1[[2]]$test_statistics_species,
                          ts_2[[2]]$test_statistics_species,
                          ts_3[[2]]$test_statistics_species,
                          ts_4[[2]]$test_statistics_species,
                          ts_5[[2]]$test_statistics_species,
                          along = 1)
  
  ts_all_site_med <- abind(ts_1[[2]]$test_statistics_site,
                           ts_2[[2]]$test_statistics_site,
                           ts_3[[2]]$test_statistics_site,
                           ts_4[[2]]$test_statistics_site,
                           ts_5[[2]]$test_statistics_site,
                           along = 1)
  
  ts_all_spp_high <- abind(ts_1[[3]]$test_statistics_species,
                           ts_2[[3]]$test_statistics_species,
                           ts_3[[3]]$test_statistics_species,
                           ts_4[[3]]$test_statistics_species,
                           ts_5[[3]]$test_statistics_species,
                           along = 1)
  
  ts_all_site_high <- abind(ts_1[[3]]$test_statistics_site,
                            ts_2[[3]]$test_statistics_site,
                            ts_3[[3]]$test_statistics_site,
                            ts_4[[3]]$test_statistics_site,
                            ts_5[[3]]$test_statistics_site,
                            along = 1)
  
  ## NA_percentage
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_low$NA_proportion[ts, model] <- proportionNA(ts_all_spp_low[ , ts, ])
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_low$NA_proportion[ts, model] <- proportionNA(ts_all_site_low[ , ts, ])
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_med$NA_proportion[ts, model] <- proportionNA(ts_all_spp_med[ , ts, ])
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_med$NA_proportion[ts, model] <- proportionNA(ts_all_site_med[ , ts, ])
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_high$NA_proportion[ts, model] <- proportionNA(ts_all_spp_high[ , ts, ])
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_high$NA_proportion[ts, model] <- proportionNA(ts_all_site_high[ , ts, ])
    
  }
  
  ## Mean
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_low$mean[ts, model] <- mean(ts_all_spp_low[ , ts, ],
                                               na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_low$mean[ts, model] <- mean(ts_all_site_low[ , ts, ],
                                               na.rm = TRUE)
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_med$mean[ts, model] <- mean(ts_all_spp_med[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_med$mean[ts, model] <- mean(ts_all_site_med[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_high$mean[ts, model] <- mean(ts_all_spp_high[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_high$mean[ts, model] <- mean(ts_all_site_high[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  ## Median
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_low$median[ts, model] <- median(ts_all_spp_low[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_low$median[ts, model] <- median(ts_all_site_low[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_med$median[ts, model] <- median(ts_all_spp_med[ , ts, ],
                                                       na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_med$median[ts, model] <- median(ts_all_site_med[ , ts, ],
                                                       na.rm = TRUE)
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_high$median[ts, model] <- median(ts_all_spp_high[ , ts, ],
                                                       na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_high$median[ts, model] <- median(ts_all_site_high[ , ts, ],
                                                       na.rm = TRUE)
    
  }
  
  ## CI_95_lower
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_low$CI_95_lower[ts, model] <- quantile(ts_all_spp_low[ , ts, ],
                                                          probs = 0.025,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_low$CI_95_lower[ts, model] <- quantile(ts_all_site_low[ , ts, ],
                                                          probs = 0.025,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_med$CI_95_lower[ts, model] <- quantile(ts_all_spp_med[ , ts, ],
                                                              probs = 0.025,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_med$CI_95_lower[ts, model] <- quantile(ts_all_site_med[ , ts, ],
                                                              probs = 0.025,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_high$CI_95_lower[ts, model] <- quantile(ts_all_spp_high[ , ts, ],
                                                              probs = 0.025,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_high$CI_95_lower[ts, model] <- quantile(ts_all_site_high[ , ts, ],
                                                              probs = 0.025,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  ## CI_95_upper
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_low$CI_95_upper[ts, model] <- quantile(ts_all_spp_low[ , ts, ],
                                                          probs = 0.975,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_low$CI_95_upper[ts, model] <- quantile(ts_all_site_low[ , ts, ],
                                                          probs = 0.975,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_med$CI_95_upper[ts, model] <- quantile(ts_all_spp_med[ , ts, ],
                                                              probs = 0.975,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_med$CI_95_upper[ts, model] <- quantile(ts_all_site_med[ , ts, ],
                                                              probs = 0.975,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_spp_names){
    
    ts_summary_condLOI_high$CI_95_upper[ts, model] <- quantile(ts_all_spp_high[ , ts, ],
                                                              probs = 0.975,
                                                              na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_condLOI_high$CI_95_upper[ts, model] <- quantile(ts_all_site_high[ , ts, ],
                                                              probs = 0.975,
                                                              na.rm = TRUE)[[1]]
    
  }
  
}

## Save to file

filename <- sprintf("outputs/test_statistics/ts_summary_JSDM_%s_condLOI_low.rds",
                    dataset_id)

saveRDS(ts_summary_condLOI_low,
        filename)

filename <- sprintf("outputs/test_statistics/ts_summary_JSDM_%s_condLOI_med.rds",
                    dataset_id)

saveRDS(ts_summary_condLOI_med,
        filename)

filename <- sprintf("outputs/test_statistics/ts_summary_JSDM_%s_condLOI_high.rds",
                    dataset_id)

saveRDS(ts_summary_condLOI_high,
        filename)

#----

# ## Conditional LOO Prediction ----
# 
# ### Create empty storage
# 
# ts_summary_condLOO <- list(NA_proportion = matrix(data = NA,    # Placeholder
#                                                   nrow = 42,    # Number of test statistics
#                                                   ncol = 6),    # Number of JSDM models
#                            mean = matrix(data = NA,
#                                          nrow = 42,
#                                          ncol = 6),
#                            median = matrix(data = NA,
#                                            nrow = 42,
#                                            ncol = 6),
#                            CI_95_lower = matrix(data = NA,
#                                                 nrow = 42,
#                                                 ncol = 6),
#                            CI_95_upper = matrix(data = NA,
#                                                 nrow = 42,
#                                                 ncol = 6))
# 
# for(i in seq_len(length(ts_summary_condLOO))){
#   
#   rownames(ts_summary_condLOO[[i]]) <- ts_all_names
#   
#   colnames(ts_summary_condLOO[[i]]) <- model_names
#   
# }
# 
# ### Fill values
# 
# for(model in model_names){
#   
#   ## Load test statistics for all folds
#   
#   ts_1 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold1_condLOO_ts.rds",
#                           model,
#                           dataset_id))
#   
#   ts_2 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold2_condLOO_ts.rds",
#                           model,
#                           dataset_id))
#   
#   ts_3 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold3_condLOO_ts.rds",
#                           model,
#                           dataset_id))
#   
#   ts_4 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold4_condLOO_ts.rds",
#                           model,
#                           dataset_id))
#   
#   ts_5 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold5_condLOO_ts.rds",
#                           model,
#                           dataset_id))
#   
#   ## Merge folds
#   
#   ts_all_spp <- abind(ts_1$test_statistics_species,
#                       ts_2$test_statistics_species,
#                       ts_3$test_statistics_species,
#                       ts_4$test_statistics_species,
#                       ts_5$test_statistics_species,
#                       along = 1)
#   
#   ts_all_site <- abind(ts_1$test_statistics_site,
#                        ts_2$test_statistics_site,
#                        ts_3$test_statistics_site,
#                        ts_4$test_statistics_site,
#                        ts_5$test_statistics_site,
#                        along = 1)
#   
#   ## NA_percentage
#   
#   for(ts in ts_spp_names){
#     
#     ts_summary_condLOO$NA_proportion[ts, model] <- proportionNA(ts_all_spp[ , ts, ])
#     
#   }
#   
#   for(ts in ts_site_names){
#     
#     ts_summary_condLOO$NA_proportion[ts, model] <- proportionNA(ts_all_site[ , ts, ])
#     
#   }
#   
#   ## Mean
#   
#   for(ts in ts_spp_names){
#     
#     ts_summary_condLOO$mean[ts, model] <- mean(ts_all_spp[ , ts, ],
#                                                na.rm = TRUE)
#     
#   }
#   
#   for(ts in ts_site_names){
#     
#     ts_summary_condLOO$mean[ts, model] <- mean(ts_all_site[ , ts, ],
#                                                na.rm = TRUE)
#     
#   }
#   
#   ## Median
#   
#   for(ts in ts_spp_names){
#     
#     ts_summary_condLOO$median[ts, model] <- median(ts_all_spp[ , ts, ],
#                                                    na.rm = TRUE)
#     
#   }
#   
#   for(ts in ts_site_names){
#     
#     ts_summary_condLOO$median[ts, model] <- median(ts_all_site[ , ts, ],
#                                                    na.rm = TRUE)
#     
#   }
#   
#   ## CI_95_lower
#   
#   for(ts in ts_spp_names){
#     
#     ts_summary_condLOO$CI_95_lower[ts, model] <- quantile(ts_all_spp[ , ts, ],
#                                                           probs = 0.025,
#                                                           na.rm = TRUE)[[1]]
#     
#   }
#   
#   for(ts in ts_site_names){
#     
#     ts_summary_condLOO$CI_95_lower[ts, model] <- quantile(ts_all_site[ , ts, ],
#                                                           probs = 0.025,
#                                                           na.rm = TRUE)[[1]]
#     
#   }
#   
#   ## CI_95_upper
#   
#   for(ts in ts_spp_names){
#     
#     ts_summary_condLOO$CI_95_upper[ts, model] <- quantile(ts_all_spp[ , ts, ],
#                                                           probs = 0.975,
#                                                           na.rm = TRUE)[[1]]
#     
#   }
#   
#   for(ts in ts_site_names){
#     
#     ts_summary_condLOO$CI_95_upper[ts, model] <- quantile(ts_all_site[ , ts, ],
#                                                           probs = 0.975,
#                                                           na.rm = TRUE)[[1]]
#     
#   }
#   
#   
# }
# 
# ## Save to file
# 
# filename <- sprintf("outputs/test_statistics/ts_summary_JSDM_%s_condLOO.rds",
#                     dataset_id)
# 
# saveRDS(ts_summary_condLOO,
#         filename)
# 
# #----

## Joint Prediction ----

### Create empty storage

ts_summary_joint <- list(NA_proportion = matrix(data = NA,
                                                nrow = length(ts_all_names),
                                                ncol = length(JSDM_model_names)),
                           mean = matrix(data = NA,
                                         nrow = length(ts_all_names),
                                         ncol = length(JSDM_model_names)),
                           median = matrix(data = NA,
                                           nrow = length(ts_all_names),
                                           ncol = length(JSDM_model_names)),
                           CI_95_lower = matrix(data = NA,
                                                nrow = length(ts_all_names),
                                                ncol = length(JSDM_model_names)),
                           CI_95_upper = matrix(data = NA,
                                                nrow = length(ts_all_names),
                                                ncol = length(JSDM_model_names)))

for(i in seq_len(length(ts_summary_joint))){
  
  rownames(ts_summary_joint[[i]]) <- ts_all_names
  
  colnames(ts_summary_joint[[i]]) <- JSDM_model_names
  
}

### Fill values

for(model in model_names){
  
  ## Load test statistics for all folds
  
  ts_1 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold1_joint_ts.rds",
                          model,
                          dataset_id))
  
  ts_2 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold2_joint_ts.rds",
                          model,
                          dataset_id))
  
  ts_3 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold3_joint_ts.rds",
                          model,
                          dataset_id))
  
  ts_4 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold4_joint_ts.rds",
                          model,
                          dataset_id))
  
  ts_5 <- readRDS(sprintf("outputs/test_statistics/%s_%s_fold5_joint_ts.rds",
                          model,
                          dataset_id))
  
  ## Merge folds
  
  ts_all_spp <- abind(ts_1$test_statistics_species,
                      ts_2$test_statistics_species,
                      ts_3$test_statistics_species,
                      ts_4$test_statistics_species,
                      ts_5$test_statistics_species,
                      along = 1)
  
  ts_all_site <- abind(ts_1$test_statistics_site,
                       ts_2$test_statistics_site,
                       ts_3$test_statistics_site,
                       ts_4$test_statistics_site,
                       ts_5$test_statistics_site,
                       along = 1)
  
  ## NA_percentage
  
  for(ts in ts_spp_names){
    
    ts_summary_joint$NA_proportion[ts, model] <- proportionNA(ts_all_spp[ , ts, ])
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_joint$NA_proportion[ts, model] <- proportionNA(ts_all_site[ , ts, ])
    
  }
  
  ## Mean
  
  for(ts in ts_spp_names){
    
    ts_summary_joint$mean[ts, model] <- mean(ts_all_spp[ , ts, ],
                                               na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_joint$mean[ts, model] <- mean(ts_all_site[ , ts, ],
                                               na.rm = TRUE)
    
  }
  
  ## Median
  
  for(ts in ts_spp_names){
    
    ts_summary_joint$median[ts, model] <- median(ts_all_spp[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_joint$median[ts, model] <- median(ts_all_site[ , ts, ],
                                                   na.rm = TRUE)
    
  }
  
  ## CI_95_lower
  
  for(ts in ts_spp_names){
    
    ts_summary_joint$CI_95_lower[ts, model] <- quantile(ts_all_spp[ , ts, ],
                                                          probs = 0.025,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_joint$CI_95_lower[ts, model] <- quantile(ts_all_site[ , ts, ],
                                                          probs = 0.025,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  ## CI_95_upper
  
  for(ts in ts_spp_names){
    
    ts_summary_joint$CI_95_upper[ts, model] <- quantile(ts_all_spp[ , ts, ],
                                                          probs = 0.975,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_joint$CI_95_upper[ts, model] <- quantile(ts_all_site[ , ts, ],
                                                          probs = 0.975,
                                                          na.rm = TRUE)[[1]]
    
  }
  
  
}

## Save to file

filename <- sprintf("outputs/test_statistics/ts_summary_JSDM_%s_joint.rds",
                    dataset_id)

saveRDS(ts_summary_joint,
        filename)

#----

## SSDM/SESAM Prediction ----

### Create empty storage

ts_summary_SSDM_bin <- list(NA_proportion = matrix(data = NA,
                                                   nrow = length(ts_all_names),
                                                   ncol = length(SSDM_model_names)),
                            mean = matrix(data = NA,
                                          nrow = length(ts_all_names),
                                          ncol = length(SSDM_model_names)),
                            median = matrix(data = NA,
                                            nrow = length(ts_all_names),
                                            ncol = length(SSDM_model_names)),
                            CI_95_lower = matrix(data = NA,
                                                 nrow = length(ts_all_names),
                                                 ncol = length(SSDM_model_names)),
                            CI_95_upper = matrix(data = NA,
                                                 nrow = length(ts_all_names),
                                                 ncol = length(SSDM_model_names)))

for(i in seq_len(length(ts_summary_SSDM_bin))){
  
  rownames(ts_summary_SSDM_bin[[i]]) <- ts_all_names
  
  colnames(ts_summary_SSDM_bin[[i]]) <- SSDM_model_names
  
}

### Fill values

for(model in SSDM_model_names){
  
  ## Load test statistics for all folds
  
  if(model == "SSDM"){
    
    ts_1 <- readRDS(sprintf("outputs/test_statistics/SSDM_bin_%s_fold1_ts.rds",
                            dataset_id))
    
    ts_2 <- readRDS(sprintf("outputs/test_statistics/SSDM_bin_%s_fold2_ts.rds",
                            dataset_id))
    
    ts_3 <- readRDS(sprintf("outputs/test_statistics/SSDM_bin_%s_fold3_ts.rds",
                            dataset_id))
    
    ts_4 <- readRDS(sprintf("outputs/test_statistics/SSDM_bin_%s_fold4_ts.rds",
                            dataset_id))
    
    ts_5 <- readRDS(sprintf("outputs/test_statistics/SSDM_bin_%s_fold5_ts.rds",
                            dataset_id))
    
  }
  
  if(model == "SESAM"){
    
    ts_1 <- readRDS(sprintf("outputs/test_statistics/SESAM_%s_fold1_ts.rds",
                            dataset_id))
    
    ts_2 <- readRDS(sprintf("outputs/test_statistics/SESAM_%s_fold2_ts.rds",
                            dataset_id))
    
    ts_3 <- readRDS(sprintf("outputs/test_statistics/SESAM_%s_fold3_ts.rds",
                            dataset_id))
    
    ts_4 <- readRDS(sprintf("outputs/test_statistics/SESAM_%s_fold4_ts.rds",
                            dataset_id))
    
    ts_5 <- readRDS(sprintf("outputs/test_statistics/SESAM_%s_fold5_ts.rds",
                            dataset_id))
    
  }
  
  ## Merge folds
  
  ts_all_spp <- abind(ts_1$test_statistics_species,
                      ts_2$test_statistics_species,
                      ts_3$test_statistics_species,
                      ts_4$test_statistics_species,
                      ts_5$test_statistics_species,
                      along = 1)
  
  ts_all_site <- abind(ts_1$test_statistics_site,
                       ts_2$test_statistics_site,
                       ts_3$test_statistics_site,
                       ts_4$test_statistics_site,
                       ts_5$test_statistics_site,
                       along = 1)
  
  ## NA_percentage
  
  for(ts in ts_spp_names){
    
    ts_summary_SSDM_bin$NA_proportion[ts, model] <- proportionNA(ts_all_spp[ , ts, ])
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_SSDM_bin$NA_proportion[ts, model] <- proportionNA(ts_all_site[ , ts, ])
    
  }
  
  ## Mean
  
  for(ts in ts_spp_names){
    
    ts_summary_SSDM_bin$mean[ts, model] <- mean(ts_all_spp[ , ts, ],
                                                na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_SSDM_bin$mean[ts, model] <- mean(ts_all_site[ , ts, ],
                                                na.rm = TRUE)
    
  }
  
  ## Median
  
  for(ts in ts_spp_names){
    
    ts_summary_SSDM_bin$median[ts, model] <- median(ts_all_spp[ , ts, ],
                                                    na.rm = TRUE)
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_SSDM_bin$median[ts, model] <- median(ts_all_site[ , ts, ],
                                                    na.rm = TRUE)
    
  }
  
  ## CI_95_lower
  
  for(ts in ts_spp_names){
    
    ts_summary_SSDM_bin$CI_95_lower[ts, model] <- quantile(ts_all_spp[ , ts, ],
                                                           probs = 0.025,
                                                           na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_SSDM_bin$CI_95_lower[ts, model] <- quantile(ts_all_site[ , ts, ],
                                                           probs = 0.025,
                                                           na.rm = TRUE)[[1]]
    
  }
  
  ## CI_95_upper
  
  for(ts in ts_spp_names){
    
    ts_summary_SSDM_bin$CI_95_upper[ts, model] <- quantile(ts_all_spp[ , ts, ],
                                                           probs = 0.975,
                                                           na.rm = TRUE)[[1]]
    
  }
  
  for(ts in ts_site_names){
    
    ts_summary_SSDM_bin$CI_95_upper[ts, model] <- quantile(ts_all_site[ , ts, ],
                                                           probs = 0.975,
                                                           na.rm = TRUE)[[1]]
    
  }
  
}

## Save to file

filename <- sprintf("outputs/test_statistics/ts_summary_SSDM_%s_bin.rds",
                    dataset_id)

saveRDS(ts_summary_SSDM_bin,
        filename)

#----

##################################
### Summarise Species Richness ###
##################################

## Model names

model_names <- c("MPR", "HPR", "LPR",
                 "DPR", "HLR_NS", "HLR_S",
                 "SSDM", "SESAM")

## Marginal ----

### Build empty array

sr_summary_marg <- matrix(NA,
                          nrow = 4,
                          ncol = length(model_names),
                          dimnames = list(c("mean",
                                            "median",
                                            "CI_95_lower",
                                            "CI_95_upper"),
                                          model_names))

for(model in model_names){
  
  ## Load test statistics for all folds
  
  sr_1 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold1_marginal_SR.rds",
                          model,
                          dataset_id))
  
  sr_2 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold2_marginal_SR.rds",
                          model,
                          dataset_id))
  
  sr_3 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold3_marginal_SR.rds",
                          model,
                          dataset_id))
  
  sr_4 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold4_marginal_SR.rds",
                          model,
                          dataset_id))
  
  sr_5 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold5_marginal_SR.rds",
                          model,
                          dataset_id))
  
  ## Merge folds
  
  sr_all <- abind(sr_1,
                  sr_2,
                  sr_3,
                  sr_4,
                  sr_5,
                  along = 1)
  
  ## Fill array of difference values
  
  sr_summary_marg["mean", model] <- mean(sr_all[ , "Difference", ],
                                         na.rm = TRUE)
  
  sr_summary_marg["median", model] <- median(sr_all[ , "Difference", ],
                                         na.rm = TRUE)
  
  sr_summary_marg["CI_95_lower", model] <- quantile(sr_all[ , "Difference", ],
                                                    probs = 0.025,
                                         na.rm = TRUE)[[1]]
  
  sr_summary_marg["CI_95_upper", model] <- quantile(sr_all[ , "Difference", ],
                                                    probs = 0.975,
                                                    na.rm = TRUE)[[1]]
  
  ## save to file
  
  filename <- sprintf("outputs/species_richness/SR_summary_JSDM_%s_marginal.rds",
                      dataset_id)
  
  saveRDS(sr_summary_marg,
          filename)
  
}

#----

## Conditional LOI ----

### Build empty array

sr_summary_condLOI_low <- matrix(NA,
                          nrow = 4,
                          ncol = length(model_names),
                          dimnames = list(c("mean",
                                            "median",
                                            "CI_95_lower",
                                            "CI_95_upper"),
                                          model_names))

sr_summary_condLOI_med <- matrix(NA,
                                 nrow = 4,
                                 ncol = length(model_names),
                                 dimnames = list(c("mean",
                                                   "median",
                                                   "CI_95_lower",
                                                   "CI_95_upper"),
                                                 model_names))

sr_summary_condLOI_high <- matrix(NA,
                                 nrow = 4,
                                 ncol = length(model_names),
                                 dimnames = list(c("mean",
                                                   "median",
                                                   "CI_95_lower",
                                                   "CI_95_upper"),
                                                 model_names))

for(model in model_names){
  
  ## Load test statistics for all folds
  
  sr_1 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold1_condLOI_SR.rds",
                          model,
                          dataset_id))
  
  sr_2 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold2_condLOI_SR.rds",
                          model,
                          dataset_id))
  
  sr_3 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold3_condLOI_SR.rds",
                          model,
                          dataset_id))
  
  sr_4 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold4_condLOI_SR.rds",
                          model,
                          dataset_id))
  
  sr_5 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold5_condLOI_SR.rds",
                          model,
                          dataset_id))
  
  ## Merge folds
  
  sr_all_low <- abind(sr_1[ , , , 1],
                      sr_2[ , , , 1],
                      sr_3[ , , , 1],
                      sr_4[ , , , 1],
                      sr_5[ , , , 1],
                      along = 1)
  
  sr_all_med <- abind(sr_1[ , , , 2],
                      sr_2[ , , , 2],
                      sr_3[ , , , 2],
                      sr_4[ , , , 2],
                      sr_5[ , , , 2],
                      along = 1)
  
  sr_all_high <- abind(sr_1[ , , , 3],
                      sr_2[ , , , 3],
                      sr_3[ , , , 3],
                      sr_4[ , , , 3],
                      sr_5[ , , , 3],
                      along = 1)
  
  ## Fill array of difference values
  
  sr_summary_condLOI_low["mean", model] <- mean(sr_all_low[ , "Difference", ],
                                         na.rm = TRUE)
  
  sr_summary_condLOI_med["mean", model] <- mean(sr_all_med[ , "Difference", ],
                                                na.rm = TRUE)
  
  sr_summary_condLOI_high["mean", model] <- mean(sr_all_high[ , "Difference", ],
                                                na.rm = TRUE)
  
  sr_summary_condLOI_low["median", model] <- median(sr_all_low[ , "Difference", ],
                                             na.rm = TRUE)
  
  sr_summary_condLOI_med["median", model] <- median(sr_all_med[ , "Difference", ],
                                                    na.rm = TRUE)
  
  sr_summary_condLOI_high["median", model] <- median(sr_all_high[ , "Difference", ],
                                                    na.rm = TRUE)
  
  sr_summary_condLOI_low["CI_95_lower", model] <- quantile(sr_all_low[ , "Difference", ],
                                                    probs = 0.025,
                                                    na.rm = TRUE)[[1]]
  
  sr_summary_condLOI_med["CI_95_lower", model] <- quantile(sr_all_med[ , "Difference", ],
                                                           probs = 0.025,
                                                           na.rm = TRUE)[[1]]
  
  sr_summary_condLOI_high["CI_95_lower", model] <- quantile(sr_all_high[ , "Difference", ],
                                                           probs = 0.025,
                                                           na.rm = TRUE)[[1]]
  
  sr_summary_condLOI_low["CI_95_upper", model] <- quantile(sr_all_low[ , "Difference", ],
                                                    probs = 0.975,
                                                    na.rm = TRUE)[[1]]
  
  sr_summary_condLOI_med["CI_95_upper", model] <- quantile(sr_all_med[ , "Difference", ],
                                                           probs = 0.975,
                                                           na.rm = TRUE)[[1]]
  
  sr_summary_condLOI_high["CI_95_upper", model] <- quantile(sr_all_high[ , "Difference", ],
                                                           probs = 0.975,
                                                           na.rm = TRUE)[[1]]
  
  ## save to file
  
  filename <- sprintf("outputs/species_richness/SR_summary_JSDM_%s_condLOI_low.rds",
                      dataset_id)
  
  saveRDS(sr_summary_condLOI_low,
          filename)
  
  filename <- sprintf("outputs/species_richness/SR_summary_JSDM_%s_condLOI_med.rds",
                      dataset_id)
  
  saveRDS(sr_summary_condLOI_med,
          filename)
  
  filename <- sprintf("outputs/species_richness/SR_summary_JSDM_%s_condLOI_high.rds",
                      dataset_id)
  
  saveRDS(sr_summary_condLOI_high,
          filename)
  
}

#----

## Joint ----

### Build empty array

sr_summary_joint <- matrix(NA,
                          nrow = 4,
                          ncol = length(model_names),
                          dimnames = list(c("mean",
                                            "median",
                                            "CI_95_lower",
                                            "CI_95_upper"),
                                          model_names))

for(model in model_names){
  
  ## Load test statistics for all folds
  
  sr_1 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold1_joint_SR.rds",
                          model,
                          dataset_id))
  
  sr_2 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold2_joint_SR.rds",
                          model,
                          dataset_id))
  
  sr_3 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold3_joint_SR.rds",
                          model,
                          dataset_id))
  
  sr_4 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold4_joint_SR.rds",
                          model,
                          dataset_id))
  
  sr_5 <- readRDS(sprintf("outputs/species_richness/%s_%s_fold5_joint_SR.rds",
                          model,
                          dataset_id))
  
  ## Merge folds
  
  sr_all <- abind(sr_1,
                  sr_2,
                  sr_3,
                  sr_4,
                  sr_5,
                  along = 1)
  
  ## Fill array of difference values
  
  sr_summary_joint["mean", model] <- mean(sr_all[ , "Difference", ],
                                         na.rm = TRUE)
  
  sr_summary_joint["median", model] <- median(sr_all[ , "Difference", ],
                                             na.rm = TRUE)
  
  sr_summary_joint["CI_95_lower", model] <- quantile(sr_all[ , "Difference", ],
                                                    probs = 0.025,
                                                    na.rm = TRUE)[[1]]
  
  sr_summary_joint["CI_95_upper", model] <- quantile(sr_all[ , "Difference", ],
                                                    probs = 0.975,
                                                    na.rm = TRUE)[[1]]
  
  ## save to file
  
  filename <- sprintf("outputs/species_richness/SR_summary_JSDM_%s_joint.rds",
                      dataset_id)
  
  saveRDS(sr_summary_joint,
          filename)
  
}

#----