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
# library(circlize)

##############################
### Load Prediction Script ###
##############################

source("scripts/prediction/MEM_predict.R")

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

model_order_binary <- c("SESAM",
                        "MPR",
                        "HPR",
                        "LPR",
                        "DPR",
                        "HLR_NS",
                        "HLR_S")

model_order_prob <- c("MPR",
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
                        "joint")

probability_predictions <- c("marginal_prob",
                             "condLOI_marg_low",
                             "condLOI_marg_med",
                             "condLOI_marg_high")

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
  # "PLR",
  # "NLR",
  # "DOR",
  # "Accuracy",
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
             # "PLR",
             # "NLR",
             # "DOR",
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
                               "SSE"))#,
                               # "PLR",
                               # "NLR",
                               # "DOR"))

## Chapter

if(length(model_options) == 2){
  
  chapter <- "Ch2"
  
} else {
  
  chapter <- "Ch3"
  
}

########################################################
### Load MEMs and create big dataframe of all values ###
########################################################

## Build empty data.frame

n_row <- length(model_order_binary) * length(binary_predictions) * length(binary_ts) +
  length(model_order_prob) * length(probability_predictions) * length(prob_ts)

plot_df <- data.frame(model = character(n_row),
                      metric = character(n_row),
                      pred_type = character(n_row),
                      pred_class = character(n_row),
                      value = numeric(n_row),
                      rel_value = numeric(n_row))

## Extract values

row_index <- 1

for(prediction in binary_predictions){
  
  for(ts in binary_ts){
    
    ## Load Model
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        prediction,
                        ts,
                        chapter)
    
    if(!file.exists(filename)){
      
      message(sprintf("No model for: %s - %s",
                      prediction,
                      ts))
      
      next()
      
    }
    
    mem_model <- readRDS(filename)
    
    ## SSDM calculation
    
    SSDM_pred <- MEM_predict(mem = mem_model,
                             model = "SSDM")
    
    SSDM_mn <- SSDM_pred$mean
    SSDM_se <- SSDM_pred$se
    
    if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
      
      quantile_fun <- qnorm
      
      SSDM_mean <- SSDM_mn
      
    }
    
    if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
      
      quantile_fun <- qlnorm
      
      SSDM_mean <- exp(SSDM_mn + (SSDM_se ^ 2 / 2))
      
    }
    
    if(mem_model$transformed == TRUE & ts %in% ts_0_1){
      
      quantile_fun <- qlogitnorm
      
      SSDM_mean <- mean(plogis(rnorm(1000000, SSDM_mn, SSDM_se)))
      
    }
    
    plot_df[row_index, ] <- list("SSDM",
                                 ts,
                                 prediction,
                                 "binary",
                                 SSDM_mean,
                                 0)
    
    row_index <- row_index + 1
    
    ## Other model calculations
    
    for(model in model_order_binary){
      
      pred <- MEM_predict(mem = mem_model,
                          model = model)
      
      mn <- pred$mean
      se <- pred$se
      
      if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
        
        quantile_fun <- qnorm
        
        mean <- mn
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
        
        quantile_fun <- qlnorm
        
        mean <- exp(mn + (se ^ 2 / 2))
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_1){
        
        quantile_fun <- qlogitnorm
        
        mean <- mean(plogis(rnorm(1000000, mn, se)))
        
      }
      
      rel_mean <- ((mean - SSDM_mean) / abs(SSDM_mean)) * 100
      
      plot_df[row_index, ] <- list(model,
                                   ts,
                                   prediction,
                                   "binary",
                                   mean,
                                   rel_mean)
      
      row_index <- row_index + 1
      
    }
  }
}

for(prediction in probability_predictions){
  
  for(ts in prob_ts){
    
    ## Load Model
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_model.rds",
                        prediction,
                        ts,
                        chapter)
    
    if(!file.exists(filename)){
      
      message(sprintf("No model for: %s - %s",
                      prediction,
                      ts))
      
      next()
      
    }
    
    mem_model <- readRDS(filename)
    
    ## SSDM calculation
    
    SSDM_pred <- MEM_predict(mem = mem_model,
                             model = "SSDM")
    
    SSDM_mn <- SSDM_pred$mean
    SSDM_se <- SSDM_pred$se
    
    if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
      
      quantile_fun <- qnorm
      
      SSDM_mean <- SSDM_mn
      
    }
    
    if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
      
      quantile_fun <- qlnorm
      
      SSDM_mean <- exp(SSDM_mn + (SSDM_se ^ 2 / 2))
      
    }
    
    if(mem_model$transformed == TRUE & ts %in% ts_0_1){
      
      quantile_fun <- qlogitnorm
      
      SSDM_mean <- mean(plogis(rnorm(1000000, SSDM_mn, SSDM_se)))
      
    }
    
    plot_df[row_index, ] <- list("SSDM",
                                 ts,
                                 prediction,
                                 "probabilistic",
                                 SSDM_mean,
                                 0)
    
    row_index <- row_index + 1
    
    ## Other model calculations
    
    for(model in model_order_binary){
      
      pred <- MEM_predict(mem = mem_model,
                          model = model)
      
      mn <- pred$mean
      se <- pred$se
      
      if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
        
        quantile_fun <- qnorm
        
        mean <- mn
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
        
        quantile_fun <- qlnorm
        
        mean <- exp(mn + (se ^ 2 / 2))
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_1){
        
        quantile_fun <- qlogitnorm
        
        mean <- mean(plogis(rnorm(1000000, mn, se)))
        
      }
      
      rel_mean <- ((mean - SSDM_mean) / abs(SSDM_mean)) * 100
      
      plot_df[row_index, ] <- list(model,
                                   ts,
                                   prediction,
                                   "probabilistic",
                                   mean,
                                   rel_mean)
      
      row_index <- row_index + 1
      
    }
  }
}







#########################
### Old circlize code ###
#########################

# for(pred in unique(plot_df$pred_type)){
#   
#   ## Subset Data
#   
#   if(pred %in% binary_predictions){
#     
#     tmp_circos_df <- plot_df[plot_df$pred_type == pred &
#                                plot_df$ts %in% binary_ts, ]
#     
#   }
#   
#   if(pred %in% probability_predictions){
#     
#     tmp_circos_df <- plot_df[plot_df$pred_type == pred &
#                                plot_df$ts %in% prob_ts, ]
#     
#   }
#   
#   ## Reformat Data
#   
#   ### Blank matrix
#   
#   if(pred %in% binary_predictions){
#     
#     plot_mat <- matrix(NA,
#                        nrow = length(unique(tmp_circos_df$model)),
#                        ncol = sum(unlist(lapply(lapply(subsets, 
#                                                        function(x) x %in% binary_ts), 
#                                                 sum))),
#                        dimnames = list(model_order_binary,
#                                        c(binary_ts[binary_ts %in% subsets[["ts_0_1_H"]]],
#                                          binary_ts[binary_ts %in% subsets[["ts_0_1_L_a"]]],
#                                          binary_ts[binary_ts %in% subsets[["ts_0_1_L_b"]]],
#                                          binary_ts[binary_ts %in% subsets[["ts_1_1"]]],
#                                          binary_ts[binary_ts %in% subsets[["ts_Inf_Inf_H"]]],
#                                          binary_ts[binary_ts %in% subsets[["ts_Inf_Inf_0"]]],
#                                          binary_ts[binary_ts %in% subsets[["ts_0_Inf_L"]]])))
#     
#   }
#   
#   if(pred %in% probability_predictions){
#     
#     plot_mat <- matrix(NA,
#                        nrow = length(unique(tmp_circos_df$model)),
#                        ncol = sum(unlist(lapply(lapply(subsets, 
#                                                        function(x) x %in% prob_ts), 
#                                                 sum))),
#                        dimnames = list(model_order_prob,
#                                        c(prob_ts[prob_ts %in% subsets[["ts_0_1_H"]]],
#                                          prob_ts[prob_ts %in% subsets[["ts_0_1_L_a"]]],
#                                          prob_ts[prob_ts %in% subsets[["ts_0_1_L_b"]]],
#                                          prob_ts[prob_ts %in% subsets[["ts_1_1"]]],
#                                          prob_ts[prob_ts %in% subsets[["ts_Inf_Inf_H"]]],
#                                          prob_ts[prob_ts %in% subsets[["ts_Inf_Inf_0"]]],
#                                          prob_ts[prob_ts %in% subsets[["ts_0_Inf_L"]]])))
#     
#   }
#   
#   ### Fill with values
#   
#   for(model_id in rownames(plot_mat)){
#     
#     for(ts_id in colnames(plot_mat)){
#       
#       plot_mat[model_id, ts_id] <- tmp_circos_df[tmp_circos_df$model == model_id &
#                                                    tmp_circos_df$ts == ts_id, "mean_diff"]
#       
#     }
#   }
#   
#   ## Make Circos plot
#   
#   col_fun <- colorRamp2(c(-100, 0, 100), c("Red", "White", "Blue"))
#   
#   factors <- rep(c("ts_0_1_H", 
#                    "ts_0_1_L_a",
#                    "ts_0_1_L_b",
#                    "ts_1_1",
#                    "ts_Inf_Inf_H", 
#                    "ts_Inf_Inf_0",
#                    "ts_0_Inf_L"),
#                times = c(sum(subsets[["ts_0_1_H"]] %in% prob_ts), 
#                          sum(subsets[["ts_0_1_L_a"]] %in% prob_ts),
#                          sum(subsets[["ts_0_1_L_b"]] %in% prob_ts), 
#                          sum(subsets[["ts_1_1"]] %in% prob_ts), 
#                          sum(subsets[["ts_Inf_Inf_H"]] %in% prob_ts), 
#                          sum(subsets[["ts_Inf_Inf_0"]] %in% prob_ts), 
#                          sum(subsets[["ts_0_Inf_L"]] %in% prob_ts)))
# 
#   
#   mat_list <- list(ts_0_1_H = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_0_1_H"]])],
#                    ts_0_1_L_a = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_0_1_L_a"]])],
#                    ts_0_1_L_b = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_0_1_L_b"]])],
#                    ts_1_1 = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_1_1"]])],
#                    ts_Inf_Inf_H = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_Inf_Inf_H"]])],
#                    ts_Inf_Inf_0 = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_Inf_Inf_0"]])],
#                    ts_0_Inf_L = plot_mat[ , which(colnames(plot_mat) %in% subsets[["ts_0_Inf_L"]])])
#   
#   dend_list <- list(ts_0_1_H = as.dendrogram(hclust(dist(t(mat_list[["ts_0_1_H"]])))),
#                     ts_0_1_L_a = as.dendrogram(hclust(dist(t(mat_list[["ts_0_1_L_a"]])))),
#                     ts_0_1_L_b = as.dendrogram(hclust(dist(t(mat_list[["ts_0_1_L_b"]])))),
#                     ts_1_1 = as.dendrogram(hclust(dist(t(mat_list[["ts_1_1"]])))),
#                     ts_Inf_Inf_H = as.dendrogram(hclust(dist(t(mat_list[["ts_Inf_Inf_H"]])))),
#                     ts_Inf_Inf_0 = as.dendrogram(hclust(dist(t(mat_list[["ts_Inf_Inf_0"]])))),
#                     ts_0_Inf_L = as.dendrogram(hclust(dist(t(mat_list[["ts_0_Inf_L"]])))))
#   
#   filename <- sprintf("outputs/test_statistics/plots/circos_%s_chapter3.pdf",
#                       pred)
#   
#   pdf(filename,
#       width = 10,
#       height = 10)
#   
#   circos.par(cell.padding = c(0, 
#                               0, 
#                               0, 
#                               0), 
#              gap.degree = 10,
#              "track.height" = 0.8)
#   
#   circos.initialize(factors = factors,
#                     xlim = cbind(0, table(factors)))
#   
#   circos.track(ylim = c(0, 6), 
#                bg.border = NA, 
#                panel.fun = function(x, y){
#                  
#                  sector.index = CELL_META$sector.index
#                  m = mat_list[[sector.index]]
#                  dend = dend_list[[sector.index]]
#                  
#                  m2 = m[, order.dendrogram(dend)]
#                  col_mat = col_fun(m2)
#                  nr = nrow(m2)
#                  nc = ncol(m2)
#                  for(i in 1:nr) {
#                    circos.rect(1:nc - 1, rep(nr - i, nc), 
#                                1:nc, rep(nr - i + 1, nc), 
#                                border = col_mat[i, ], col = col_mat[i, ])
#                    
#                  }
#                  circos.axis(h = "top",
#                              major.tick = FALSE,
#                              labels = colnames(mat_list[[sector.index]]),
#                              labels.cex = 0.6,
#                              labels.facing = "bending.inside",
#                              labels.away.percentage = 0.3)
#                  circos.yaxis(labels = model_order_prob, 
#                               labels.cex = 0.6)
#                })
#   
#   circos.clear()
#   
#   dev.off()
# }
