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
                   #"HPR",
                   #"LPR",
                   #"DPR",
                   #"HLR_NS",
                   #"HLR_S",
                   "SSDM")#,
#"SESAM")

model_order <- c("SSDM",
                 #"SESAM",
                 "MPR")#,
#                 "HPR",
#                 "LPR",
#                 "DPR",
#                 "HLR_NS",
#                "HLR_S")

JSDM_models <- model_options[1]#:6]

SSDM_models <- model_options[2]#7:8]

## Datasets

dataset_options <- c("frog")#,
#"eucalypt")#,
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

#################################################
### Create empty dataframe to store plot info ###
#################################################

n_row <- length(dataset_options) * length(model_options) *
  length(pred_sets) * length(c(ts_site, ts_species))

plot_df <- data.frame(ts = character(n_row),
                      mean_diff = numeric(n_row),
                      upper = numeric(n_row),
                      lower = numeric(n_row),
                      dataset = character(n_row),
                      pred_type = character(n_row),
                      stringsAsFactors = FALSE)

###################################################
### Extract values to create plotting dataframe ###
###################################################

row_index <- 1

for(dataset in dataset_options){
  
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
      
      filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_%4$s_model.rds",
                          dataset,
                          names(pred_sets[prediction]),
                          ts,
                          chapter)
      
      if(!file.exists(filename)){
        
        message(sprintf("No model for: %s - %s - %s",
                        dataset,
                        names(pred_sets[prediction]),
                        ts))
        
        next()
        
      }
      
      mem_model <- readRDS(filename)
      
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
                 coef_matrix["modelMPR", "Value"])
      
      SSDM_mean <- coef_matrix["modelSSDM", "Value"]
      
      SSDM_se <- coef_matrix["modelSSDM", "Std.Error"]
      
      covar_matrix <- matrix(c(vc_matrix["modelSSDM", "modelSSDM"],
                               vc_matrix["modelSSDM", "modelMPR"],
                               vc_matrix["modelSSDM", "modelMPR"],
                               vc_matrix["modelMPR", "modelMPR"]),
                             nrow = 2, 
                             ncol = 2)
      
      draws <- mvrnorm(n = 1E6, 
                       means, 
                       covar_matrix)
      
      SSDM_orig_scale <- draws[ , 1]
      
      MPR_orig_scale <- draws[, 2]
      
      if(mem_model$transformed == FALSE | ts %in% ts_Inf_Inf){
        
        SSDM_coef <- SSDM_mean
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_Inf){
        
        SSDM_orig_scale <- exp(SSDM_orig_scale)
        
        MPR_orig_scale <- exp(MPR_orig_scale)
        
        SSDM_coef <- exp(SSDM_mean + (SSDM_se ^ 2 / 2)) 
        
      }
      
      if(mem_model$transformed == TRUE & ts %in% ts_0_1){
        
        SSDM_orig_scale <- plogis(SSDM_orig_scale)
        
        MPR_orig_scale <- plogis(MPR_orig_scale)
        
        SSDM_coef <- mean(plogis(rnorm(1000000, SSDM_mean, SSDM_se)))
      }
      
      diff_samples <- SSDM_orig_scale - MPR_orig_scale
      
      diff_samples <- (diff_samples / SSDM_coef) * 100 - 100
      
      ## Divide by backtransformed regression coefficient here
      
      plot_df[row_index, ] <- list(ts,
                                   mean(diff_samples),
                                   quantile(diff_samples, 0.975),
                                   quantile(diff_samples, 0.025),
                                   dataset,
                                   names(pred_sets)[prediction])
      
      row_index <- row_index + 1
      
    }
  }
}

plot_df <- plot_df[plot_df$ts != "", ]

################
### Plotting ###
################

dodge <- position_dodge(width = 0.5)

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
                               "SSE"))

## Plot

for(dataset in unique(plot_df$dataset)){
  
  for(pred_type in unique(plot_df$pred_type)){
    
    for(subset in names(subsets)){
    
    tmp_df <- plot_df[plot_df$dataset == dataset &
                        plot_df$pred_type == pred_type &
                        plot_df$ts %in% subsets[[subset]], ]
    
    if(nrow(tmp_df) == 0){
      
      next()
      
    }
    
    tmp_df$ts <- factor(tmp_df$ts,
                        levels = sort(subsets[[subset]]))
    
    ## Make plot and save to file
    
    filename <- sprintf("outputs/test_statistics/plots/relative_performance_%1$s_%2$s_%3$s_%4$s_MME.pdf",
                        dataset,
                        pred_type,
                        subset,
                        chapter)
    
    pdf(filename)
    
    ### Start plot
    
    tmp_plot <- ggplot(tmp_df,
                       aes(x = ts,
                           y = mean_diff)) +
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
}
