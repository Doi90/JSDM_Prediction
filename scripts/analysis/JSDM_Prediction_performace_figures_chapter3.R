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

##################
##################
### MAKE PLOTS ###
##################
##################

dodge <- position_dodge(width = 0.5)

#########################
### Set colour scheme ###
#########################

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

###########################################
### Test Statistic-Specific Axis Limits ###
###########################################

graph_customisation <- read.csv("scripts/analysis/test_statistic_axis_limits.csv",
                                stringsAsFactors = FALSE)

#####################################################
### Set Prediction "Pairs" to put JSDMs and SSDMs ###
###       on same plot where applicable           ###
#####################################################

pred_pairs <- list(marginal_bin = c("marginal_bin", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Marginal"),
                   marginal_prob = c("marginal_prob", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Marginal"),
                   condLOI_low = c("condLOI_low", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (low)"),
                   condLOI_med = c("condLOI_med", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (med)"),
                   condLOI_high = c("condLOI_high", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Conditional (high)"),
                   condLOI_marg_low = c("condLOI_marg_low", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Conditional Marginal (low)"),
                   condLOI_marg_med = c("condLOI_marg_med", "SSDM_prob", "", "Probabilistic predictions", "JSDM prediction: Conditional Marginal (med)"),
                   condLOI_marg_high = c("condLOI_marg_high", "SSDM_bin", "SESAM", "Probabilistic predictions", "JSDM prediction: Conditional Marginal (high)"),
                   joint = c("joint", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Joint"))

################################################
### Generate plot for different combinations ###
################################################

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
    
    ## Extract values
    
    model_summary <- summary(mem_model)
    
    coefs <- model_summary$tTable
    
    ## Dataframe to store values
    
    n_row <- length(model_order)
    
    plot_df <- data.frame(model = factor(model_order,
                                         levels = model_order),
                          mean = numeric(n_row),
                          upper = numeric(n_row),
                          lower = numeric(n_row))
    
    for(i in seq_len(length(model_order))){
      
      model <- model_order[i]
      
      if(model != "SSDM" & model %nin% rownames(coefs)){
        
        plot_df[i, ] <- list(model,
                             NA,
                             NA,
                             NA)
      }
      
      if(model != "SSDM" & model %nin% rownames(coefs)){
        
        # plot_df[i, ] <- list(model,                                         # Model
        #                      coefs[paste0("model", model), "Value"],        # Mean
        #                      coefs[paste0("model", model), "Value"] +       # Upper
        #                        qnorm(0.975) * coefs[paste0("model", model), "Std.Error"],
        #                      coefs[paste0("model", model), "Value"] +       # Lower
        #                        qnorm(0.025) * coefs[paste0("model", model), "Std.Error"])
        
        pred_coef <- MEM_predict(mem = mem_model,
                                 model = model)
        
        mn <- pred_coef$mean
        se <- pred_coef$se
        
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
        
        plot_df[i, ] <- list(model,
                             mean,
                             quantile_fun(0.975, mn, se),
                             quantile_fun(0.025, mn, se))
        
      }
      
      if(model == "SSDM"){
        
        # plot_df[i, ] <- list(model,                                         # Model
        #                      coefs[paste0("model", model), "Value"],        # Mean
        #                      coefs[paste0("model", model), "Value"] +       # Upper
        #                        qnorm(0.975) * coefs[paste0("model", model), "Std.Error"],
        #                      coefs[paste0("model", model), "Value"] +       # Lower
        #                        qnorm(0.025) * coefs[paste0("model", model), "Std.Error"])
        
        pred_coef <- MEM_predict(mem = mem_model,
                                 model = model)
        
        mn <- pred_coef$mean
        se <- pred_coef$se
        
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
        
        plot_df[i, ] <- list(model,
                             mean,
                             quantile_fun(0.975, mn, se),
                             quantile_fun(0.025, mn, se))
        
      }
      
    }
    
    ## Make plot and save to file
    
    filename <- sprintf("outputs/test_statistics/plots/%1$s_%2$s_%3$s_MME.pdf",
                        names(pred_sets[prediction]),
                        ts,
                        chapter)
    pdf(filename)
    
    tmp_plot <- tryCatch(expr = ggplot(plot_df,
                                       aes(x = model,
                                           y = mean,
                                           colour = model)) +
                           geom_hline(yintercept = plot_df[plot_df$model == "SSDM", "mean"],
                                      linetype = 11) +
                           geom_hline(yintercept = plot_df[plot_df$model == "SSDM", "upper"],
                                      linetype = 11,
                                      colour = "lightgrey") +
                           geom_hline(yintercept = plot_df[plot_df$model == "SSDM", "lower"],
                                      linetype = 11,
                                      colour = "lightgrey") + 
                           geom_point(position = dodge,
                                      size = 2) + 
                           geom_errorbar(aes(ymax = upper,
                                             ymin = lower,),
                                         position = dodge,
                                         width = 0.1) +
                           ylim(c(graph_customisation[graph_customisation$test_statistic == ts, "y_min"],
                                  graph_customisation[graph_customisation$test_statistic == ts, "y_max"])) +
                           xlab("Model") +
                           ylab(eval(parse(text = graph_customisation[graph_customisation$test_statistic == ts, "graph_label"]))) +
                           scale_colour_manual(values = colour,
                                               breaks = levels(plot_df$model)) +
                           theme_bw() +
                           theme(legend.position = "none",
                                 panel.grid.minor = element_blank(),
                                 panel.grid.major = element_blank()) +
                           ggtitle(label = pred_pairs[[prediction]][4],
                                   subtitle = pred_pairs[[prediction]][5]),
                         error = function(err){
                           
                           return(NA)
                           
                         })
    
    if(chapter == "Ch3"){
      
      tmp_plot <- tmp_plot + geom_vline(xintercept = 2.5,
                                        lwd = 1.1)
      
    }
    
    tryCatch(expr = print(tmp_plot),
             error = function(err){
               
               message(sprintf("Plot failed for: %s - %s",
                               names(pred_sets)[prediction],
                               ts))
               
             }
    )
    
    dev.off()
  }
}
