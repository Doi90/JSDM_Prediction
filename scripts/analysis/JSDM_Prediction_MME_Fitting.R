


#####################
### Load Packages ###
#####################

library(nlme)
library(ggplot2)
library(RColorBrewer)
library(PassButter)

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
                             "SSDM_prob")


## Define prediction combinations for subsetting

pred_sets <- list(marginal_bin = c("marginal_bin", "SSDM_bin", "SESAM"),
                  marginal_prob = c("marginal_prob", "SSDM_prob"),
                  condLOI_low = c("condLOI_low", "SSDM_bin", "SESAM"),
                  condLOI_med = c("condLOI_med", "SSDM_bin", "SESAM"),
                  condLOI_high = c("condLOI_high", "SSDM_bin", "SESAM"),
                  condLOI_marg_low = c("condLOI_marg_low", "SSDM_prob"),
                  condLOI_marg_med = c("condLOI_marg_med", "SSDM_prob"),
                  condLOI_marg_high = c("condLOI_marg_high", "SSDM_prob"),
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
                "TP",
                "FP",
                "TN",
                "FN",
                "TPR",
                "FPR",
                "TNR",
                "FNR",
                "PLR",
                "NLR",
                "DOR",
                "Prevalence",
                "Accuracy",
                "PPV",
                "FOR",
                "FDR",
                "NPV",
                "F_1",
                "Youden_J",
                "Kappa")

### By Site

ts_site <- c("Binomial",
             "Bray",
             "Canberra",
             "Euclidean",
             "Gower",
             "Gower_alt",
             "Horn",
             "Jaccard",
             "Kulczynski",
             "Mahalanobis",
             "Manhattan",
             "Mountford",
             "Raup")

## Test statistic / prediction type compatibility

binary_ts <- c("TP",
               "FP",
               "TN",
               "FN",
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
               "Binomial",
               "Bray",
               "Canberra",
               "Euclidean",
               "Gower",
               "Gower_alt",
               "Horn",
               "Jaccard",
               "Kulczynski",
               "Mahalanobis",
               "Manhattan",
               "Mountford",
               "Raup")

prob_ts <- c("AUC",
             "bias",
             "MSE",
             "R2",
             "RMSE",
             "SSE",
             "Pearson",
             "Spearman",
             "Kendall",
             "TP",
             "FP",
             "TN",
             "FN",
             "TPR",
             "FPR",
             "TNR",
             "FNR",
             "PLR",
             "NLR",
             "DOR",
             "Prevalence",
             "Accuracy",
             "PPV",
             "FOR",
             "FDR",
             "NPV",
             "F_1",
             "Youden_J")

## Chapter

if(length(model_options) == 2){
  
  chapter <- "Ch2"
  
} else {
  
  chapter <- "Ch3"
  
}

#######################
#######################
### TEST STATISTICS ###
#######################
#######################

#################
### Load Data ###
#################

ts_df_species <- readRDS(sprintf("outputs/test_statistics/test_statistics_species_summary_MME_%s.rds",
                                 chapter))

ts_df_site <- readRDS(sprintf("outputs/test_statistics/test_statistics_site_summary_MME_%s.rds",
                              chapter))

############################
### Mixed Effects Models ###
############################

## Empty dataframe to store Kolmogorov-Smirnov test outputs

n_row <- length(unique(ts_df_species$dataset)) *
  (length(unique(ts_df_species$test_statistic)) + 
     length(unique(ts_df_site$test_statistic))) *
  length(pred_sets)

ks_df <- data.frame(dataset = factor(character(n_row),
                                     levels = unique(ts_df_species$dataset)),
                    pred_type = factor(character(n_row),
                                       levels = names(pred_sets)),
                    test_statistic = factor(character(n_row),
                                            levels = c(unique(as.character(ts_df_species$test_statistic)),
                                                       unique(as.character(ts_df_site$test_statistic)))),
                    p_value = numeric(n_row),
                    ks_D = numeric(n_row),
                    model_fit_success = numeric(n_row))

## Loop over different model iterations

row_index <- 1

## Species-based statistics

for(dataset in unique(ts_df_species$dataset)){
  
  for(prediction in seq_len(length(pred_sets))){
    
    for(ts in ts_species){
      
      ## Check ts/prediction compatibility
      
      ### Binary test statistics for binary predictions
      
      if(ts %in% binary_ts & names(pred_sets)[prediction] %nin% binary_predictions){
        
        next()
        
      }
      
      ### Probabilistic test statistics for probabilistic predictions
      
      if(ts %in% prob_ts & names(pred_sets)[prediction] %nin% probability_predictions){
        
        next()
        
      }
      
      ## Subset dataset
      
      tmp_df <- ts_df_species[ts_df_species$dataset == dataset &
                                ts_df_species$prediction_type %in% pred_sets[[prediction]] &
                                ts_df_species$test_statistic == ts, ]
      
      tmp_df <- na.omit(tmp_df)
      
      ## Fit mixed effects model
      
      mme_model <- tryCatch(expr = nlme::lme(mean ~ -1 + model + fold, 
                                             random = ~ 1|species,
                                             weights = varIdent(form = ~1|model),
                                             data = tmp_df,
                                             control = lmeControl(msMaxIter = 1000,
                                                                  opt = "optim")),
                            error = function(err){
                              
                              message(sprintf("Model fit failed for: %s - %s - %s",
                                              dataset,
                                              names(pred_sets[prediction]),
                                              ts))
                              
                              return(NA)
                              
                            }
      )
      
      if(class(mme_model) != "lme"){
        
        ks_df[row_index, ] <- list(dataset,
                                   names(pred_sets[prediction]),
                                   ts,
                                   NA,
                                   NA,
                                   0)
        
        row_index <- row_index + 1
        
        next()
        
      }
      
      ## Extract residuals
      
      residuals <- resid(mme_model, type = "pearson")
      
      ## Plot residuals and save to file
      
      ### Set filename
      
      if(length(model_options) == 2){
        
        chapter <- "Ch2"
        
      } else {
        
        chapter <- "Ch3"
        
      }
      
      filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s_%4$s.pdf",
                          dataset,
                          names(pred_sets[prediction]),
                          ts,
                          chapter)
      
      ### Plot to pdf
      
      pdf(filename)
      
      ### Residuals boxplot
      
      plot(residuals ~ tmp_df$model,
           ylab = "Pearson residuals",
           xlab = "Model")
      
      ### Residuals qq plot
      
      qqplot(residuals, 
             rnorm(n = 1000, 
                   mean(residuals), 
                   sd(residuals)),
             ylab = "N(mean(residuals), sd(residuals))",
             xlab = "Model Pearson residuals",
             main = "QQ-plot")
      abline(0,1)
      
      dev.off()
      
      ### Check normality with Kolmogorov-Smirnov test
      
      ks <- ks.test(residuals, 
                    pnorm, 
                    mean(residuals), 
                    sd(residuals))
      
      ks_df[row_index, ] <- list(dataset,
                                 names(pred_sets[prediction]),
                                 ts,
                                 ks$p.value,
                                 ks$statistic,
                                 1)
      
      ### Save model to file
      
      filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_%4$s_model.rds",
                          dataset,
                          names(pred_sets[prediction]),
                          ts,
                          chapter)
      
      saveRDS(mme_model,
              filename)
      
      ### Statusbar
      
      statusbar(run = row_index,
                max.run = n_row,
                info = sprintf("%s: %s - %s - %-15s",
                               row_index,
                               dataset,
                               names(pred_sets[prediction]),
                               ts))
      
      row_index <- row_index + 1
      
    }
  }
}

## Site-based statistics

for(dataset in unique(ts_df_site$dataset)){
  
  for(prediction in seq_len(length(pred_sets))){
    
    for(ts in ts_site){
      
      ## Check ts/prediction compatibility
      
      ### Binary test statistics for binary predictions
      
      if(ts %in% binary_ts & names(pred_sets)[prediction] %nin% binary_predictions){
        
        next()
        
      }
      
      ### Probabilistic test statistics for probabilistic predictions
      
      if(ts %in% prob_ts & names(pred_sets)[prediction] %nin% probability_predictions){
        
        next()
        
      }
      
      ## Subset dataset
      
      tmp_df <- ts_df_site[ts_df_site$dataset == dataset &
                             ts_df_site$prediction_type %in% pred_sets[[prediction]] &
                             ts_df_site$test_statistic == ts, ]
      
      tmp_df <- na.omit(tmp_df)
      
      ## Fit mixed effects model
      
      mme_model <- tryCatch(expr = nlme::lme(mean ~ -1 + model + fold, 
                                             random = ~ 1|site,
                                             weights = varIdent(form = ~1|model),
                                             data = tmp_df,
                                             control = lmeControl(msMaxIter = 1000,
                                                                  opt = "optim")),
                            error = function(err){
                              
                              message(sprintf("Model fit failed for: %s - %s - %s",
                                              dataset,
                                              names(pred_sets[prediction]),
                                              ts))
                              
                              return(NA)
                              
                            }
      )
      
      if(class(mme_model) != "lme"){
        
        ks_df[row_index, ] <- list(dataset,
                                   names(pred_sets[prediction]),
                                   ts,
                                   NA,
                                   NA,
                                   0)
        
        row_index <- row_index + 1
        
        next()
        
      }
      
      ## Extract residuals
      
      residuals <- resid(mme_model, type = "pearson")
      
      ## Plot residuals and save to file
      
      ### Set filename
      
      if(length(model_options) == 2){
        
        chapter <- "Ch2"
        
      } else {
        
        chapter <- "Ch3"
        
      }
      
      filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s_%4$s.pdf",
                          dataset,
                          names(pred_sets[prediction]),
                          ts,
                          chapter)
      
      ### Plot to pdf
      
      pdf(filename)
      
      ### Residuals boxplot
      
      plot(residuals ~ tmp_df$model,
           ylab = "Pearson residuals",
           xlab = "Model")
      
      ### Residuals qq plot
      
      qqplot(residuals, 
             rnorm(n = 1000, 
                   mean(residuals), 
                   sd(residuals)),
             ylab = "N(mean(residuals), sd(residuals))",
             xlab = "Model Pearson residuals",
             main = "QQ-plot")
      abline(0,1)
      
      dev.off()
      
      ### Check normality with Kolmogorov-Smirnov test
      
      ks <- ks.test(residuals, 
                    pnorm, 
                    mean(residuals), 
                    sd(residuals))
      
      ks_df[row_index, ] <- list(dataset,
                                 names(pred_sets[prediction]),
                                 ts,
                                 ks$p.value,
                                 ks$statistic,
                                 1)
      
      ### Save model to file
      
      filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_%4$s_model.rds",
                          dataset,
                          names(pred_sets[prediction]),
                          ts,
                          chapter)
      
      saveRDS(mme_model,
              filename)
      
      ### Statusbar
      
      statusbar(run = row_index,
                max.run = n_row,
                info = sprintf("%s: %s - %s - %-15s",
                               row_index,
                               dataset,
                               names(pred_sets[prediction]),
                               ts))
      
      row_index <- row_index + 1
      
    }
  }
}

########################
########################
### SPECIES RICHNESS ###
########################
########################

#################
### Load Data ###
#################

sr_df <- readRDS(sprintf("outputs/species_richness/species_richness_summary_%s.rds",
                         chapter))

############################
### Mixed Effects Models ###
############################

## Empty dataframe to store Kolmogorov-Smirnov test outputs

n_row <- length(unique(sr_df$dataset)) *
  length(unique(sr_df$test_statistic)) *
  length(pred_sets)

ks_sr <- data.frame(dataset = factor(character(n_row),
                                     levels = unique(sr_df$dataset)),
                    pred_type = factor(character(n_row),
                                       levels = names(pred_sets)),
                    test_statistic = factor(character(n_row),
                                            levels = unique(as.character(sr_df$test_statistic))),
                    p_value = numeric(n_row),
                    ks_D = numeric(n_row),
                    model_fit_success = numeric(n_row))

## Loop over different different model iterations

row_index <- 1

for(dataset in unique(sr_df$dataset)){
  
  for(prediction in seq_len(length(pred_sets))){
    
    ## Subset dataset
    
    tmp_df <- sr_df[sr_df$dataset == dataset &
                      sr_df$prediction_type %in% pred_sets[[prediction]], ]
    
    tmp_df <- na.omit(tmp_df)
    
    ## Fit mixed effects model
    
    mme_model <- tryCatch(expr = nlme::lme(mean ~ -1 + model + fold, 
                                           random = ~ 1|site,
                                           weights = varIdent(form = ~1|model),
                                           data = tmp_df,
                                           control = lmeControl(msMaxIter = 1000,
                                                                opt = "optim")),
                          error = function(err){
                            
                            message(sprintf("Model fit failed for: %s - %s - %s",
                                            dataset,
                                            names(pred_sets[prediction]),
                                            "species_richness_difference"))
                            
                            return(NA)
                            
                          }
    )
    
    if(class(mme_model) != "lme"){
      
      ks_sr[row_index, ] <- list(dataset,
                                 names(pred_sets[prediction]),
                                 "species_richness_difference",
                                 NA,
                                 NA,
                                 0)
      
      row_index <- row_index + 1
      
      next()
      
    }
    
    ## Extract residuals
    
    residuals <- resid(mme_model, type = "pearson")
    
    ## Plot residuals and save to file
    
    ### Set filename
    
    if(length(model_options) == 2){
      
      chapter <- "Ch2"
      
    } else {
      
      chapter <- "Ch3"
      
    }
    
    filename <- sprintf("outputs/test_statistics/plots/test_statistics_MME_%1$s_%2$s_%3$s_%4$s.pdf",
                        dataset,
                        names(pred_sets[prediction]),
                        "SR",
                        chapter)
    
    ### Plot to pdf
    
    pdf(filename)
    
    ### Residuals boxplot
    
    plot(residuals ~ tmp_df$model,
         ylab = "Pearson residuals",
         xlab = "Model")
    
    ### Residuals qq plot
    
    qqplot(residuals, 
           rnorm(n = 1000, 
                 mean(residuals), 
                 sd(residuals)),
           ylab = "N(mean(residuals), sd(residuals))",
           xlab = "Model Pearson residuals",
           main = "QQ-plot")
    abline(0,1)
    
    dev.off()
    
    ### Check normality with Kolmogorov-Smirnov test
    
    ks <- ks.test(residuals, 
                  pnorm, 
                  mean(residuals), 
                  sd(residuals))
    
    ks_sr[row_index, ] <- list(dataset,
                               names(pred_sets[prediction]),
                               "species_richness_difference",
                               ks$p.value,
                               ks$statistic,
                               1)
    
    ### Save model to file
    
    filename <- sprintf("outputs/test_statistics/models/%1$s_%2$s_%3$s_%4$s_model.rds",
                        dataset,
                        names(pred_sets[prediction]),
                        "species_richness_difference",
                        chapter)
    
    saveRDS(mme_model,
            filename)
    
    ### Statusbar
    
    statusbar(run = row_index,
              max.run = n_row,
              info = sprintf("%s: %s - %s - SR",
                             row_index,
                             dataset,
                             names(pred_sets[prediction])))
    
    row_index <- row_index + 1
    
  }
}

ks_df <- rbind(ks_df,
               ks_sr)











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
                   joint = c("joint", "SSDM_bin", "SESAM", "Binary predictions", "JSDM prediction: Joint"))

################################################
### Generate plot for different combinations ###
################################################

for(dataset in unique(ts_df_species$dataset)){
  
  for(prediction in seq_len(length(pred_sets))){
    
    for(ts in c(ts_species, ts_site, "species_richness_difference")){
      
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
      
      mme_model <- readRDS(filename)
      
      ## Extract values
      
      model_summary <- summary(mme_model)
      
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
        
        if(paste0("model", model) %nin% rownames(coefs)){
          
          plot_df[i, ] <- list(model,
                               NA,
                               NA,
                               NA)
        }
        
        if(paste0("model", model) %in% rownames(coefs)){
          
          plot_df[i, ] <- list(model,                                         # Model
                               coefs[paste0("model", model), "Value"],        # Mean
                               coefs[paste0("model", model), "Value"] +       # Upper
                                 coefs[paste0("model", model), "Std.Error"],
                               coefs[paste0("model", model), "Value"] - 
                                 coefs[paste0("model", model), "Std.Error"])  # Lower
          
        }
      }
      
      ## Make plot and save to file
      
      filename <- sprintf("outputs/test_statistics/plots/%1$s_%2$s_%3$s_%4$s_MME.pdf",
                          dataset,
                          names(pred_sets[prediction]),
                          ts,
                          chapter)
      pdf(filename)
      
      tmp_plot <- ggplot(plot_df,
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
                subtitle = pred_pairs[[prediction]][5])
      
      if(chapter == "Ch3"){
        
        tmp_plot <- tmp_plot + geom_vline(xintercept = 2.5,
                                          lwd = 1.1)
        
      }
      
      print(tmp_plot)
      
      dev.off()
    }
  }
}
