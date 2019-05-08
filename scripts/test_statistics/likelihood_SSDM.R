######################################################
######################################################
###                                                ###
###      SSDM/SESAM LOG LIKELIHOOD WORKFLOW        ###
###                                                ###
###    This script calculates the log-likelihood   ###
### statistics used in this analysis for SSDMs and ###
### SESAM models.                                  ###
###                                                ###
######################################################
######################################################

#################
### Load Data ###
#################

## Likelihood is calculated using the data used to fit the model
## so we need to load in the TRAINING data

## Pres/Abs data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y <- eval(parse(text = command))                        # Evaluate command to read in data

y <- y[ , -1]                                           # Remove rownames

########################
### Load Predictions ###
########################

## SSDM - probabilities

filename <- sprintf("outputs/predictions/%s_prob_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

SSDM_predictions_prob <- readRDS(filename)

# ## SSDM - binary
# 
# filename <- sprintf("outputs/predictions/%s_bin_%s_fold%s.rds",
#                     model_id,
#                     dataset_id,
#                     fold_id)
# 
# SSDM_predictions_bin <- readRDS(filename)

# ## SESAM
# 
# filename <- sprintf("outputs/predictions/SESAM_%s_fold%s.rds",
#                     dataset_id,
#                     fold_id)
# 
# SESAM_predictions <- readRDS(filename)

################################
### Calculate Log-Likelihood ###
################################ 

# SSDM - probabilities

## Independent likelihood

### Vector to store likelihoods

SSDM_log_lik_vector <- vector(length = dim(SSDM_predictions_prob)[2])

### Calculate likelihoods

for(i in seq_len(dim(SSDM_predictions_prob)[2])){
  
  ### Get predicted species probabilities for site
  
  pred_probs <- SSDM_predictions_prob[ , i, 1]
  
  ### Get observed species state for site
  
  obs_spp <- y[ , i]
  
  ### Calculate log likelihood for site
  
  log_lik <- sum(dbinom(x = unlist(obs_spp),
                        size = 1,
                        prob = pred_probs,
                        log = TRUE))
  
  ### Add value to vector
  
  SSDM_log_lik_vector[i] <- log_lik
  
}

# ## Obtain single likelihood value for whole model
# 
# SSDM_log_likelihood_estimate <- sum(SSDM_log_lik_vector,
#                                     na.rm = TRUE)

## #Save to file

filename <- sprintf("outputs/likelihood/SSDM_prob_%s_fold%s_independent_likelihood.rds",
                    dataset_id,
                    fold_id)

saveRDS(SSDM_log_lik_vector,
        filename)

## Joint likelihood

### Vector to store likelihoods

SSDM_log_lik_vector <- vector(length = dim(SSDM_predictions_prob)[1])

### Calculate likelihoods

for(i in seq_len(dim(SSDM_predictions_prob)[1])){
  
  ### Get predicted species probabilities for site
  
  pred_probs <- SSDM_predictions_prob[i, , 1]
  
  ### Get observed species state for site
  
  obs_spp <- y[i, ]
  
  ### Mu and Sigma
  
  mu <- qnorm(pred_probs)
  
  sigma <- diag(length(pred_probs))
  
  ### Define probability distribution thresholds
  
  #### lower / upper to limit integral of density for lielihood
  
  lower <- rep(-Inf, length(pred_probs))  # default vector of -Inf lower limits
  upper <- rep(+Inf, length(pred_probs))  # default vector of +Inf upper limits
  
  for(k in seq_len(length(pred_probs))){  # set actual lower/upper limits based on known occurrence states
    
    if(obs_spp[k] == 0){     # if species is absent
      
      upper[k] <- 0            # species absent when z<0
      
    } 
    
    if(obs_spp[k] == 1){     # if species is present
      
      lower[k] <- 0            # species present when z>0
      
    } 
    
  }
  
  #### Prediction for species assemblage at site i using values from slice a
  
  likelihood_tmp <- lik_probit(obs = obs_spp,
                               mu = mu,
                               R = sigma,
                               log.p = FALSE,
                               niter = 1000)
  
  # # If the default GenzBretz algorithm fails, redo with slower Miwa algorithm
  # 
  # if(likelihood_tmp[1] == 0 & attr(likelihood_tmp, "error") == 0){
  #   
  #   likelihood_tmp <- pmvnorm(mean = mu,
  #                             sigma = sigma,
  #                             lower = lower,
  #                             upper = upper,
  #                             algorithm = "Miwa")
  #   
  # }
  # 
  # if(likelihood_tmp[1] != 0){
  #   
  #   likelihood <- likelihood_tmp[1]
  #   
  # }
  # 
  # if(likelihood_tmp[1] == 0){
  #   
  #   likelihood <- attr(likelihood_tmp, "error")
  #   
  # }

  SSDM_log_lik_vector[i] <- log(likelihood_tmp)
  
}

## #Save to file

filename <- sprintf("outputs/likelihood/SSDM_prob_%s_fold%s_joint_likelihood.rds",
                    dataset_id,
                    fold_id)

saveRDS(SSDM_log_lik_vector,
        filename)


# # SSDM - binary
# 
# ## Vector to store likelihoods
# 
# SSDM_log_lik_array <- array(NA,
#                              dim = c(dim(SSDM_predictions_bin)[1],
#                                      dim(SSDM_predictions_bin)[3]))
# 
# ## Calculate likelihoods
# 
# for(i in seq_len(dim(SSDM_predictions_bin)[1])){
#   
#   for(j in seq_len(dim(SSDM_predictions_bin)[3])){
#     
#   ## Get predicted species probabilities for site
#   
#   pred_probs <- SSDM_predictions_bin[i, , j]
#   
#   ## Get observed species state for site
#   
#   obs_spp <- y[i, ]
#   
#   ## Calculate log likelihood for site
#   
#   log_lik <- sum(dbinom(x = unlist(obs_spp),
#                         size = 1,
#                         prob = pred_probs,
#                         log = TRUE))
#   
#   ## Add value to vector
#   
#   SSDM_log_lik_array[i, j] <- log_lik
#   
#   }
# }
# 
# ## Obtain single likelihood value for whole model
# 
# SSDM_log_likelihood_estimate <- sum(SSDM_log_lik_array,
#                                     na.rm = TRUE)
# 
# ## Save to file
# 
# filename <- sprintf("outputs/likelihood/SSDM_bin_%s_fold%s_likelihood.rds",
#                     dataset_id,
#                     fold_id)
# 
# saveRDS(SSDM_log_likelihood_estimate,
#         filename)
# 
# # SESAM
# 
# ## Vector to store likelihoods
# 
# SESAM_log_lik_vector <- vector(length = dim(SESAM_predictions)[1])
# 
# ## Calculate likelihoods
# 
# for(i in seq_len(dim(SESAM_predictions)[1])){
#   
#   ## Get predicted species probabilities for site
#   
#   pred_probs <- SESAM_predictions[i, , 1]
#   
#   ## Get observed species state for site
#   
#   obs_spp <- y[i, ]
#   
#   ## Calculate log likelihood for site
#   
#   log_lik <- sum(dbinom(x = unlist(obs_spp),
#                         size = 1,
#                         prob = unlist(pred_probs),
#                         log = TRUE))
#   
#   ## Add value to vector
#   
#   SESAM_log_lik_vector[i] <- log_lik
#   
# }
# 
# # ## Obtain single likelihood value for whole model
# # 
# # SESAM_log_likelihood_estimate <- sum(SESAM_log_lik_vector,
# #                                      na.rm = TRUE)
# 
# ## Save to file
# 
# filename <- sprintf("outputs/likelihood/SESAM_%s_fold%s_likelihood.rds",
#                     dataset_id,
#                     fold_id)
# 
# saveRDS(SESAM_log_lik_vector,
#         filename)
