#############################################################
#############################################################
###                                                       ### 
###              DEFINE PREDICTION CODE                   ###
###                                                       ###
### This script defines the different types of prediction ###
### we use in our JSDM prediction comparison.             ###
###                                                       ###
#############################################################
#############################################################

### Marginal Prediction ----

# Marginal prediction involves predicting the probability
# of occurrence of a species WITHOUT reference to the
# other target species.
# Environment-only prediction. Ignores species interactions.

predict.marginal <- function(Beta = NULL,
                             X = NULL,
                             n_species = NULL,
                             n_iter = NULL){
  
  ## Tests to make sure correct inputs supplied
  
  if(is.null(Beta)){
    stop("Beta not supplied.")
  }
  
  if(is.null(X)){
    stop("X not supplied.")
  }
  
  if(is.null(n_species)){
    stop("n_species not supplied.")
  }
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  }
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(nrow(X),     # Number of sites in test data
                               n_species,   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       colnames(Beta),
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  for(i in seq_len(dim(Beta)[3])){
    
    predictions[ , , i] <- pnorm(X * Beta[ , , i])
    
  }
  
  if(any(is.na(predictions))){
    warning("Some predictions returned NAs")
  }
  
  return(predictions)
  
}

### Conditional prediction ----

# Conditional prediction involves predicting the probability
#  of species occurrence conditional on the occurrence state
#  of other target species.
# Can be defined J-1 ways, where J = # species
# For this study we consider leave-one-in and leave-one-out
#  style conditional predictions.

#### Leave-one-in ----

# Leave-one-in style conditional prediction involves predicting
#  the probability of occurrence for J-1 target species when we
#  know the occurrence state of one

predict.conditional.LOI <- function(Beta = NULL,
                                    X = NULL,
                                    y = NULL,
                                    R = NULL,
                                    n_species = NULL,
                                    n_iter = NULL){
  
  ## Tests to make sure correct inputs supplied
  
  if(is.null(Beta)){
    stop("Beta not supplied.")
  }
  
  if(is.null(X)){
    stop("X not supplied.")
  }
  
  if(is.null(y)){
    stop("y not supplied.")
  } 
  
  if(is.null(R)){
    stop("R not supplied.")
  } 
  
  if(is.null(n_species)){
    stop("n_species not supplied.")
  } 
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  } 
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- Beta * array(data = X,
                              dim = c(nrow(X),
                                      ncol(X),
                                      dim(Beta)[3]))
  
  ## Create a prediction array full of NAs. Unlike other predictions,
  ## this one needs a 4D array for predictions as we need to predict
  ## all species at a site once for each species left out. Thus,
  ## J-1 predictions per species per site
  
  predictions <- array(NA,
                       dim = c(nrow(X),     # Number of sites in test data
                               n_species,
                               n_iter,      # 1:1 Prediction slice:Posterior slice
                               n_species),  # Need to predict for each species left in   
                       dimnames = list(rownames(X),
                                       colnames(Beta),
                                       NULL,
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  ### For each 4th dimension/left-in species
  
  for(j in seq_len(n_species)){
    
    ### For each slice of array
    
    for(a in seq_len(dim(Beta)[3])){
      
      ### For each site
      
      for(i in seq_len(nrow(y))){
        
        occ_state <- y[i, ]         # observed occurrence state  at site i
        
        #### Define probability distribution thresholds
        
        ## lower / upper to truncate distribution (known pres/abs)
        
        lower <- rep(-Inf, n_species)  # default vector of -Inf lower limits
        upper <- rep(+Inf, n_species)  # default vector of +Inf upper limits
        
        if(occ_state[j] == 0){         # If species j (left-in) is absent
          upper[j] <- 0                # set upper threshold
        } 
        
        if(occ_state[j] == 1){         # If species j (left-in) is present
          lower[j] <- 0                # set lower threshold
        } 
        
        #### Perform prediction / random draws from multivariate normal  
        
        spp_pred <- rtmvnorm(n = 1,
                             mean = colSums(mean_values[ , , a]),
                             sigma = R[ , , a],
                             lower = lower, 
                             upper = upper)
        
        #### Convert latent variable draws to pres/abs
        
        spp_pred <- ifelse(spp_pred > 0, 1, 0) 
        
        #### Fill predictions array with value
        
        predictions[i, , a, j] <- spp_pred
        
      } 
    }
  }
  
  return(predictions)  
  
}

#### Leave-one-out ----

# Leave-one-out style conditional prediction involves predicting
#  the probability of occurrence of one species when we know the
#  occurrence state of J-1 target species.

predict.conditional.LOO <- function(Beta = NULL,
                                    X = NULL,
                                    y = NULL,
                                    R = NULL,
                                    n_species = NULL,
                                    n_iter = NULL){
  
  ## Tests to make sure correct inputs supplied
  
  if(is.null(Beta)){
    stop("Beta not supplied.")
  }
  
  if(is.null(X)){
    stop("X not supplied.")
  }
  
  if(is.null(y)){
    stop("y not supplied.")
  } 
  
  if(is.null(R)){
    stop("R not supplied.")
  } 
  
  if(is.null(n_species)){
    stop("n_species not supplied.")
  } 
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  } 
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- Beta * array(data = X,
                              dim = c(nrow(X),
                                      ncol(X),
                                      dim(Beta)[3]))
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(nrow(X),     # Number of sites in test data
                               n_species,   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       colnames(Beta),
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  ### For each slice of array
  
  for(a in seq_len(dim(Beta)[3])){
    
    ### For each site
    
    for(i in seq_len(nrow(y))){
      
      occ_state <- y[i, ]         # observed occurrence state  at site i
      
      ### For each species
      
      for(j in seq_len(n_species)){      
        
        #### Define probability distribution thresholds
        
        ## lower / upper to truncate distribution (known pres/abs)
        
        lower <- rep(-Inf, n_species)  # default vector of -Inf lower limits
        upper <- rep(+Inf, n_species)  # default vector of +Inf upper limits
        
        for(k in seq_len(n_species)){  # set actual lower/upper limits based on known occurrence states
          
          if(k != j){                  # do this for all species we aren't predicting
            
            if(occ_state[k] == 0){     # if species is absent
              upper[k] <- 0            # species absent when z<0
            } 
            
            if(occ_state[k] == 1){     # if species is present
              lower[k] <- 0            # species present when z>0
            } 
            
          } 
        }
        
        ## lowerx / upperx to define probability to calculate (target species Pr(z>0))
        
        lowerx <- lower
        upperx <- upper
        
        lowerx[j] <- 0
        
        #### Prediction for species j at site i using values from slice a
        
        spp_pred <- ptmvnorm(mean = colSums(mean_values[ , , a]),
                             sigma = R[ , , a],
                             lower = lower,
                             upper = upper,
                             lowerx = lowerx,
                             upperx = upperx)
        
        #### Fill predictions array with value
        
        predictions[i, j, a] <- spp_pred
        
      } 
    } 
  } 
  
  return(predictions)
  
} 

### Joint prediction ----

# Joint prediction involves predicting the probability of 
#  occurrence of all species while simultaneously accounting
#  for both the measured environmental variables and species
#  interactions

predict.joint <- function(Beta = NULL,
                          X = NULL,
                          y = NULL,
                          R = NULL,
                          n_species = NULL,
                          n_iter = NULL){
  
  ## Tests to make sure correct inputs supplied
  
  if(is.null(Beta)){
    stop("Beta not supplied.")
  } 
  
  if(is.null(X)){
    stop("X not supplied.")
  } 
  
  if(is.null(y)){
    stop("y not supplied.")
  }  
  
  if(is.null(R)){
    stop("R not supplied.")
  }  
  
  if(is.null(n_species)){
    stop("n_species not supplied.")
  }  
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  }  
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- Beta * array(data = X,
                              dim = c(nrow(X),
                                      ncol(X),
                                      dim(Beta)[3]))
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(nrow(X),     # Number of sites in test data
                               n_species,   # Number of species in test data   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       "Prediction",
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  ### For each slice of array
  
  for(a in seq_len(dim(Beta)[3])){
    
    ### For each site
    
    for(i in seq_len(nrow(y))){
      
      occ_state <- y[i, ]         # observed occurrence state  at site i
      
      #### Define probability distribution thresholds
      
      ## lower / upper to truncate distribution (known pres/abs)
      
      lower <- rep(-Inf, n_species)  # default vector of -Inf lower limits
      upper <- rep(+Inf, n_species)  # default vector of +Inf upper limits
      
      #### Prediction for species assemblage at site i using values from slice a
      
      spp_pred <- rtmvnorm(n = 1,
                           mean = colSums(mean_values[ , , a]),
                           sigma = R[ , , a], 
                           lower = lower, 
                           upper = upper)
      
      #### Convert latent variable draws to pres/abs
      
      spp_pred <- ifelse(spp_pred > 0, 1, 0)
      
      #### Fill predictions array with value
      
      predictions[i, , a] <- spp_pred
      
    } 
  } 
  
  return(predictions)
  
}  