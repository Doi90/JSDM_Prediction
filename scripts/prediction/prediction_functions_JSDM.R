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
# TWO MARGINAL PROBABILITY FUNCTIONS
# predict.marginal.probability returns probabilities. Integrates over distribution
# predict.marginal.binary returns 1/0s. Samples distribution

predict.marginal.probability <- function(Beta = NULL,
                                         X = NULL,
                                         n_species = NULL,
                                         n_sites = NULL,
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
  
  if(is.null(n_sites)){
    stop("n_sites not supplied.")
  }
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  }
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(n_sites,     # Number of sites in test data
                               n_species,   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       colnames(Beta),
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  for(i in seq_len(n_sites)){
    
    for(j in seq_len(n_species)){
      
      for(s in seq_len(n_iter)){
        
        predictions[i, j, s] <- pnorm(sum(X[i, ] * Beta[ , j, s]))
        
      }
    }
  }
  
  
  if(any(is.na(predictions))){
    warning("Some predictions returned NAs")
  }
  
  return(predictions)
  
}

predict.marginal.binary <- function(Beta = NULL,
                                    X = NULL,
                                    n_species = NULL,
                                    n_sites = NULL,
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
  
  if(is.null(n_sites)){
    stop("n_sites not supplied.")
  }
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  }
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(n_sites,     # Number of sites in test data
                               n_species,   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       colnames(Beta),
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  for(i in seq_len(n_sites)){
    
    for(j in seq_len(n_species)){
      
      for(s in seq_len(n_iter)){
        
        predictions[i, j, s] <- ifelse(rnorm(n = 1,
                                             mean = sum(X[i, ] * Beta[ , j, s]),
                                             sd = 1) > 0,
                                       yes = 1,
                                       no = 0)
        
      }
    }
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
                                    n_sites = NULL,
                                    n_iter = NULL,
                                    dataset_id = NULL){
  
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
  
  if(is.null(n_sites)){
    stop("n_sites not supplied.")
  } 
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  } 
  
  if(is.null(dataset_id)){
    stop("dataset_id not supplied.")
  }
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- array(data = NA,
                       dim = c(n_sites,
                               n_species,
                               n_iter))
  
  for(s in seq_len(n_iter)){
    
    mean_values[, , s] <- X %*% Beta[ , , s]
    
  }
  
  ## Create a prediction array full of NAs. Unlike other predictions,
  ## this one needs a 4D array for predictions as we need to predict
  ## all species at a site once for three different species left in
  ## scenarios. Thus, 3 predictions per species per site
  
  predictions <- array(NA,
                       dim = c(n_sites,     # Number of sites in test data
                               n_species,
                               n_iter,      # 1:1 Prediction slice:Posterior slice
                               3),          # Need to predict for each species left in   
                       dimnames = list(rownames(X),
                                       colnames(y),
                                       NULL,
                                       NULL))
  
  
  ## Identify the species being left in
  
  if(dataset_id == "frog"){
    
    species_left_in_IDs <- c(9, 4, 6) # Lit_rani, Lim_tas, Lit_ewing
    
  }
  
  if(dataset_id == "eucalypt"){
    
    species_left_in_IDs <- c(8, 9, 3) # OVA, WIL, BAX
    
  }
  
  if(dataset_id == "bird"){
    
    species_left_in_IDs <- c(16, 19, 289) # Common_Goldeneye, Canada_Goose, American_Robin
    
  }
  
  if(dataset_id == "butterfly"){
    
    species_left_in_IDs <- c(5, 6, 43) # high_brown_fritillary, dark.green_fritillary, common_blue
    
  }
  
  if(dataset_id %in% c("sim1random",
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
                       "sim10spatial")){
    
    spp_LOI <- readRDS("data/simulated_datasets_species_left_in.rds")
    
    species_left_in_IDs <- as.numeric(spp_LOI[dataset_id, ])
    
  }
  ## Make predictions. Fill predictions array with values as we go
  
  ### For each 4th dimension/left-in species
  
  for(j in species_left_in_IDs){
    
    j_array_id <- which(species_left_in_IDs == j)
    
    ### For each slice of array
    
    for(s in seq_len(n_iter)){
      
      ### For each site
      
      for(i in seq_len(n_sites)){
        
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
                             mean = mean_values[ i, , s],
                             sigma = R[ , , s],
                             lower = lower, 
                             upper = upper,
                             algorithm = "gibbs")
        
        #### Convert latent variable draws to pres/abs
        
        spp_pred <- ifelse(spp_pred > 0, 1, 0) 
        
        #### Fill predictions array with value
        
        predictions[i, , s, j_array_id] <- spp_pred
        
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
                                    n_sites = NULL,
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
  
  if(is.null(n_sites)){
    stop("n_sites not supplied.")
  }
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  } 
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- array(data = NA,
                       dim = c(n_sites,
                               n_species,
                               n_iter))
  
  for(s in seq_len(n_iter)){
    
    mean_values[, , s] <- X %*% Beta[ , , s]
    
  }
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(n_sites,     # Number of sites in test data
                               n_species,   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       colnames(y),
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  ### For each slice of array
  
  for(s in seq_len(n_iter)){
    
    ### For each site
    
    for(i in seq_len(n_sites)){
      
      occ_state <- y[i, ]         # observed occurrence state  at site i
      
      ### For each species
      
      for(j in seq_len(n_species)){      
        
        #### Define probability distribution thresholds
        
        ## lower / upper to truncate distribution (known pres/abs)
        
        lower <- rep(-Inf, n_species)  # default vector of -Inf lower limits
        upper <- rep(+Inf, n_species)  # default vector of +Inf upper limits
        
        for(jj in seq_len(n_species)){  # set actual lower/upper limits based on known occurrence states
          
          if(jj != j){                  # do this for all species we aren't predicting
            
            if(occ_state[jj] == 0){     # if species is absent
              upper[jj] <- 0            # species absent when z<0
            } 
            
            if(occ_state[jj] == 1){     # if species is present
              lower[jj] <- 0            # species present when z>0
            } 
            
          } 
        }
        
        ## lowerx / upperx to define probability to calculate (target species Pr(z>0))
        
        lowerx <- lower
        upperx <- upper
        
        lowerx[j] <- 0
        
        #### Prediction for species j at site i using values from slice a
        
        spp_pred <- ptmvtnorm.new(mean = mean_values[ i, , s],
                                  sigma = R[ , , s],
                                  lower = lower,
                                  upper = upper,
                                  lowerx = lowerx,
                                  upperx = upperx,
                                  algorithm = "GenzBretz")
        
        ## If default algorithm returns NA, use slower Miwa algorithm
        
        if(is.na(spp_pred)){
          
          spp_pred <- ptmvtnorm.new(mean = mean_values[ i, , s],
                                    sigma = R[ , , s],
                                    lower = lower,
                                    upper = upper,
                                    lowerx = lowerx,
                                    upperx = upperx,
                                    algorithm = "Miwa")
          
        }
        
        #### Fill predictions array with value
        
        predictions[i, j, s] <- spp_pred[[1]]
        
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
                          n_sites = NULL,
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
  
  if(is.null(n_sites)){
    stop("n_sites not supplied.")
  } 
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  }  
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- array(data = NA,
                       dim = c(n_sites,
                               n_species,
                               n_iter))
  
  for(s in seq_len(n_iter)){
    
    mean_values[, , s] <- X %*% Beta[ , , s]
    
  }
  
  ## Create a prediction array full of NAs
  
  predictions <- array(NA,
                       dim = c(n_sites,     # Number of sites in test data
                               n_species,   # Number of species in test data   
                               n_iter),     # 1:1 Prediction slice:Posterior slice
                       dimnames = list(rownames(X),
                                       colnames(y),
                                       NULL))
  
  ## Make predictions. Fill predictions array with values as we go
  
  ### For each slice of array
  
  for(s in seq_len(n_iter)){
    
    ### For each site
    
    for(i in seq_len(n_sites)){
      
      occ_state <- y[i, ]         # observed occurrence state  at site i
      
      #### Define probability distribution thresholds
      
      ## lower / upper to truncate distribution (known pres/abs)
      
      lower <- rep(-Inf, n_species)  # default vector of -Inf lower limits
      upper <- rep(+Inf, n_species)  # default vector of +Inf upper limits
      
      #### Prediction for species assemblage at site i using values from slice a
      
      spp_pred <- rtmvnorm(n = 1,
                           mean = mean_values[ i, , s],
                           sigma = R[ , , s], 
                           lower = lower, 
                           upper = upper,
                           algorithm = "gibbs")
      
      #### Convert latent variable draws to pres/abs
      
      spp_pred <- ifelse(spp_pred > 0, 1, 0)
      
      #### Fill predictions array with value
      
      predictions[i, , s] <- spp_pred
      
    } 
  } 
  
  return(predictions)
  
}  