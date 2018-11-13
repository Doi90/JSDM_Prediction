#############################################################
#############################################################
###                                                       ### 
###              DEFINE LIKELIHOOD CODE                   ###
###                                                       ###
### This script contains the function for calculating the ###
### log likelihood value of the models.                   ###
###                                                       ###
#############################################################
#############################################################

log_likelihood <- function(Beta = NULL,
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
  
  if(is.null(n_iter)){
    stop("n_iter not supplied.")
  }  
  
  ## Create an array of distribution mean values. Beta * X values
  
  mean_values <- array(data = NA,
                       dim = c(n_sites,
                               n_species,
                               n_iter))
  
  for(i in seq_len(n_sites)){
    
    for(j in seq_len(n_species)){
      
      for(s in seq_len(n_iter)){
        
        mean_values[i, j, s] <- sum(X[i, ] * Beta[ , j, s])
        
      }
    }
  }
  
  ## Create a log_likelihood matrix full of NAs
  
  log_lik <- matrix(NA,
                    nrow = n_sites,
                    ncol = n_iter)
  
  ## Calculate log likelihood values. Fill matrix with values as we go
  
  approx_counter <- 0
  
  ### For each slice of array
  
  for(s in seq_len(n_iter)){
    
    ### For each site
    
    for(i in seq_len(n_sites)){
      
      occ_state <- y[i, ]         # observed occurrence state  at site i
      
      #### Define probability distribution thresholds
      
      ## lower / upper to limit integral of density for lielihood
      
      lower <- rep(-Inf, n_species)  # default vector of -Inf lower limits
      upper <- rep(+Inf, n_species)  # default vector of +Inf upper limits
      
      for(k in seq_len(n_species)){  # set actual lower/upper limits based on known occurrence states
        
        if(occ_state[k] == 0){     # if species is absent
            
          upper[k] <- 0            # species absent when z<0
        
          } 
          
        if(occ_state[k] == 1){     # if species is present
            
          lower[k] <- 0            # species present when z>0
          
        } 
        
      } 

      #### Prediction for species assemblage at site i using values from slice a
      
      likelihood_tmp <- pmvnorm(mean = colSums(mean_values[ , , s]),
                                sigma = R[ , , s],
                                lower = lower,
                                upper = upper)
      
      if(likelihood_tmp[1] != 0){
        
        likelihood <- likelihood_tmp[1]
        
      }
      
      if(likelihood_tmp[1] == 0){
        
        likelihood <- attr(likelihood_tmp, "error")
        
        approx_counter <- approx_counter + 1
        
      }
      
      #### Fill predictions array with value
      
      log_lik[i, s] <- log(likelihood)
      
    } 
  } 
  
  ## Calculate single likelihood value for whole model
  
  ### Take the product of site-level likelihoods within a single sample
  
  sum_log_lik <- colSums(log_lik,
                          na.rm = TRUE)
  
  ### Take the mean likelihood across samples
  
  mean_log_lik <- mean(sum_log_lik,
                       na.rm = TRUE)
  
  ## Generate single output
  
  log_lik_out <- list(likelihood = mean_log_lik,
                      approximation_rate = approx_counter / (n_sites * n_iter))
  
  return(log_lik_out)
  
}  
