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

log_likelihood <- function(Beta = Beta,
                           X,
                           y,
                           R,
                           n_species = n_species,
                           n_sites = n_sites,
                           n_iter = n_iter){
  
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
                              dim = c(n_sites,
                                      n_species,
                                      n_iter))
  
  ## Create a log_likelihood matrix full of NAs
  
  log_lik <- matrix(NA,
                    nrow = n_sites,
                    ncol = n_iter)
  
  ## Calculate log likelihood values. Fill matrix with values as we go
  
  ### For each slice of array
  
  for(a in seq_len(n_iter)){
    
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
      
      likelihood_tmp <- pmvnorm(mean = colSums(mean_values[ , , a]),
                                sigma = R[ , , a],
                                lower = lower,
                                upper = upper)
      
      approx_counter <- 0
      
      if(likelihood_tmp[1] != 0){
        
        likelihood <- likelihood_tmp[1]
        
      }
      
      if(likelihood_tmp[1] == 0){
        
        likelihood <- likelihood_tmp[2]
        
        approx_counter <- approx_counter + 1
        
      }
      
      #### Fill predictions array with value
      
      log_lik[i, a] <- log(likelihood)
      
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
                      approximation_rate = approx_counter / n_sites)
  
  return(log_lik_out)
  
}  
