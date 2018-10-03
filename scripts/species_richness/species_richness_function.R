######################################################
######################################################
###                                                ###
###           Species Richness Function            ###
###                                                ###
###  This script calculates the species richness   ###
### of sites                                       ###
###                                                ###
######################################################
######################################################

species_richness <- function(observed = NULL,
                             predicted = NULL){
  
  ## Tests to make sure correct inputs supplied
  
  if(is.null(observed)){
    stop("observed not supplied.")
  } 
  
  if(is.null(predicted)){
    stop("predicted not supplied.")
  }
  
  ## Create empty array to store results
  
  spp_richness <- array(NA,
                        dim = c(nrow(observed),       # sites
                                3,                    # Observed + Predicted + Difference
                                dim(predicted)[3]),   # Samples
                        dimnames = c(NULL,
                                     c("Observed",
                                       "Predicted",
                                       "Difference"),
                                     NULL))   
  
  ## Observed richness
  
  obs_richness <- apply(observed,
                        MARGIN = 1,
                        sum,
                        na.rm = TRUE)
  
  spp_richness[ , 1, ] <- obs_richness
  
  ## Predicted richness
  
  for(i in seq_len(dim(predicted)[3])){
    
    pred_richness <- apply(predicted[ , , i],
                           MARGIN = 1,
                           sum,
                           na.rm = TRUE)
    
    spp_richness[ , 2, i] <- pred_richness
    
  }
  
  ## Richness difference
  
  spp_richness[ , 3, ] <- spp_richness[ , 1, ] - spp_richness[ , 2, ] 
  
  ## Output
  
  return(spp_richness)
  
}