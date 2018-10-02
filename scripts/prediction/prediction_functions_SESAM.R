##########################################################
##########################################################
###                                                    ###
###               DEFINE SESAM FUNCTION                ###
###                                                    ###
###    This script defines the SESAM function used in  ###
### this analysis. It takes predicted probabilities of ###
### presence for all species/sites and predicted       ###
### estimates (here, sum of probabilities at sites)    ###
### and returns binary predictions capped to species   ###
### richness estimates.                                ###
###                                                    ###
##########################################################
##########################################################

SESAM <- function(probabilities,
                  species_richness){
  
  ## Round species richness estimates to whole numbers
  
  projected_SR <- round(as.vector(species_richness))
  
  ## Create new probabilities object
  
  new_prediction <- probabilities
  
  ## For each site
  
  for(i in seq_len(nrow(probabilities))){
    
    ## Site-specific species richness estimate
    
    SR <- projected_SR[i]
    
    if(SR > 0){
      
      ## Site-specific species probabilities
      
      predicted_composition <- probabilities[i, ]
      
      ## Get order of probabilities
      
      order <- order(predicted_composition,
                     decreasing = TRUE)
      
      ## Only most probable species (up to SR limit) are predicted
      
      present_species <- order[1:SR]
      
      predicted_composition[, present_species] <- 1
      
      predicted_composition[, -present_species] <- 0
      
    } else {
      
      predicted_composition[, ] <- 0
      
    }
    
    new_prediction[i, ] <- predicted_composition
    
  }
  
  return(new_prediction)
  
}
