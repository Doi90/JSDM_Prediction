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

species_richness <- function(observed,
                             predicted){
  
  ## Create empty array to store results
  
  metric_names <- c("Observed",
                    "Chao",
                    "Chao_se",
                    "Jackknife_1",
                    "Jackknife_1_se",
                    "Jackknife_2",
                    "Bootstrap",
                    "Bootstrap_se")
  
  spp_richness <- array(NA,
                        dim = c(nrow(observed),
                                length(metric_names),
                                dim(predictions)[3]),
                        dimnames = list(NULL,
                                        metric_names,
                                        NULL))
  
  
  
}