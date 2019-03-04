############################################
############################################
###                                      ###
###          SIMULATE JSDM DATA          ###
###                                      ###
###   This script will simulate datasets ###
### for use in the JSDM prediction       ###
### comparison. Prediction is partially  ###
### manual and partially within the HMSC ###
### packages communitySimul() function.  ###
###                                      ###
############################################
############################################

#####################
### Load Packages ###
#####################

library(HMSC)
library(blockCV)

########################
### Define Constants ###
########################

n_species <- 10

n_sites <- 100

n_iter <- 10

#########################################
### Set Latitude/Longitude Boundaries ###
#########################################

latitude_range <- c(-10,10)

longitude_range <- c(-10,10)

########################
### Custom FUnctions ###
########################

rcoef <- function(nsp,
                  ncoef = 4,
                  variance = 1){
  
  matrix(rnorm(nsp * ncoef,
               0,
               variance),
         ncoef,nsp)
  
}

rcoef_hier <- function(nsp,
                       ncoef = 4,
                       variance = 1,
                       hierarchical_variance = 1){
  
  means <- rnorm(ncoef, 0, hierarchical_variance)
  
  baseline <- rcoef(nsp, ncoef, variance)
  
  sweep(baseline, 1, means, "+")
  
}

rlv <- function(nsp,
                nlv = 3){
  
  Lambda <- matrix(rnorm(nsp * nlv), nsp, nlv)
  
  Sigma <- Lambda %*% t(Lambda) + diag(nsp)
  
  cov2cor(Sigma)
  
}

################
### Simulate ###
################

for(i in seq_len(n_iter)){
  
  for(j in c("random", "spatial")){
    
    set.seed <- 28041948 + 2*i
    
    dataset_ok <- FALSE
    
    while(dataset_ok == FALSE){
      
      ###################################
      ### Simulate Latitude/Longitude ###
      ###################################
      
      latlon <- cbind(runif(n = n_sites,
                            min = latitude_range[1],
                            max = latitude_range[2]),
                      runif(n = n_sites,
                            min = longitude_range[1],
                            max = longitude_range[2]))
      
      ###################################
      ### Simulate Measured Variables ###
      ###################################
      
      meas_var <- cbind(1,
                        scale(rnorm(n = n_sites,
                                    mean = sample(seq(1, 5, .1),
                                                  size = 1,
                                                  replace = TRUE),
                                    sd = 1)),
                        scale(rnorm(n = n_sites,
                                    mean = sample(seq(1, 5, .1),
                                                  size = 1,
                                                  replace = TRUE),
                                    sd = 1)),
                        scale(rnorm(n = n_sites,
                                    mean = sample(seq(-5, -1, .1),
                                                  size = 1,
                                                  replace = TRUE),
                                    sd = 1)),
                        rbinom(n = n_sites,
                               size = 1,
                               prob = sample(seq(0.1, 0.5, .1),
                                             size = 1,
                                             replace = TRUE)),
                        rbinom(n = n_sites,
                               size = 1,
                               prob = sample(seq(0.5, 0.9, .1),
                                             size = 1,
                                             replace = TRUE)))
      
      ##################################
      ### Simulate Species Responses ###
      ##################################
      
      spp_resp <- rcoef_hier(nsp = n_species,
                             ncoef = 6,
                             variance = 1,
                             hierarchical_variance = 1)
      
      spp_resp <- t(spp_resp)  
      
      #######################################
      ### Simulate Random/Spatial Effects ###
      #######################################
      
      ## Random
      
      random <- data.frame(as.factor(seq_len(n_sites)))
      
      ## Spatial
      
      LL <- as.data.frame(latlon)
      
      LF_code <- seq_len(nrow(LL))
      
      spatial <- data.frame(level = as.factor(LF_code),
                            LL)
      
      ###################################################
      ### Simulate LF Approximated Correlation Matrix ###
      ###################################################
      
      corr <- rlv(nsp = n_species,
                  nlv = 3)
      
      ##########################
      ### Simulate Community ###
      ##########################
      
      message(sprintf("%s_%s", i, j))
      
      if(j == "random"){
        
        data_random <- communitySimul(X = meas_var,
                                      paramX = spp_resp,
                                      Random = random,
                                      nsp = n_species,
                                      family = "probit")
        
        if(all(colSums(data_random$data$Y) >= 5 &
               colSums(data_random$data$Y) <= 95)){
          
          dataset_ok <- TRUE
          
        }
      }
      
      if(j == "spatial"){
        
        data_spatial <- communitySimul(X = meas_var,
                                       paramX = spp_resp,
                                       Auto = spatial,
                                       nsp = n_species,
                                       family = "probit")
        
        if(all(colSums(data_random$data$Y) >= 5 &
               colSums(data_random$data$Y) <= 90)){
          
          dataset_ok <- TRUE
          
        }
      }
      
    }
    
    ##########################
    ### Write Data To File ###
    ##########################
    
    ## Create Directory
    
    dir_name <- sprintf("data/sim%s%s",
                        i,
                        j)
    
    if(!dir.exists(dir_name)){
      
      dir.create(dir_name)
      
    }
    
    ## Write to file
    
    ### Random
    
    if(j == "random"){
      
      filename <- sprintf("data/sim%1$s%2$s/sim%1$srandom_communitySimul.rds",
                          i,
                          j)
      
      saveRDS(object = data_random,
              file = filename)
      
    }
    
    ### Spatial
    
    if(j == "spatial"){
      
      filename <- sprintf("data/sim%1$s%2$s/sim%1$sspatial_communitySimul.rds",
                          i,
                          j)
      
      saveRDS(object = data_spatial,
              file = filename)
      
    }
    
    ## Coordinates
    
    filename <- sprintf("data/sim%1$s%2$s/sim%1$s%2$s_coordinates.rds",
                        i,
                        j)
    
    saveRDS(object = latlon,
            file = filename)
    
  }
}