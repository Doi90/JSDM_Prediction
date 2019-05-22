######################################################
######################################################
###                                                ###
###       Create Siulated Dataset Folds            ###
###                                                ###
###  This script generates the test/train splits   ###
### for random and spatial cross-validation folds. ###
### Needs to be done manually before model fitting ###
### to ensure they are identical between model     ###
### runs on Spartan.                               ###
###                                                ###
######################################################
######################################################

#####################
### Load Packages ###
#####################

library(caret)
library(blockCV)
library(sp)
library(sf)

###############################
### Create Custom Functions ###
###############################

## blockCV_multi_spp_test() checks the generated folds
## meet our requirements. Training data in each CV split
## must have at least one presence for all species.

blockCV_multi_spp_test <- function(dataset,
                                   blockCV,
                                   n_sites){
  
  train_tests <- vector(length = 5)
  
  for(i in seq_len(5)){
    
    train_id <- blockCV$folds[[i]][[1]]
    
    y_train <- dataset[train_id, ]
    
    train_tests[i] <- all(colMax(y_train) > 0)
    
  }
  
  test_sites <- c(blockCV$folds[[1]][[2]],
                  blockCV$folds[[2]][[2]],
                  blockCV$folds[[3]][[2]],
                  blockCV$folds[[4]][[2]],
                  blockCV$folds[[5]][[2]])
  
  
  if(all(train_tests == TRUE) & length(unique(test_sites)) == n_sites){
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
  }
  
}

#####################
### Create Splits ###
#####################

for(i in seq_len(10)){
  
  for(j in c("random", "spatial")){
    
    message(sprintf("%s_%s",
                    i, j))

    ## Load data
    
    filename <- sprintf("data/sim%1$s%2$s/sim%1$s%2$s_communitySimul.rds",
                        i,
                        j)
    
    data <- readRDS(filename)
    
    filename <- sprintf("data/sim%1$s%2$s/sim%1$s%2$s_coordinates.rds",
                        i,
                        j)
    
    data_coords <- readRDS(filename)
    
    ## Manipulate data
    
    data_Y <- data$data$Y
    
    data_X <- data$data$X
    
    data_LL <- data_coords
    
    ## Generate folds
    
    ### Define control parameters
    
    dataset_ok <- FALSE
    
    iteration <- 0
    
    ### Generate folds until meet our test condition
    
    while(dataset_ok == FALSE){
      
      ## Increase counter
      
      iteration <- iteration + 1
      
      ## Prevent infinite loop
      
      if(iteration >= 1000){
        
        break()
        
      }
      
      ## Generate folds
      
      set.seed(28041948 + iteration)                    # Creator's Birthday
      
      data_SP <- st_as_sf(as.data.frame(data_LL), coords = 1:2)
      
      data_blockCV <- spatialBlock(speciesData = st_jitter(data_SP, 0.05),
                                   rows = 5,
                                   cols = 5,
                                   k = 5,
                                   iteration = 100)
      
      dataset_ok <- blockCV_multi_spp_test(dataset = data_Y,
                                           blockCV = data_blockCV,
                                           n_sites = nrow(data_Y))
      
    }
    
    ## Write to file
    
    for(k in seq_len(5)){
      
      fold_id <- paste0("fold", k, collapse = "")          # set fold id
      
      train_id <- data_blockCV$folds[[k]][[1]]            # Training site ids
      test_id <- data_blockCV$folds[[k]][[2]]             # Testing site ids
      
      tmp_y_test <- data_Y[test_id, ]                     # y test data
      tmp_y_train <- data_Y[train_id, ]                   # y train data
      
      tmp_X_test <- data_X[test_id, -1]                     # X test data
      tmp_X_train <- data_X[train_id, -1]                   # X train data
      
      tmp_LL_test <- data_LL[test_id, ]                   # Lat/long test data
      tmp_LL_train <- data_LL[train_id, ]                 # Lat/long train data
      
      site_ids <- list(data_blockCV$folds[[1]][[2]],
                       data_blockCV$folds[[2]][[2]],
                       data_blockCV$folds[[3]][[2]],
                       data_blockCV$folds[[4]][[2]],
                       data_blockCV$folds[[5]][[2]])
      
      y_test_filename <- sprintf('data/sim%1$s%2$s/y_sim%1$s%2$s_%3$s_test.csv',
                                 i, j, fold_id)     # test filename
      y_train_filename <- sprintf('data/sim%1$s%2$s/y_sim%1$s%2$s_%3$s_train.csv',
                                  i, j, fold_id)   # train filename
      
      X_test_filename <- sprintf('data/sim%1$s%2$s/X_sim%1$s%2$s_%3$s_test.csv',
                                 i, j, fold_id)     # test filename
      X_train_filename <- sprintf('data/sim%1$s%2$s/X_sim%1$s%2$s_%3$s_train.csv',
                                  i, j,fold_id)   # train filename
      
      LL_test_filename <- sprintf('data/sim%1$s%2$s/LL_sim%1$s%2$s_%3$s_test.csv',
                                  i, j, fold_id)   # test filename
      LL_train_filename <- sprintf('data/sim%1$s%2$s/LL_sim%1$s%2$s_%3$s_train.csv',
                                   i, j, fold_id) # train filename
      
      site_id_filename <- sprintf("data/sim%1$s%2$s/site_ids.rds",
                                  i, j)
      
      write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
      write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
      
      write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
      write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
      
      write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
      write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
      
      saveRDS(site_ids,
              site_id_filename)
      
    }
  }
}


