###################################################
###################################################
###                                             ###
###     EXTRACT SITE IDS FROM EXISTING DATA     ###
###                                             ###
###   This script extracts the site ids from    ###
### already generated dataet folds.             ###
###                                             ###
###################################################
###################################################

#####################
### Load Packages ###
#####################

library(PassButter)

###############
### Options ###
###############

dataset_options <- c("frog",
                     "eucalypt",
                     #"bird",
                     "sim1random",
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
                     "sim10spatial")

####################################################
### Extract Site IDs from existing dataset folds ###
####################################################

for(dataset in dataset_options){
  
  if(dataset %in% c("frog", "eucalypt")){
    
    df_1 <- read.csv(sprintf("data/%1$s/X_%1$s_fold1_test.csv", dataset))
    df_2 <- read.csv(sprintf("data/%1$s/X_%1$s_fold2_test.csv", dataset))
    df_3 <- read.csv(sprintf("data/%1$s/X_%1$s_fold3_test.csv", dataset))
    df_4 <- read.csv(sprintf("data/%1$s/X_%1$s_fold4_test.csv", dataset))
    df_5 <- read.csv(sprintf("data/%1$s/X_%1$s_fold5_test.csv", dataset))
    
    site_ids <- list(df_1[ , "X"],
                     df_2[ , "X"],
                     df_3[ , "X"],
                     df_4[ , "X"],
                     df_5[ , "X"])
    
    saveRDS(site_ids,
            sprintf("data/%1$s/site_ids.rds", dataset))
  }
  
  if(dataset == "bird"){
    
    bird <- read.csv("data/bird/Birds_Compiled.csv")[ , 371:378]
    bird <- as.matrix(scale(bird))
    
    df_1 <- read.csv(sprintf("data/%1$s/X_%1$s_fold1_test.csv", dataset))[, -1]
    df_2 <- read.csv(sprintf("data/%1$s/X_%1$s_fold2_test.csv", dataset))[, -1]
    df_3 <- read.csv(sprintf("data/%1$s/X_%1$s_fold3_test.csv", dataset))[, -1]
    df_4 <- read.csv(sprintf("data/%1$s/X_%1$s_fold4_test.csv", dataset))[, -1]
    df_5 <- read.csv(sprintf("data/%1$s/X_%1$s_fold5_test.csv", dataset))[, -1]
    
    site_ids <- vector("list", 5)
    
    for(i in 1:5){
      
      command <- sprintf("tmp <- df_%s", i)
      
      eval(parse(text = command))
      
      tmp <- as.matrix(tmp)
      
      counter <- 0
      
      for(j in seq_len(nrow(tmp))){
        
        for(k in seq_len(nrow(bird))){
          
          test <- all.equal(bird[k, ], tmp[j, ])
          
          if(isTRUE(test)){
            
            site_ids[[i]] <- c(site_ids[[i]], k)
            
          }
        
          counter <- counter + 1
          
          statusbar(run = counter,
                    max.run = (5 * nrow(df_1) * nrow(bird)),
                    percent.max = 20)
          }
      }
    }
    
    saveRDS(site_ids,
            sprintf("data/%1$s/site_ids.rds", dataset))
    
  }
  
  if(isTRUE(grep("sim", dataset) == 1)){
    
    df_1 <- read.csv(sprintf("data/%1$s/X_%1$s_fold1_test.csv", dataset))
    df_2 <- read.csv(sprintf("data/%1$s/X_%1$s_fold2_test.csv", dataset))
    df_3 <- read.csv(sprintf("data/%1$s/X_%1$s_fold3_test.csv", dataset))
    df_4 <- read.csv(sprintf("data/%1$s/X_%1$s_fold4_test.csv", dataset))
    df_5 <- read.csv(sprintf("data/%1$s/X_%1$s_fold5_test.csv", dataset))
    
    site_ids <- list(gsub("site", "", df_1$X),
                     gsub("site", "", df_2$X),
                     gsub("site", "", df_3$X),
                     gsub("site", "", df_4$X),
                     gsub("site", "", df_5$X))
    
    saveRDS(site_ids,
            sprintf("data/%1$s/site_ids.rds", dataset))
  }
  
}
