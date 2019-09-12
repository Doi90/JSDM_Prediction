MEM_predict <- function(mem = NULL,
                        model = NULL){
  
  if(is.null(mem)){
    stop("MEM not supplied")
  }
  
  if(is.null(model)){
    stop("Model not supplied")
  }
  
  ## Packages
  
  require(AICcmodavg)
  
  ## Source new prediction function
  
  source("scripts/prediction/new_predictSE.lme.R")
  
  ## Data from model
  
  old_df <- mem$data
  
  ## Create base prediction new data
  
  if("species" %in% colnames(old_df)){
    
    new_df <- expand.grid(dataset = factor(unique(old_df$dataset)),
                          fold = factor(1:5),
                          species = factor(unique(old_df$species)),
                          model = factor(unique(old_df$model)),
                          MPR = 0,
                          HPR = 0,
                          LPR = 0,
                          DPR = 0,
                          HLR_NS = 0,
                          HLR_S = 0,
                          SESAM = 0)
    
  }
  
  if("site" %in% colnames(old_df)){
    
    new_df <- expand.grid(dataset = factor(unique(old_df$dataset)),
                          fold = factor(1:5),
                          site = factor(unique(old_df$site)),
                          model = factor(unique(old_df$model)),
                          MPR = 0,
                          HPR = 0,
                          LPR = 0,
                          DPR = 0,
                          HLR_NS = 0,
                          HLR_S = 0,
                          SESAM = 0)
    
  }
  
  ## Change base prediction data to suit model prediction
  
  if(model == "SSDM"){
    
    ## Leave unchanged as SSDM is the intercept
    
  }
  
  if(model == "MPR"){
    
    new_df$MPR <- 1
    
  }
  
  if(model == "HPR"){
    
    new_df$HPR <- 1
    
  }
  
  if(model == "LPR"){
    
    new_df$LPR <- 1
    
  }
  
  if(model == "DPR"){
    
    new_df$DPR <- 1
    
  }
  
  if(model == "HLR_NS"){
    
    new_df$HLR_NS <- 1
    
  }
  
  if(model == "HLR_S"){
    
    new_df$HLR_S <- 1
    
  }
  
  if(model == "SESAM"){
    
    new_df$SESAM <- 1
    
  }
  
  ## Perform prediction
  
  preds <- new_predictSE.lme(mod = mem,
                             newdata = new_df,
                             se.fit = TRUE,
                             type = "response",
                             print.matrix = TRUE)
  
  ## Extract mean value of model response
  
  sim <- function(preds){
    
    sd(rnorm(n = nrow(preds), 
             mean = preds[ , 1], 
             sd = preds[ , 2]))
    
  }
  
  p_mean <- mean(preds[ , 1])
  
  p_se <- mean(replicate(10, sim(preds)))
  
  #p_se <- sqrt(sum(preds[ , 2] ^ 2) * (1 / nrow(preds)) ^ 2)
  
  ## Return values
  
  return(list(mean = p_mean,
              se = p_se))
  
}