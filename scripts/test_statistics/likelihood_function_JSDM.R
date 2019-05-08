#############################################################
#############################################################
###                                                       ### 
###              DEFINE LIKELIHOOD CODE                   ###
###                                                       ###
### This script contains the different functions for      ###
### calculating the log likelihood value of the models.   ###
### Kept the old version (log_likelihood_old) and Nick's  ###
### replacements for reference                            ###
###                                                       ###
#############################################################
#############################################################

####################################
### Old version no longer in use ###
####################################

log_likelihood_old <- function(Beta = NULL,
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
      
      likelihood_tmp <- pmvnorm(mean = mean_values[ i, , s],
                                sigma = R[ , , s],
                                lower = lower,
                                upper = upper)
      
      # If the default GenzBretz algorithm fails, redo with slower Miwa algorithm
      
      if(likelihood_tmp[1] == 0 & attr(likelihood_tmp, "error") == 0){
        
        likelihood_tmp <- pmvnorm(mean = mean_values[ i, , s],
                                  sigma = R[ , , s],
                                  lower = lower,
                                  upper = upper,
                                  algorithm = "Miwa")
        
      }
      
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

#######################################
### New Version: Multiple functions ###
#######################################

### Nick's four internal functions

logapb <- function(loga = NULL,
                   logb = NULL){
  
  # calculate log(a + b) stably from log(a) and log(b)
  # log(a + b)  = log(a * (1 + b / a))
  #             = log(a) + log(1 + b / a)
  #             = log(a) + log(1 + exp(log(b) - log(a)))
  #             = log(a) + log1p(exp(log(b) - log(a)))
  
  ans <- loga + log1p(exp(logb - loga))
  
  return(ans)
  
}

logmean <- function(lx = NULL){
  
  # given a vector `lx` containing the logs of the elements of a vector x,
  # stably compute log(mean(x))
  
  # sum the logs
  
  l_sum <- lx[1]
  
  for(i in 2:length(lx)){
    
    l_sum <- logapb(l_sum, lx[i])
    
  }
  
  # divide by the number
  
  ans <- l_sum - log(length(lx))
  
  return(ans)
  
}

logint_mvprobit <- function(lower, 
                            upper, 
                            R, 
                            niter){
  
  # log-integral of a zero-centred multivariate probit with correlation matrix `R`, 
  # between the vectors `lower` and `upper`. Estimation by a Monte-Carlo-like 
  # simulation using the method of Botev (2015), with `niter` draws. Inspired by
  # (and in plces borrowing heavily from) the implementation in the
  # TruncatedNormal package, but optimised for this problem. All bugs my own
  # etc.
  
  require(TruncatedNormal)
  
  # dimension of MVN (number of species)
  
  m <- length(lower)
  
  # truncated Cholesky decomposition
  
  chol_result <- TruncatedNormal::cholperm(R, lower, upper)
  
  L <- chol_result$L
  
  D <- diag(L)
  
  if(any(D < 10 ^ -10)){
    
    warning("Method may fail as covariance matrix is singular!")
    
  }
  
  # scaling
  
  l <- chol_result$l / D
  
  u <- chol_result$u / D
  
  L <- L / D - diag(m)
  
  # magic
  
  xmu <- TruncatedNormal::nleq(l, u, L)
  
  mu <- xmu[m:(2 * m - 2)]
  
  # loop through dimensions adding to the log integral  
  
  mu[m] <- 0
  
  Z <- matrix(0, m, niter)
  
  lp <- 0
  
  for(k in 1:(m - 1)){
    
    col <- t(L[k, 1:k]) %*% Z[1:k, ]
    
    tl <- l[k] - mu[k] - col
    
    tu <- u[k] - mu[k] - col
    
    Z[k, ] <- mu[k] + TruncatedNormal::trandn(tl, tu)
    
    lp <- lp + TruncatedNormal::lnNpr(tl, tu) + 0.5 * mu[k]^2 - mu[k] * Z[k, ]
    
  }
  
  col <- L[m, ] %*% Z
  
  tl <- l[m] - col
  
  tu <- u[m] - col
  
  lp <- lp + TruncatedNormal::lnNpr(tl, tu)
  
  # get log mean across all iterations and return
  
  lp <- logmean(lp)
  
  lp
  
}

lik_probit <- function(obs = NULL, 
                       mu = NULL, 
                       R = NULL, 
                       log.p = FALSE, 
                       niter = 1000){
  
  # given a binary matrix `obs` of 1s and 0s (rows being sites and columns being
  # species), a matching matrix `mu`, each row giving the mean of the
  # multivariate normal across species for that site, and a species-species
  # correlation matrix `R`, use `niter` Monte Carlo samples to estimate the
  # probability of observing that vector of observations, under the multivariate
  # probit model. Returns a vector of probabilities (by default on the log
  # scale), corresponding to each site.
  
  # check and get dimensions
  
  stopifnot(all(dim(obs) == dim(mu)))
  
  m <- ncol(obs)
  
  n <- nrow(obs)
  
  # set up truncation matrices
  
  lower <- matrix(-Inf, nrow = n, ncol = m)
  
  upper <- matrix(Inf, nrow = n, ncol = m)
  
  # The subsequent code assumes a zero-mean MVN, so need to adjust truncation to
  # account for mean, hence minus signs on mu
  
  present <- obs == 1
  
  lower[present] <- -mu[present]
  
  upper[!present] <- -mu[!present]
  
  # get probabilities for each row
  
  ans <- rep(NA, n)
  
  for(i in 1:n){
    
    ans[i] <- logint_mvprobit(lower[i, ],
                              upper[i, ],
                              R,
                              niter)
    
  }
  
  # on probability scale if requested
  
  if(!log.p){
    
    ans <- exp(ans)
    
  }
  
  return(ans)
  
}

### My new wrapper around Nick's functions to fit the overall workflow

independent_log_likelihood <- function(Beta = NULL,
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
  
  ### For each slice of array
  
  for(s in seq_len(n_iter)){
    
    #### Prediction for species assemblage at site i using values from slice a
    
    likelihood <- lik_probit(obs = y,
                             mu = mean_values[ , , s],
                             R = diag(n_species),
                             log.p = FALSE,
                             niter = 1000)
    
    #### Fill predictions array with value
    
    log_lik[ , s] <- log(likelihood)
    
    
  } 
  
  # ## Calculate single likelihood value for whole model
  # 
  # ### Take the product of site-level likelihoods within a single sample
  # 
  # sum_log_lik <- colSums(log_lik,
  #                        na.rm = TRUE)
  # 
  # ### Take the mean likelihood across samples
  # 
  # mean_log_lik <- mean(sum_log_lik,
  #                      na.rm = TRUE)
  # 
  # return(mean_log_lik)
  # 
  
  return(log_lik)
  
}  


joint_log_likelihood <- function(Beta = NULL,
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
  
  ### For each slice of array
  
  for(s in seq_len(n_iter)){
    
    #### Prediction for species assemblage at site i using values from slice a
    
    likelihood <- lik_probit(obs = y,
                             mu = mean_values[ , , s],
                             R = R[ , , s],
                             log.p = FALSE,
                             niter = 1000)
    
    #### Fill predictions array with value
    
    log_lik[ , s] <- log(likelihood)
    
    
  } 
  
  # ## Calculate single likelihood value for whole model
  # 
  # ### Take the product of site-level likelihoods within a single sample
  # 
  # sum_log_lik <- colSums(log_lik,
  #                        na.rm = TRUE)
  # 
  # ### Take the mean likelihood across samples
  # 
  # mean_log_lik <- mean(sum_log_lik,
  #                      na.rm = TRUE)
  # 
  # return(mean_log_lik)
  # 
  
  return(log_lik)
  
}  
