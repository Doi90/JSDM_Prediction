###############################################################
###############################################################
###                                                         ###
###             NEW ptmvtnorm() FUNCTION                    ###
###                                                         ###
###  This script rebuilds the ptmvtnorm() function from the ###
### tmvtnorm package so that it exposes the algorithm       ###
### argument in the underlying pmvnorm() function.          ###
###                                                         ###
###############################################################
###############################################################

## Define the two tmvtnorm() functions that the function
## requires but for some reason aren't loaded with the
## package. Copied straight from the GitHub page

checkSymmetricPositiveDefinite <- function(x, name="sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  
  if (any(diag(x) <= 0)) {
    stop(sprintf("%s all diagonal elements must be positive", name))
  }
  
  if (det(x) <= 0) {
    stop(sprintf("%s must be positive definite", name))
  }
}

# Uses partly checks as in mvtnorm:::checkmvArgs!
checkTmvArgs <- function(mean, sigma, lower, upper)
{
  if (is.null(lower) || any(is.na(lower))) 
    stop(sQuote("lower"), " not specified or contains NA")
  if (is.null(upper) || any(is.na(upper))) 
    stop(sQuote("upper"), " not specified or contains NA")
  if (!is.numeric(mean) || !is.vector(mean)) 
    stop(sQuote("mean"), " is not a numeric vector")
  if (is.null(sigma) || any(is.na(sigma))) 
    stop(sQuote("sigma"), " not specified or contains NA")
  
  if (!is.matrix(sigma)) {
    sigma <- as.matrix(sigma)
  }
  
  if (NCOL(lower) != NCOL(upper)) {
    stop("lower and upper have non-conforming size")
  }
  
  checkSymmetricPositiveDefinite(sigma)
  
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  if (length(lower) != length(mean) || length(upper) != length(mean)) {
    stop("mean, lower and upper must have the same length")
  }
  
  if (any(lower>=upper)) {
    stop("lower must be smaller than or equal to upper (lower<=upper)")
  }
  
  # checked arguments
  cargs <- list(mean=mean, sigma=sigma, lower=lower, upper=upper)
  return(cargs)
}

## Define the ptmvtnorm.new() function

ptmvtnorm.new <- function (lowerx,
                           upperx,
                           mean = rep(0, length(lowerx)),
                           sigma,
                           lower = rep(-Inf, length = length(mean)),
                           upper = rep(Inf, length = length(mean)),
                           maxpts = 25000,
                           abseps = 0.001,
                           releps = 0,
                           algorithm = "GenzBretz"){
  
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  
  mean <- cargs$mean
  
  sigma <- cargs$sigma
  
  lower <- cargs$lower
  
  upper <- cargs$upper
  
  if(is.null(lowerx) || any(is.na(lowerx))){ 
    stop(sQuote("lowerx"), " not specified or contains NA")
  }
  
  if(is.null(upperx) || any(is.na(upperx))){ 
    stop(sQuote("upperx"), " not specified or contains NA")
  }
  
  if(!is.numeric(lowerx) || !is.vector(lowerx)){ 
    stop(sQuote("lowerx"), " is not a numeric vector")
  }
  
  if(!is.numeric(upperx) || !is.vector(upperx)){ 
    stop(sQuote("upperx"), " is not a numeric vector")
  }
  
  if(length(lowerx) != length(lower) || length(lower) != length(upperx)){ 
    stop("lowerx an upperx must have the same length as lower and upper!")
  }
  
  if(any(lowerx >= upperx)){ 
    stop("lowerx must be smaller than or equal to upperx (lowerx<=upperx)")
  }
  
  f <- pmvnorm(lower = pmax(lowerx, lower),
               upper = pmin(upperx, upper),
               mean = mean,
               sigma = sigma,
               maxpts = maxpts, 
               abseps = abseps,
               releps = releps,
               algorithm = algorithm) / pmvnorm(lower = lower,
                                                upper = upper, 
                                                mean = mean, 
                                                sigma = sigma, 
                                                maxpts = maxpts, 
                                                abseps = abseps, 
                                                releps = releps, 
                                                algorithm = algorithm)
  
  return(f)

}