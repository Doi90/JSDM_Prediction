proportionNA <- function(a){
  
  count <- sum(is.na(a))
  
  length <- length(a)
  
  proportion <- count / length
  
  return(proportion)
  
}

