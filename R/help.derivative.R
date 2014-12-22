.help.digamma <- function(xx,alpha){
  Tm <- dim(xx)
  TT <- Tm[1]
  m <- Tm[2]
  
  res <- matrix(0,TT,1)
  
  for (ii in 1:m){
    res <- res + digamma(xx[,ii] + alpha)
  }
  return(res)
}

.help.trigamma <- function(xx,alpha){
  Tm <- dim(xx)
  if (is.null(Tm)){
    TT <- length(xx)
    m <- 1
  }else{
    TT <- Tm[1]
    m <- Tm[2]
  }
  res <- matrix(0,TT,1)
  
  for (ii in 1:m){
    res <- res + trigamma(xx[,ii] + alpha)
  }
  return(res)
  
}