.help.postprob <- function(dxx,dyy,dnn,xx,yy,alpha,beta){
  
  N <- length(alpha)
  Tm <- dim(xx)
  if (is.null(Tm)){
    TT <- length(xx)
    m <- 1
  }else{
    TT <- Tm[1]
    m <- Tm[2]
  }
  res <- matrix(0,TT,N)
  
  for (ii in 1:m){
    dnx <- as.matrix(dxx[,ii])
    dny <- as.matrix(dyy[,ii])
    dn <- as.matrix(dnn[,ii])
    x <- as.matrix(xx[,ii])
    y <- as.matrix(yy[,ii])
    res <- res + (dn-dnx-dny) %*% matrix(1,1,N) + lgamma(x %*% matrix(1,1,N) + matrix(1,TT) %*% alpha) - 
      lgamma(matrix(1,TT) %*% (alpha+beta) + (x+y) %*% matrix(1,1,N)) +
      lgamma(y %*% matrix(1,1,N) + matrix(1,TT) %*% beta) + lgamma(matrix(1,TT) %*% (alpha+beta)) - 
      lgamma(matrix(1,TT) %*% alpha) - lgamma(matrix(1,TT) %*% beta)
  }
  res <- exp(res)
}