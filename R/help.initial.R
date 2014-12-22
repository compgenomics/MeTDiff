.help.initial <- function(xx,yy,postprob){
  
  N <- ncol(postprob)
  Tm <- dim(xx)
  if (is.null(Tm)){
    m <- 1
  }else{
    m <- Tm[2]
  }
  wght <- apply(postprob,2,sum)
  alpha <- matrix(0,m,N)
  beta <- matrix(0,m,N)
  for (ii in 1:m){
    x <- as.matrix(xx[,ii])
    y <- as.matrix(yy[,ii])
    ep <- apply(postprob * cbind(x/(x+y),x/(x+y)),2,sum) /  wght
    ep_2 <- apply(postprob * cbind(x/(x+y),x/(x+y))**2,2,sum) / wght
    s <- ep_2 / ep
    alpha[ii,] <- (s-1) / (1-s/ep)
    beta[ii,] <- alpha[ii,] / ep - alpha[ii,]
  }
  alpha0 <- colMeans(alpha)
  beta0 <- colMeans(beta)
  list(alpha=alpha0,beta=beta0)
}
