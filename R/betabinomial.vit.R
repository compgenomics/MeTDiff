.betabinomial.vit <- function(x,y,trans,alpha,beta)
  # use viterbi algorithm to deconvolute the status
{
  if (any(is.na(alpha))|any(is.na(beta))){
    alpha <- c(5,1)
    beta <- c(1,5)
  }
  N <- 2 #two states in our model
  T <- nrow(x)
  dx <- .help.factorial(x)
  dy <- .help.factorial(y)
  n  <- x + y
  dn  <- .help.factorial(n)
  dens <- .help.postprob(dx,dy,dn,x,y,alpha,beta)
  logdens <- log(dens + .Machine$double.xmin) # to avoid inf
  
  #dynamic programming recursion
  trans <- log(trans + .Machine$double.xmin)
  plogl <- matrix(0,T,N)
  bcktr <- matrix(0,T-1,N)
  
  plogl[1,] <- logdens[1,] - log(N)
  for (t in 2:T){
    tmp <- t(plogl[t-1,,drop=F]) %*% matrix(1,1,N) + trans
    plogl[t,] <- apply(tmp,2,max)
    bcktr[t-1,] <- max.col(t(tmp))
    plogl[t,] <- plogl[t,] + logdens[t,]
    
  }
  
  class <- matrix(0,nrow=T)
  class[T] <- which.max(plogl[T,])
  
  for (t in (T-1):1){
    class[t] <- bcktr[t,class[t+1]]
  }
  list(class=class, logl=max(plogl[T,]))
}

