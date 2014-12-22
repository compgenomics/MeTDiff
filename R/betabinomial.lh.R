.betabinomial.lh <- function(x,y,Nit=40,Npara=1e-9){
#   x <- as.matrix(x[peak,])
#   y <- as.matrix(y[peak,])
  N <- 2 # number of states
  J <- matrix(0,N,1)
  H <- matrix(0,N,N)
  T <- nrow(x)
  n <- x + y
  m <- ncol(x)
  rr <- x/n
  # parameters initialization
#   exx <- colMeans(rr)
#   vxx <- apply(rr,2,sd)^2
#   alpha <- mean( exx*(exx*(1-exx)/vxx-1) )
#   beta <- mean( (1-exx)*(exx*(1-exx)/vxx-1) )

  # use another method to initialize
  p1_e <- exp(sum( log(rr) )/(T*m))
  p2_e <- exp(sum( log(1-rr)/(T*m) ))
  alpha <- 1/2 *(1-p2_e)/(1-p1_e-p2_e ) # to avoid 0
  beta <-  1/2 *(1-p1_e)/(1-p1_e-p2_e )
  c = rbind(alpha,beta)

# add break condition to avoid alpha is na   
 if ( !any(is.finite(beta)) | !is.finite(alpha) | any(beta <= 0) | any(alpha<= 0) ){
   return(list(logl=rnorm(1)*10000,alpha=c(1,1),beta=c(1,1))) 
 }


  
  for (nit in 1:Nit){
    J[1] <- T*digamma(sum(c))*m - sum( .help.digamma(as.matrix(n),sum(c)) ) + sum( .help.digamma(as.matrix(x),c[1]) ) - T*digamma(c[1])*m
    J[2] <- T*digamma(sum(c))*m - sum( .help.digamma(as.matrix(n),sum(c)) ) + sum( .help.digamma(as.matrix(y),c[2]) ) - T*digamma(c[2])*m

    H[1,1] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c))) + sum(.help.trigamma(as.matrix(x),c[1])) - T*trigamma(c[1])*m
    H[2,2] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c))) + sum(.help.trigamma(as.matrix(y),c[2]))  - T*trigamma(c[2])*m    
    H[1,2] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c)))
    H[2,1] <- H[1,2]
    eigvalue <- eigen(H)$values
    
    if ( (any(beta < Npara)) | (any(alpha < Npara)) 
        | abs(eigvalue[1]/eigvalue[2]) > 1e12 | abs(eigvalue[1]/eigvalue[2]) < 1e-12){   break  }
    
#     tmp_step <- -solve(H,tol=1e-20) %*% J
    tmp_step <- -solve(H, J) # using newton smoothing
    tmp <- c + tmp_step
    while(any(tmp <= 0)){
#       warning(sprintf("Could not update the Newton step ...\n"))
      tmp_step <- tmp_step / 20
      tmp <- c + tmp_step
    }
    c <- tmp
    
  }
#   caculate the likelihood
  alpha <- c[1]
  beta <- c[2]
  dnx <- .help.factorial(x)
  dny <- .help.factorial(y)
  dn <- .help.factorial(n)
  prob <- .help.postprob(dnx,dny,dn,x,y,alpha,beta)
  return(list(logl=sum(log(prob)),alpha=alpha,beta=beta))

}