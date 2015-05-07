
.betabinomial.hmm  <- function(x,y,alpha=c(5,1),beta=c(1,5),trans=matrix(0.5,2,2),Nit=39,Nerr=1e-3,Npara=1e-6)
{
  
  n <- x+y
  dnx <- .help.factorial(x)
  dny <- .help.factorial(y)
  dn <- .help.factorial(n)
  T <- nrow(n)
  m <- ncol(n)
  #initialization
  logl <- numeric(Nit)
  logl[c(1,2)] <- -1e5
  N <- 2 # number of states
  dens <- matrix(0,T,N)
  J <- matrix(0,N,1)
  H <- matrix(0,N,N)
  forwrd <- matrix(0,T,N)
  bckwrd <- matrix(0,T,N)
  scale <- rbind(1,matrix(0,T-1))
  postprob <- cbind(matrix(0,T,1),matrix(1,T,1))
  
  alpha0 <- alpha
  beta0 <- beta
  trans0 <- trans
  postprob0 <- postprob #before the first step caculation, initial postprob0 and nit0
  nit0 <- 3 #initial t
  
  for (nit in 3:Nit){
    # 1: E-step to compute density
    dens <- .help.postprob(dnx,dny,dn,x,y,alpha,beta)
    #     dens <- apply(dens,2,function(x) {x[x==0] <- .Machine$double.xmin; return (x)})
    # 2: E-step forward recursion
    forwrd[1,] <- dens[1,]/N
    for (t in 2:T){
      forwrd[t,] <- (forwrd[t-1,] %*% trans) * dens[t,]
      #system scaling
      scale[t] <- sum(forwrd[t,])
      forwrd[t,] <- forwrd[t,] / scale[t]
    }
    #compute log-likelyhood
    logl[nit] <- log(sum(forwrd[T,])) + sum(log(scale))
    #         message(sprintf('Iteration %d:\t%.3f', (nit-1), logl[nit]))
    
    # break condition
    if (logl[nit] > -20 | T < 10 | is.na(logl[nit]) 
        | (logl[nit] - logl[nit-1]) < -50 ) {break}
    
    #3: E-step backward recursion
    bckwrd[T,] <- matrix(1,1,N)
    for (t in (T-1):1){
      bckwrd[t,] <- (bckwrd[t+1,] * dens[t+1,]) %*% t(trans)
      bckwrd[t,] <- bckwrd[t,] / scale[t]
    }
    
    #4: M-step reestimate transition matrix
    trans <- trans * (t(forwrd[1:(T-1),]) %*% (dens[2:T,] * bckwrd[2:T,] ))
    trans <- trans / (apply(trans,1,sum)) %*% matrix(1,1,N) 
    #5: M-step reestimate parameters
    postprob <- forwrd * bckwrd
    postprob <- postprob / (apply(postprob,1,sum)) %*% matrix(1,1,N)
    
    wght <- apply(postprob,2,sum)
    par <- .help.initial(x,y,postprob)
    alpha <- par$alpha
    beta <- par$beta
    
    
    c <- rbind(alpha[1],beta[1])
    if(  any(is.na(beta)) | any(is.na(alpha)) | any(beta <= 0) | any(alpha<= 0)  ){break}
    #6: M-step reestimate parameters using newton-raphson method for the first component
    J[1] <- wght[1]*digamma(sum(c))*m - t(postprob[,1]) %*% (.help.digamma(as.matrix(n),sum(c))) + t(postprob[,1]) %*% (.help.digamma(as.matrix(x),c[1])) - wght[1]*digamma(c[1])*m
    J[2] <- wght[1]*digamma(sum(c))*m - t(postprob[,1]) %*% (.help.digamma(as.matrix(n),sum(c))) + t(postprob[,1]) %*% (.help.digamma(as.matrix(y),c[2])) - wght[1]*digamma(c[2])*m
    
    #     J[1] <- wght[1]*digamma(sum(c)) - t(postprob[,1]) %*% (digamma(n+sum(c))) + t(postprob[,1]) %*% (digamma(x+c[1])) - wght[1]*digamma(c[1])
    #     J[2] <- wght[1]*digamma(sum(c)) - t(postprob[,1]) %*% (digamma(n+sum(c))) + t(postprob[,1]) %*% (digamma(y+c[2])) - wght[1]*digamma(c[2])  
    #     J <- matrix(1,N) %*% (wght[1]*digamma(sum(c)) - t(postprob[,1]) %*% (digamma(n+sum(c)))) +
    #       t(cbind(digamma(x+c[1]),digamma(y+c[2]))) %*% postprob[,1] - wght[1] * rbind(digamma(c[1]),digamma(c[2]))
    
    H[1,1] <- wght[1]*trigamma(sum(c))*m - t(postprob[,1]) %*% (.help.trigamma(as.matrix(n),sum(c))) + t(postprob[,1]) %*% (.help.trigamma(as.matrix(x),c[1])) - wght[1]*trigamma(c[1])*m
    H[2,2] <- wght[1]*trigamma(sum(c))*m - t(postprob[,1]) %*% (.help.trigamma(as.matrix(n),sum(c))) + t(postprob[,1]) %*% (.help.trigamma(as.matrix(y),c[2]))  - wght[1]*trigamma(c[2])*m    
    H[1,2] <- wght[1]*trigamma(sum(c))*m - t(postprob[,1]) %*% (.help.trigamma(as.matrix(n),sum(c)))
    H[2,1] <- H[1,2]
    eigvalue <- eigen(H)$values
    if ( (any(beta < Npara)) | (any(alpha < Npara))
         | (abs(eigvalue[1]/eigvalue[2]) > 1e12) | (abs(eigvalue[1]/eigvalue[2]) < 1e-12)
         | any(eigvalue==0) ) {   break  }
    
    #   tmp_step <- -solve(H,tol=1e-20) %*% J 
    tmp_step <- -solve(H,J) #using gaussian smoothing
    tmp <- c + tmp_step
    while(any(tmp <= 0)){
      warning(sprintf("Could not update the Newton step ...\n"))
      tmp_step <- tmp_step / 20
      tmp <- c + tmp_step
    }
    alpha[1] <- tmp[1]
    beta[1] <- tmp[2]
    #7: M-step reestimate parameters for the second component
    c <- rbind(alpha[2],beta[2])
    J[1] <- wght[2]*digamma(sum(c))*m - t(postprob[,2]) %*% (.help.digamma(as.matrix(n),sum(c))) + t(postprob[,2]) %*% (.help.digamma(as.matrix(x),c[1])) - wght[2]*digamma(c[1])*m
    J[2] <- wght[2]*digamma(sum(c))*m - t(postprob[,2]) %*% (.help.digamma(as.matrix(n),sum(c))) + t(postprob[,2]) %*% (.help.digamma(as.matrix(y),c[2])) - wght[2]*digamma(c[2])*m
    H[1,1] <- wght[2]*trigamma(sum(c))*m - t(postprob[,2]) %*% (.help.trigamma(as.matrix(n),sum(c))) + t(postprob[,2]) %*% (.help.trigamma(as.matrix(x),c[1])) - wght[2]*trigamma(c[1])*m
    H[2,2] <- wght[2]*trigamma(sum(c))*m - t(postprob[,2]) %*% (.help.trigamma(as.matrix(n),sum(c))) + t(postprob[,2]) %*% (.help.trigamma(as.matrix(y),c[2]))  - wght[2]*trigamma(c[2])*m
    H[1,2] <- wght[2]*trigamma(sum(c))*m - t(postprob[,2]) %*% (.help.trigamma(as.matrix(n),sum(c)))
    H[2,1] <- H[1,2]
    # if the H is nearly sigular then break
    eigvalue <- eigen(H)$values 
    
    if ( (any(beta < Npara)) | (any(alpha < Npara))
         | (abs(eigvalue[1]/eigvalue[2]) > 1e12) | (abs(eigvalue[1]/eigvalue[2]) < 1e-12)
         | any(eigvalue==0) ) {   break  }
    
    #     tmp_step <- -solve(H,tol=1e-20) %*% J
    tmp_step <- -solve(H,J)
    tmp <- c + tmp_step
    while(any(tmp <= 0)){
      warning(sprintf("Could not update the Newton step ...\n"))
      tmp_step <- tmp_step / 20
      tmp <- c + tmp_step
    }
    alpha[2] <- tmp[1]
    beta[2] <- tmp[2]
    
    # parameters backup
    alpha0 <- alpha
    beta0 <- beta
    postprob0 <- postprob
    trans0 <- trans
    nit0 <- nit
    #8: Break if conditions are satisfied
    if ( ( (abs(logl[nit]-logl[nit-1]) + abs(logl[nit]-logl[nit-2]) ) <= 5*Nerr ) 
         | (any(beta < Npara)) | (any(alpha < Npara))  ){
      
      break; #if likelihood decrease less than nerr and alpha and beta are approximately to 0
    }
    
  }
  
  list(alpha=alpha0,beta=beta0,trans=trans0,logl=logl[3:nit0],postprob=postprob0)
  
}