.help.factorial <- function(count){
  #compute the log(count!)
  cm = max(count)
  if (is.null(ncol(count))){
    D <- 1
  }else{
    D <- ncol(count)
  }
  
  if(cm > 50000){
    dnorm <- as.matrix(lgamma(data.matrix(count+1)))
  }
  else{
    tmp  <- cumsum(rbind(0,log(as.matrix(1:max(count)))))
    dnorm <- matrix(tmp[data.matrix(count+1)],ncol=D)
  }
}
