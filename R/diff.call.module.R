
# merge and compare two conditions
.diff.call.module <- function(PEAK,READS_COUNT,SAMPLE_ID,PARAMETERS){
  
  no_peak=length(PEAK$loci2peak_merged[,1])
  pvalues <- rep(1,no_peak)
  log.fc <- rep(0,no_peak)
  for (ipeak in 1:no_peak) {
    print(ipeak)
    temp=PEAK$loci2peak_merged[ipeak,]
    x <- as.matrix(READS_COUNT[temp[1]:temp[2],SAMPLE_ID$untreated_ip])
    y <- as.matrix(READS_COUNT[temp[1]:temp[2],SAMPLE_ID$untreated_input])
    xx <- as.matrix(READS_COUNT[temp[1]:temp[2],SAMPLE_ID$treated_ip])
    yy <- as.matrix(READS_COUNT[temp[1]:temp[2],SAMPLE_ID$treated_input])
    xxx = cbind(x,xx)
    yyy = cbind(y,yy)
    #BBtest
    logl1 <- .betabinomial.lh(x,y+1)
    logl2 <- .betabinomial.lh(xx,yy+1)
    logl3 <- .betabinomial.lh(xxx,yyy+1)
    tst <- (logl1$logl+logl2$logl-logl3$logl)*2
    pvalues[ipeak] <- 1 - pchisq(tst,2)
    log.fc[ipeak] <- log( (sum(xx)+1)/(1+sum(yy)) * (1+sum(y))/(1+sum(x)) ) 
    
  }
  
  log.p <- log(pvalues)
  log.fdr <- log(p.adjust(pvalues,method='fdr'))
  log.fdr <- pmax(log.fdr,-1000)
  log.p <- pmax(log.p,-1000)
  DIFF <- list(log.fdr=log.fdr,log.p=log.p,log.fc=log.fc)
  # result
  result =list()
  result$DIFF = DIFF
  #   result$peak_reads_count = list(untreated_ip=untreated_ip,
  #                                  untreated_input=untreated_input,
  #                                  treated_ip=treated_ip,
  #                                  treated_input=treated_input,
  #                                  untreated_ip_total=untreated_ip_total,
  #                                  untreated_input_total=untreated_input_total,
  #                                  treated_ip_total=treated_ip_total,
  #                                  treated_input_total=treated_input_total,
  #                                  sample_peak_reads_count=peak_reads_count,
  #                                  sample_total_reads_count = sample_total,
  #                                  sample_id=SAMPLE_ID)
  return(result)
  
}