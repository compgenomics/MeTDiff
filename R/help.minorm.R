.help.minorm <- function(READS_COUNT,SAMPLE_ID)
# using median to norm the samples depth
{
  
  TR <- apply(READS_COUNT[,1:length(SAMPLE_ID$sample_names)],2,sum)
  tmp <- prod(TR)^(1/length(SAMPLE_ID$sample_names))/TR
  for (i in 1:length(SAMPLE_ID$sample_names)){
    READS_COUNT[,i]  <- round(READS_COUNT[,i]*tmp[i])
  }
  return(READS_COUNT)
  
}
