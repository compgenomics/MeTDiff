# get bam info
.get.bam.read.length<- function(file){
  
  # index
  if (! file.exists(paste(file,'.bai',sep=""))) {
    print(paste("Stage: index bam file", file))
    date()
    indexBam(file)
  }
  
  param <- ScanBamParam(what="qwidth")
  bam <- scanBam(file, param=param)
  read_length=round(median(bam[[1]]$qwidth))
  
  # save result
  return(read_length)
  
}
