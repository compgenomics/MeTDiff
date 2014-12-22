# get bam info
.get.bam.chrs<- function(file){
  
  # index
  if (! file.exists(paste(file,'.bai',sep=""))) {
    print(paste("Stage: index bam file", file))
    date()
    indexBam(file)
  }
  
  param <- ScanBamParam(what="rname")
  bam <- scanBam(file, param=param)
  chr=levels(bam[[1]]$rname)
  
  # save result
  return(chr)
  
}
