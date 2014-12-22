
# get blocks
.get.block.from.rna.peak <- function(rna_peak,anno){
  
  # get peak
  rna_start=which(anno$DNA2RNA==rna_peak[1])
  rna_stop=which(anno$DNA2RNA==rna_peak[2])
  
  # extract peak region
  rna_peak=anno$DNA2RNA[rna_start:rna_stop]
  
  # get information
  left=anno$left+rna_start-1
  right=anno$left+rna_stop-1
  
  # get block information
  dna_length=length(rna_peak)
  rna_peak=(rna_peak>0)
  block_start=which((rna_peak-c(0,rna_peak[1:(dna_length-1)]))>0)
  block_end=which((rna_peak-c(rna_peak[2:dna_length],0))>0)
  
  # get block info
  blockSizes=character(0)
  blockStarts=character(0)
  block_length=block_end-block_start+1
  for (i in 1:length(block_start)) {
    blockSizes=paste(blockSizes,as.character(block_length[i]),',',sep="")
    if (i < length(block_start)) {blockStarts=paste(blockStarts,as.character(block_start[i]-1),',',sep="")}
    if (i == length(block_start)) {blockStarts=paste(blockStarts,as.character(block_start[i]-1),sep="")}
  }
  
  
  # get result
  block=list(
  chromStart=left-1,
  chromEnd=right,
  thickStart=left-1,
    thickEnd=right,
    blockCount=length(block_start),
    blockSizes=blockSizes,
    blockStarts=blockStarts)
  
  return(block)
}