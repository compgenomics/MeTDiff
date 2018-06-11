
# get check point reads count
.get.check.points.reads.count<- function(ibam,anno,bam,check_points,PARAMETERS){
  
  # prepare bam parameters
  which <- IRangesList(chr_unclear=IRanges(anno$left, anno$right))
  names(which)=anno$chr
  what <- c("strand", "pos","mapq")
  param <- ScanBamParam(which=which, what=what)
  
  # read bam file
  ba <- scanBam(bam[ibam], param=param)
  pos=ba[[1]]$pos-anno$left+1
  strand=ba[[1]]$strand
  mapq=ba[[1]]$mapq
  
  # fileter mapq
  mapq[which(is.na(mapq))]=255
  ID=which(mapq>=PARAMETERS$MINIMAL_MAPQ)
  pos=pos[ID]
  strand=strand[ID]
  
  # filter nagative pos
  ID=which(pos>0)
  pos=pos[ID]
  strand=strand[ID]
  
  # convert pos into rna
  rna_pos=anno$DNA2RNA[pos]
  on_rna_id=which( ((rna_pos>0) + !is.na(rna_pos)) ==2)
  rna_pos=rna_pos[on_rna_id]
  strand=strand[on_rna_id]
  
  # divide into strand
  pos_ID=which(strand=="+")
  neg_ID=which(strand=="-")
  pos_pos=rna_pos[pos_ID]
  neg_pos=rna_pos[neg_ID]
  
  # shift
  pos_pos=pos_pos+round(PARAMETERS$FRAGMENT_LENGTH/2);
  neg_pos=neg_pos+PARAMETERS$READ_LENGTH-round(PARAMETERS$FRAGMENT_LENGTH/2)
  
  # merge the two
  pos=c(pos_pos,neg_pos)
  pos=pos[pos > 0]
  pos=pos[pos < anno$exome_length]
  
  # get direct count
  check_points_count=check_points
  pos_table=table(pos)
  
  # smooth
  if (PARAMETERS$REMOVE_LOCAL_TAG_ANOMALITIES==TRUE) {pos_table=.remove.local.anomalities(pos_table)}
  
  # get count
  pos_mapped = as.numeric(names(pos_table))
  for (i in 1:length(check_points)) { 
      ID=which(abs(check_points[i]-pos_mapped)*2 < PARAMETERS$WINDOW_WIDTH)
      check_points_count[i]=sum(pos_table[ID])
  }
  
  # return result
  return(check_points_count)
}

