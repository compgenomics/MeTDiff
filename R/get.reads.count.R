

# get reads count
.get.reads.count <- function(ibatch,PARAMETERS,ANNOTATION,ANNOTATION_BATCH_ID,BAM_CHRS,no_bam_files,bam){
  
  # get annotation
  anno = .get.gene.anno(ibatch,ANNOTATION,ANNOTATION_BATCH_ID)
  check_points=.get.check.points(anno,PARAMETERS)
  no_check_points=length(check_points)
  
  # get counts
  reads_count=matrix(as.integer(0),nrow=no_check_points,ncol=no_bam_files)
  if (sum(BAM_CHRS==anno$chr)>0) {
    for (ibam in 1:no_bam_files){
      reads_count[,ibam]=.get.check.points.reads.count(ibam,anno,bam,check_points,PARAMETERS)
    }
    # compile result
    batch_id=rep(ibatch,no_check_points)
    report=cbind(reads_count,batch_id,check_points)
  } else {
    batch_id=rep(ibatch,no_check_points)
    report=cbind(reads_count,batch_id,check_points)
    report=report[character(0),]
  }
  return(report)
}


