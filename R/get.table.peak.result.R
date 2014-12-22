
# get table peak result
.get.table.peak.result <- function(PEAK,ANNOTATION,READS_COUNT,SAMPLE_ID,PARAMETERS,
                                  ANNOTATION_BATCH_ID,LOCI2PEAK) {  
  
  
  
  # initializa the peak reporting
  no_peak=length(LOCI2PEAK[,1])
  peak_report=data.frame()
  if (no_peak ==0) {return(peak_report)} else {
  # get peak
  for (i in 1:no_peak) {
    peak_row_id=LOCI2PEAK[i,]
    rna_peak=READS_COUNT$check_points[peak_row_id]
    
    # batch id
    batch_id=unique(READS_COUNT$batch_id[peak_row_id])
    lg.p=min(PEAK$PW$log.p[peak_row_id[1]:peak_row_id[2]])/log(10)
    lg.fdr=min(PEAK$PW$log.fdr[peak_row_id[1]:peak_row_id[2]])/log(10)
    fold_enrchment=exp(max(PEAK$PW$log.fc[peak_row_id[1]:peak_row_id[2]]))
    
    # get sig digits
    lg.p=signif(lg.p, digits = 3)
    lg.fdr=signif(lg.fdr, digits = 3)
    fold_enrchment=signif(fold_enrchment, digits = 3)
    
    # get annotation
    anno = .get.gene.anno(batch_id,ANNOTATION,ANNOTATION_BATCH_ID)
    
    # get blocks
    block=.get.block.from.rna.peak(rna_peak,anno)
    
    # save result
    xls=data.frame(chr=anno$chr,
             chromStart=block$chromStart,
               chromEnd=block$chromEnd,
               name=anno$gene,
                score=signif(10^lg.p,digits=2),
             strand=anno$strand,
             thickStart=block$thickStart,
               thickEnd=block$thickEnd,
               itemRgb=0,
               blockCount=block$blockCount,
               blockSizes=block$blockSizes,
               blockStarts=block$blockStarts,
               lg.p=lg.p,
             lg.fdr=lg.fdr,
             fold_enrchment=fold_enrchment)
  
    # append result
    peak_report=rbind(peak_report,xls)
  }
  return(peak_report)}
}