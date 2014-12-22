.peak.call.module <- function(READS_COUNT,SAMPLE_ID,PARAMETERS,ANNOTATION){
  
  
  PW = .metpeak.process(as.matrix(READS_COUNT[,SAMPLE_ID$ip]),as.matrix(READS_COUNT[,SAMPLE_ID$input]),READS_COUNT$batch_id) 
  if (PARAMETERS$PEAK_CUTOFF_TYPE =="FDR") {ID_stat= (PW$log.fdr < log(PARAMETERS$PEAK_CUTOFF_FDR))}
  if (PARAMETERS$PEAK_CUTOFF_TYPE =="PVALUE") {ID_stat= (PW$log.p < log(PARAMETERS$PEAK_CUTOFF_PVALUE))}
  ID_fc= (PW$log.fc > log(PARAMETERS$FOLD_ENRICHMENT))
  global.sig.check.point.id=ID_stat & ID_fc
  global.pw=PW
  
  # get sample consistency
  consistent.sid = global.sig.check.point.id
  
  # output peaks
  PEAK=list(PW=global.pw, Merged=global.sig.check.point.id, Consistent=consistent.sid)
  ID=PEAK$Merged
  PEAK$loci2peak_merged=.get.peak.from.loci(READS_COUNT,ID,PARAMETERS)
  ID=PEAK$Consistent
  PEAK$loci2peak_consistent=PEAK$loci2peak_merged
  
  return(PEAK)}
