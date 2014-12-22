


# get peak reads count
.get.peak.reads.count <- function(PEAK,READS_COUNT,SAMPLE_ID,PARAMETERS) {
  
  # get number of peaks
  no_peak=length(PEAK$loci2peak_merged[,1])
  
  # peak_reads_count
  peak_reads_count = READS_COUNT[1:no_peak,]
  no_sample=length(peak_reads_count[1,])-2
  
  # cut the unnecessary information
  peak_reads_count = READS_COUNT[1:no_peak,1:no_sample]
  
  # count
  for (ipeak in 1:no_peak) {
    temp=PEAK$loci2peak_merged[ipeak,]
    temp2=colSums(READS_COUNT[temp[1]:temp[2],1:no_sample])
    peak_reads_count[ipeak,1:no_sample]=temp2
  }
  
  # remove the overlapping window effects
  peak_reads_count = round (peak_reads_count * PARAMETERS$SLIDING_STEP / PARAMETERS$WINDOW_WIDTH);
  
return(peak_reads_count)
}