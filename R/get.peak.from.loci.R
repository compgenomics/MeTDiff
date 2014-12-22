.get.peak.from.loci <- function(READS_COUNT,ID,PARAMETERS){
  # get no_line
  no_line=length(ID)
  
  # start id
  start_id = which((ID[2:no_line]-ID[1:no_line-1]==1) |   
    ((READS_COUNT$batch_id[2:no_line]!=READS_COUNT$batch_id[1:no_line-1]) & (ID[2:no_line] ==TRUE)) )
  start_id=start_id + 1
  if ( ID[1]==TRUE ) { start_id = c(1,start_id) }
  
  # end id 
  end_id = which((ID[1:no_line-1]-ID[2:no_line]==1) |   
                     ((READS_COUNT$batch_id[1:no_line-1]!=READS_COUNT$batch_id[2:no_line]) & (ID[1:no_line-1] ==TRUE)) )
  if (ID[no_line]==TRUE) {end_id=c(end_id,no_line)}  
  
  # label peaks
  PEAK_LENGTH=READS_COUNT$check_points[end_id]-READS_COUNT$check_points[start_id]
  good_peak_id=which(PEAK_LENGTH > PARAMETERS$MINIMAL_PEAK_LENGTH)
  start_id=start_id[good_peak_id]
  end_id=end_id[good_peak_id]
  
  
  # result
  peak=cbind(start_id,end_id)

return(peak)
}