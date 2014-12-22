.get.check.points <- function(anno,PARAMETERS){
  
  # number of points
  no_points=1+ceiling(anno$exome_length/PARAMETERS$SLIDING_STEP)
  
  # generate points
  points=seq(from = 1, to = anno$exome_length, length.out = no_points)
  
  # round
  points=as.integer(round(points))
  
  # save result
  return(points)
}