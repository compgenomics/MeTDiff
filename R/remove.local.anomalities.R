.remove.local.anomalities <- function(pos_table){
  
  # parameters
  max_background_fold_increase = 4
  background_window= 200
  
  
  # get position
  pos_mapped = as.numeric(names(pos_table))
  no_pos_mapped=length(pos_mapped)
  
  # prepare new table
  new_table=pos_table
  
  # filter
  for (i in 1:no_pos_mapped) {
    ID=which(abs(pos_mapped[i]-pos_mapped) < (background_window/2))
    else_count = sum(pos_table[ID])-sum(pos_table[i])
    max_background = round(else_count*max_background_fold_increase/background_window)
    new_table[i] = max(1,min(pos_table[i],max_background))
  }
  
  return(new_table)
}
