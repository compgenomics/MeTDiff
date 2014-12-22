
.divide.anno.into.batches <- function(anno){
print("Divide transcriptome into chr-gene-batch sections ...")

# divide into chr-gene
chr_gene = paste(anno$chr,anno$gene,anno$strand)
unique_chr_gene = unique(chr_gene)
ID = match(chr_gene, unique_chr_gene)
gene_chr_ID=ID

# group
noGroups=ceiling(max(ID)/100);
group_bar=round(seq(from = 1, to = max(ID)+1,length.out=noGroups+1))

for (igroup in 1:noGroups) {
  # print(igroup)
  id_selected=which(((ID >= group_bar[igroup]) + (ID < group_bar[igroup+1]))==2)
  anno_small=anno[id_selected,]
  ID_small=ID[id_selected]-group_bar[igroup]+1
  
  # return result
  ID_small_new = .divide.anno.into.batches.small(anno_small,ID_small)
  ID_small_new[ID_small_new > max(ID_small)]= 
    ID_small_new[ID_small_new > max(ID_small)] - max(ID_small) + max(ID);
  ID_small_new[ID_small_new <= max(ID_small)]=
    ID_small_new[ID_small_new <= max(ID_small)] + group_bar[igroup] - 1
  ID[id_selected]=ID_small_new
}

return(as.integer(ID))

}

.divide.anno.into.batches.small <- function(anno,ID){
  # divide into chr-gene-batch
  no_batch=max(ID)
  no_gene_Chr_batch=no_batch
  for (ibatch in 1:no_gene_Chr_batch){ 
    gene_chr_batch=anno[ID==ibatch,]
    gene_chr_transcript_unique= unique(as.character(gene_chr_batch$transcript))
    no_transcript=length(gene_chr_transcript_unique)
    stop_points=sort(unique(gene_chr_batch$stop))
    
    if ((no_transcript > 1)&(length(stop_points)>1)) {
      
      # if there are more than 1 transcript annotation
      # find the breaking points
      possible_break_points=stop_points[1:(length(stop_points)-1)]
      valid_cut=possible_break_points*0
      
      # check possible_break_points
      for (ip in 1:length(possible_break_points)) {
        left_id = (gene_chr_batch$stop <= possible_break_points[ip])
        right_id = (gene_chr_batch$stop > possible_break_points[ip])
        batch_transcript=as.character(gene_chr_batch$transcript)
        after_break=length(unique(batch_transcript[left_id]))+length(unique(batch_transcript[right_id]))
        valid_cut[ip]=(after_break==no_transcript)
      }
      
      if (sum(valid_cut)>0){
        # cut into pieces
        valid_cut_points=possible_break_points[valid_cut>0]
        current_ID= gene_chr_batch$stop *0 + ibatch
        for (i in 1:length(current_ID)) {
          piece_number = sum(gene_chr_batch$stop[i] > valid_cut_points)
          if (piece_number>0) {
            current_ID[i]=piece_number+no_gene_Chr_batch
          } 
        }
        
        # update
        ID[ID==ibatch]=current_ID
        no_gene_Chr_batch=max(ID)
      } 
    }
  }
  return(ID)
  
}