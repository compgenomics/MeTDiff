

.bed2grangeslist <- function(content) {
# read bed file
# a = read.table(filepath,sep="\t",header=FALSE,stringsAsFactors =FALSE)
  op <- options(warn = (-1))
  
a = content
# get transcripts
no_tx = length(a[,1])
tx_id = 1:no_tx;
tx_name = as.character(1:no_tx)
tx_chrom = a[,1]
tx_strand = a[,6]
tx_start = a[,2]+1
tx_end = a[,3]
transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
head(transcripts)

# get splicings
splicing = data.frame()
for (i in 1:no_tx) {
  tx = a[i,]
  tx_id = i
  exon_rank=1:as.integer(tx[10])
  
  # get start
  temp = as.integer(strsplit(as.character(tx[12][1,1]), ",")[[1]]) + tx_start[i]
  exon_start=temp
  
  # get end
  temp = as.integer(strsplit(as.character(tx[11][1,1]), ",")[[1]])
  temp2 = temp + exon_start - 1
  exon_end=temp2
  
  # get CDS
  cds_start = exon_start
  cds_end = exon_end
  
  # get data frame
  splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
  
  # collect result
  splicing = rbind(splicing, splicing_tx)
}

# get genes
tx_name = tx_name
gene_id = as.character(a[,4])
gene_id[is.na(gene_id)]="NA"
gene=data.frame(tx_name,gene_id)

# 

peaks = makeTranscriptDb(transcripts=transcripts, splicings=splicing,
                 genes=gene)
peaks= exonsBy(peaks, by="tx")
return(peaks)
  
  options(op)
}


