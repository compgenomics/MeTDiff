
# read gtf file
.read.gtf <- function(PARAMETERS){
  # download the annotation
  if ((!is.na(PARAMETERS$GENOME)) & (!is.na(PARAMETERS$UCSC_TABLE_NAME)) & is.na(PARAMETERS$GENE_ANNO_GTF) & is.na(PARAMETERS$TXDB)) {
    op <- options(warn = (-1))
    txdb =makeTxDbFromUCSC(genome=PARAMETERS$GENOME,
                                   tablename=PARAMETERS$UCSC_TABLE_NAME)
    options(op)
  }
  
  # use provided annotation data file
  if (!is.na(PARAMETERS$GENE_ANNO_GTF) & is.na(PARAMETERS$TXDB) ) {
    op <- options(warn = (-1))
    txdb=makeTxDbFromGFF(PARAMETERS$GENE_ANNO_GTF,format="gtf")
    options(op)
  }
  
  # use provided annotation data file
  if (!is.na(PARAMETERS$TXDB) ) {
    txdb=PARAMETERS$TXDB
  }
  
  # try the internet method
  op <- options(warn = (-1))
  ID = keys(txdb, "TXID")
  # temp = select(txdb, ID , c(cols(txdb))[c(9:12,13,16)], "TXID")
  temp = select(txdb, ID , c(columns(txdb))[c(7:8,12:14)], "TXID")
  options(op)

  # get the anno
  temp = cbind(temp,"exon")
  temp=temp[,c(2,7,4,5,3,6,1)]
  colnames(temp)=c("chr","feature","start","stop","strand","gene","transcript")
  gtf=temp
  
  # return data
  return(gtf)
}
