

DE_result_filter <- function(res, dds,
                               symbol,
                               logFC = 1,
                               pvalue=0.05){
  # check the object
  
  if ( ! isClass(res,"DEseqResults2") ) stop("input should from DESeq2::results")
  
  # extract the normalized counts
  condition <- strsplit(res@elementMetadata$description[2],
                        split =" ", fixed = TRUE)[[1]][c(6,8)]
  select_cols <- rownames(colData(dds)[colData(dds)$condition %in% condition, ])
  counts_df <- as.data.frame(counts(dds, normalized=TRUE)[,select_cols])
  
  
  #for debug 
  #print(dim(counts_df))
  x <- as.data.frame(res)
  x <- x[,c(2,5,6)]
  x <- cbind(geneID=rownames(x), x)
  x <- cbind(x, counts_df)
  
  # add the symbol information
  symbol$SYMBOL <- paste0("\"", symbol$SYMBOL, "\"")
  x <- merge(x, symbol, by.x="geneID", by.y="geneID", all.x=TRUE)
  #print(head(x))
  x <- x[ ! is.na(x$padj) & ! is.na(x$log2FoldChange), ]
  res_flt <- x[x$padj < pvalue & abs(x$log2FoldChange) > logFC ,]
  #print(head(res_flt))
  
  res_flt

}



DE_result_filter_download <- function(outdir = "."){
  
  
  write.csv(x, file=file.path(outdir, "DESeq2-Analysis-results.csv"), row.names = FALSE )
  for ( i in c(logFC1, logFC2)){
    for (j in c(pvalue1, pvalue2)){
      file_name_up <- paste0("DESeq2-logFC-",i,"-pvalue-",j,"-up.csv")
      file_name_down <- paste0("DESeq2-logFC-",i,"-pvalue-",j,"-down.csv")
      flt_up <- subset(x, pvalue < j & log2FoldChange > i)
      flt_down <- subset(x, pvalue < j & log2FoldChange < -i)
      write.csv(flt_up, file=file.path(outdir,file_name_up), 
                row.names = FALSE)
      write.csv(flt_down, file=file.path(outdir,file_name_down), 
                row.names = FALSE)
      
    }
  }
  
  
}
