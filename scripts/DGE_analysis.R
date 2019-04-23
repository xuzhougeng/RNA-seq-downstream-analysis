
# DGE analysis

DGE_analysis <- function(dds, 
                         lfcThreshold = 0,
                         pAdjustMethod = "fdr",
                         ...){
  
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, 
                         lfcThreshold = lfcThreshold,
                         pAdjustMethod = pAdjustMethod)
  return(list(dds, res))
  
}