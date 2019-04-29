
# Build DESeqDataSeq object from matrix with given contrl and case group
# and remove rows in which there are very few reads

BuildDdsFromMatrix <- function(mt, 
                               control, case,
                               min_count = 10){
  
  contrl_rep <- length(control)
  case_rep <- length(case)
  condition <- rep(c("control","case"), times = c(contrl_rep, case_rep))
  condition <- factor(condition, levels = c("control", "case"))
  coldata <- data.frame(condition = condition)
  rownames(coldata) <- c(control, case)
  matrix <- mt[,c(control, case)]
  
  
  dds <- DESeq2::DESeqDataSetFromMatrix(matrix, 
                                        coldata, 
                                        design = ~ condition)
  
  # Pre-filtering
  keep <- rowSums(DESeq2::counts( dds )) >= min_count
  dds <-  dds[keep, ]
  
  
  return(dds)
}