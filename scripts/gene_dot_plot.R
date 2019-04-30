gene_dot_plot <- function(dds, gene){
  
  d <- plotCounts(dds, 
                  gene=gene, 
                  intgroup="condition",
                  normalized = TRUE,
                  transform = TRUE,
                  returnData = TRUE)
  p <- ggplot(d, 
              aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks=c(25,100,400)) +
    ylab("normalized and log10 transformed count")
  
  return(p)
  
}