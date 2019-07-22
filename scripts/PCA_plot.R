# Principal component plot of the samples

PCA_plot <- function(dds){
  vsd <- DESeq2::vst(dds, blind = FALSE)
  p <- DESeq2::plotPCA(vsd,  intgroup = "condition" ) +
    ggrepel::geom_label_repel(ggplot2::aes(label= name)) +
    ggplot2::theme_bw()
  p
}