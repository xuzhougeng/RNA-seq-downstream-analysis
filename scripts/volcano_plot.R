##' volcano plot 
##'
##' 
##' @title volcano_plot
##' @param df differential analysis result
##' @param selectgenes, data.frame with gene ID and gene symbol
##' @param circlesize the circle size
##' @param strocksize the strock size
##' @param log2FC1 log2 Fold Change cut-off 1
##' @param log2FC2 log2 Fold change cut-off 2
##' @param pval1 p value cut-off 1, soft
##' @param pval2 p value cut-off 2, middle
##' @param pval3 p value cut-off 3, hard
##' @export
##' @author zhougeng xu, haitao
volcano_plot <- function(df, 
                             selectgenes = "", 
							 circlesize=2.0,
							 strocksize=1.0,
                             log2FC1 = 1.5, log2FC2 = 2.5, log2FC3 = 5,
                             pval1 = 0.05, pval2 = 1e-4, pval3 = 1e-5
                             ){
  require(grid)
  mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13",
             "#088247","#58CDD9","#7A142C","#5D90BA","#431A3D",
             "#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
  
  x <- subset(df, ! is.na(df$pvalue))
  
  # set the reang of x-axis and y-axis
  xmin <- (range(x$log2FoldChange)[1]- (range(x$log2FoldChange)[1]+ 10))
  xmax <- (range(x$log2FoldChange)[1]+ (10-range(x$log2FoldChange)[1]))
  ymin <- 0
  ymax <- -log10(sort(x$pvalue))[3] * 1.1
  
  # set dot color for different log fold change
  n1 <- length(x[, 1])
  cols <- rep("grey", n1)
  names(cols)<- rownames(x)
  
  cols[x$pvalue < pval1 & x$log2FoldChange > log2FC1]<- "#FB9A99"
  cols[x$pvalue < pval2 & x$log2FoldChange > log2FC2]<- "#ED4F4F"
  cols[x$pvalue < pval1 & x$log2FoldChange < -log2FC1]<- "#B2DF8A"
  cols[x$pvalue < pval2 & x$log2FoldChange < -log2FC2]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # set dot size  for different log fold change
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  size[x$pvalue < pval1 & x$log2FoldChange > log2FC1]<- 2
  size[x$pvalue < pval2 & x$log2FoldChange > log2FC2]<- 4
  size[x$pvalue < pval3 & x$log2FoldChange > log2FC3]<- 6
  size[x$pvalue < pval1 & x$log2FoldChange < -log2FC1]<- 2
  size[x$pvalue < pval2 & x$log2FoldChange < -log2FC2]<- 4
  size[x$pvalue < pval3 & x$log2FoldChange < -log2FC3]<- 6
  
  
  # ggplot2: backgroud
  p1 <- ggplot2::ggplot(data=x, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
    ggplot2::geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
    ggplot2::labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") +
    ggplot2::ylim(c(ymin,ymax)) + 
    ggplot2::scale_x_continuous(
      breaks = c(-10, -5, -log2FC1, 0, log2FC1, 5, 10), 
      labels = c(-10, -5, -log2FC1, 0, log2FC1, 5, 10),
      limits = c(xmin, xmax) 
  )  +
  ggplot2::geom_vline(xintercept = c(-log2FC1, log2FC1), color="grey40", 
             linetype="longdash", lwd = 0.5) + 
  ggplot2::geom_hline(yintercept = -log10(pval1), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  ggplot2::theme_bw(base_size = 12, base_family = "Times") +
  ggplot2::theme(panel.grid=ggplot2::element_blank())
  # add the threshold line
  p1 <- p1 + 
  ggplot2::geom_vline(xintercept = c(-log2FC2, log2FC2), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  ggplot2::geom_hline(yintercept = -log10(pval2), color="grey40", 
             linetype="longdash", lwd = 0.5)
  
  #ggplot2: highlight the important genes
  if (nchar(selectgenes) < 1){
 	return(p1) 
  } else{
    genes <- strsplit(selectgenes, ",")[[1]]
  	selectgenes <- x[x$geneID == genes,]
  	p2 <- p1 + 
  	# add black circle in the interest gene
  	ggplot2::geom_point(data = selectgenes, 
             alpha = 1, size = circlesize, shape = 1, 
             stroke = strocksize, #size of circle
             color = "black") +
  	 # show the gene name
   	 ggplot2::scale_color_manual(values = mycol) + 
   	 ggrepel::geom_text_repel(data = selectgenes, ggplot2::aes(label=geneID),
                  show.legend = FALSE, #不显示图例
                  size = 5, box.padding = unit(0.35, "lines"), 
   	     point.padding = unit(0.3, "lines")) +
  	ggplot2::guides(color=ggplot2::guide_legend(title = NULL)) 
  
  	return(p2)
  }
  
}
