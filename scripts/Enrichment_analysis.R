

# GO analysis -------------------------------------------------------------
enrich_analysis <- function(res,
                               gene_type = "all",
                               ont = "BP",
                               pAdjustMethod = "BH",
                               minGSSize = 10, 
                               maxGSSize = 500
                               
){
  if (gene_type == "up"){
    geneID <- res[res$log2FoldChange > 0, "geneID"]
  } else if( gene_type == "down" ){
    geneID <- res[res$log2FoldChange < 0, "geneID"]
  } else{
    geneID = res$geneID
  }
  
  
  if (ont %in% c("BP", "CC", "MF")){
    eout <- enrichGO( geneID ,
                     OrgDb = org.At.tair.db::org.At.tair.db,
                     keyType = "TAIR",
                     ont = "BP")
  } else{
    eout <- enrichKEGG( geneID ,
                     organism = "ath",
                     keyType = "kegg")
  }
  #print(eout)
  eout
  
}


# Enrichment Plot ---------------------------------------------------------
enrich_plot <- function(eout , showCategory = 20){
  p <- dotplot(eout, showCategory = showCategory)
  p
}


# GSEA Enrichment ----------------------------------------------------------

GSEA_analysis <- function(res, 
                          rank_type,
                          ont = "BP",
                          pAjustMethod = "BH"){
  
  res <- as.data.frame(res)
  geneList <- res[, rank_type]
  names(geneList) <- as.character(row.names(res))
  
  
  geneList <- geneList[!is.na(geneList)]
  geneList <- sort(geneList, decreasing = T)
  
  if (ont %in% c("BP", "CC", "MF")) {
    
    gseaout <- gseGO(geneList,
                     ont = ont, 
                     OrgDb = org.At.tair.db::org.At.tair.db,
                     keyType = "TAIR",
                     pAdjustMethod = "BH"
                     )  
    
  } else {
    gseaout <- gseKEGG(geneList, 
                       organism = "ath",
                       keyType = "kegg")    
  }
  #print(gseaout)
  gseaout
}


GSEA_plot <- function(gseaout, selected_row){
  
  df <- as.data.frame(gseaout)
  
  geneSetID <- rownames(gseaout[selected_row,])
  
  p <- enrichplot::gseaplot2(gseaout, geneSetID)
  p
}

