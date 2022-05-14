Beautify_Heatmap_Seurat <- function(scRNA.SeuObj, PBMC.markers, topN = 7,
                                    Path = PathCluster, 
                                    projectName = ProjectName
                                    ) {
  
  
  
  # https://github.com/satijalab/seurat/issues/2960
  # Filter the top markers and plot the heatmap
  top_NSet = topN
  PBMC.markers %>%
    group_by(cluster) %>%
    top_n(n = top_NSet, wt = avg_log2FC) -> top_N
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  DoHeatmap(scRNA.SeuObj, features = top_N$gene) + NoLegend()
  
  
  ## Export pdf
  pdf(
    file = paste0(Path, "/",projectName,"_Heatmap_Cluster_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  HM <- DoHeatmap(scRNA.SeuObj, features = top_N$gene,size = 2,angle = 60) +
        scale_fill_gradient2(low="#5283ff",mid ="white", high ="#ff5c5c") +
        theme(axis.text.y = element_text(size  = 5)) +
        theme(legend.position = "bottom" )
  
  print(HM)
  dev.off()
  
  ## Export txt
  write.table(data.frame(Num = row.names(top_N),top_N), 
              file=paste0(Path, "/",projectName,"_ClusterMarker_top",top_NSet,"Gene.txt"),sep="\t", 
              row.names=T, quote = FALSE)
  write.table(data.frame(Gene = row.names(PBMC.markers), PBMC.markers), 
              file=paste0(Path, "/",projectName,"_ClusterMarker_AllGene.txt"),sep="\t", 
              row.names=F,quote = FALSE)
  
  return(scRNA.SeuObj)
}
