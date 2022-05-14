Export_All_DRPlot <- function(scRNA.SeuObj, 
                              Path1 = PathCellType, Path2 = PathCluster,
                              projectName = ProjectName){
  ## UMAP tSNE
  DimPlot(scRNA.SeuObj, label = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5))
  
  DimPlot(scRNA.SeuObj,group.by = "celltype",label.size = 7,label = TRUE,
          pt.size =2) %>% BeautifyUMAP(FileName = paste0(Path1,"/",ProjectName,"_nlDR_CellType"))
  DimPlot(scRNA.SeuObj,group.by = colnames(list_files.df)[2],
          pt.size =0.5) %>% BeautifyUMAP(FileName = paste0(Path1,"/",ProjectName,"_nlDR_",colnames(list_files.df)[2]))
  DimPlot(scRNA.SeuObj,group.by = "seurat_clusters",label.size = 7, label = TRUE,
          pt.size =1) %>% BeautifyUMAP(FileName = paste0(Path2, "/",ProjectName,"_nlDR_Clusters"))
  
  
  pdf(
    file = paste0(Path1,"/",ProjectName,"_nlDR_CellType_Sup.pdf"),
    width = 12,  height = 8
  )
  
  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>%
    BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)
  
  
  for (i in 1:(length(Ori_Meta.set)-3)) {
    print(DimPlot(scRNA.SeuObj, reduction = "umap", group.by = Ori_Meta.set[i+3]) %>%
            BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+
            theme(plot.title = element_text(vjust = 0.85)))
    print(DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2, split.by = Ori_Meta.set[i+3], label = TRUE, label.size = 4) %>%
            BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                           SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))
    
    print(DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = Ori_Meta.set[i+3]) %>%
            BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+
            theme(plot.title = element_text(vjust = 0.85)))
    print(DimPlot(scRNA.SeuObj, reduction = "tsne", ncol = 2, split.by = Ori_Meta.set[i+3], label = TRUE, label.size = 4) %>%
            BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                           SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))
    
  }
  rm(i)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
}