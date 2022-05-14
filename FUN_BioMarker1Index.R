BioMarker1Index <- function(scRNA.SeuObj, 
                      Path = PathBiomarkers, projectName = ProjectName,
                      sampletype = Sampletype, cellType.list = CellType.list,
                      classSet2 = ClassSet2, 
                      Type = paste0("celltype.",classSet2) ){

  ##### 08_3 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  
  Idents(scRNA.SeuObj) <- Type
  Sep_Cla3_FMar.Path <- paste0(Sampletype,"_",ProjectName,"_Pooled_FindMarkers")
  dir.create(paste0(Path,"/",Sep_Cla3_FMar.Path))
  
  CCMarker_SPA.lt <- list()
  for(i in c(1:length(cellType.list))){
    try({
      CCMarker_SPA.lt[[i]] <- Find_Markers(scRNA.SeuObj,
                                           paste0(cellType.list[i],"_",classSet2.set[1]),
                                           paste0(cellType.list[i],"_",classSet2.set[2]),
                                           cellType.list[i],
                                           Path = Path,
                                           ResultFolder =  Sep_Cla3_FMar.Path,
                                           ProjectTitle = ProjectName)
      
      # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",cellType.list[i])
      names(CCMarker_SPA.lt)[[i]] <- paste0(cellType.list[i])
    })
  }
  rm(i,Sep_Cla3_FMar.Path)
  
  CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]
  
  
  ## Generate pdf and tif file for VolcanoPlot
  Sep_Cla3_Volcano.Path <- paste0(Sampletype,"_",ProjectName,"_Pooled_","_VolcanoPlot")
  dir.create(paste0(Path,"/",Sep_Cla3_Volcano.Path))
  pdf(file = paste0(Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,".pdf"),width = 7, height = 7 )
  for (i in 1:length(cellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S")]],
                        CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Pos_List")]],
                        CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Neg_List")]], ShowGeneNum = 6)+
              ggtitle(paste0(Sampletype,"_",cellType.list[i]))
      )
    })
  }
  # graphics.off()
  dev.off()
  rm(i)
  
  for (i in 1:length(cellType.list)) {
    try({
      tiff(file = paste0(Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,"_",cellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S")]],
                        CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Pos_List")]],
                        CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Neg_List")]])+
              ggtitle(paste0(Sampletype,"_",cellType.list[i]))
      )
      
      graphics.off()
    })
  }
  rm(i,Sep_Cla3_Volcano.Path)
  
return(CCMarker.lt)

}


# ##### Export marker gene from specific cluster #####
#   # For performing differential expression after integration, we switch back to the original data
#   set.seed(1) # Fix the seed
#   DefaultAssay(scRNA.SeuObj) <- "RNA"
#
#   # nk.markers <- FindConservedMarkers(scRNA.SeuObj, ident.1 = 6, grouping.var = "sample", verbose = FALSE)
#   library(BiocManager)
#   library(multtest)
#   nk.markers <- FindConservedMarkers(scRNA.SeuObj, ident.1 = 'NK', grouping.var = "sample", verbose = FALSE)
#   head(nk.markers)
#
#   rm(nk.markers)

# ##### Identify differential expressed genes across conditions  #####
#   library(ggplot2)
#   library(cowplot)
#   theme_set(theme_cowplot())
#   CD4T.cells <- subset(scRNA.SeuObj, idents = "CD4+T")
#   Idents(CD4T.cells) <- "Cachexia"
#   avg.CD4T.cells <- as.data.frame(log1p(AverageExpression(CD4T.cells, verbose = FALSE)$RNA))
#   avg.CD4T.cells$gene <- rownames(avg.CD4T.cells)
#
#   MacrophageM2 <- subset(scRNA.SeuObj, idents = "Mac2")
#   Idents(MacrophageM2) <- "Cachexia"
#   avg.MacrophageM2 <- as.data.frame(log1p(AverageExpression(MacrophageM2, verbose = FALSE)$RNA))
#   avg.MacrophageM2$gene <- rownames(avg.MacrophageM2)
#
#   genes.to.label = c("Sox17", "Mrpl15", "Lypla1", "Tcea1", "Rgs20", "Atp6v1h", "Rb1cc1", "4732440D04Rik", "St18")
#   p1 <- ggplot(avg.CD4T.cells, aes(EO, LO)) + geom_point() + ggtitle("Cachexia T.cells")
#   p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
#   p2 <- ggplot(avg.MacrophageM2, aes(EO, LO)) + geom_point() + ggtitle("Cachexia Macrophage")
#   p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
#   p1 + p2
#   rm(p1 , p2 ,CD4T.cells, MacrophageM2, avg.CD4T.cells, avg.MacrophageM2)

