BioMarker <- function(scRNA.SeuObj, 
                      Path = PathBiomarkers, projectName = ProjectName,
                      sampletype = Sampletype,
                      classSet2 = ClassSet2, classSet3 = ClassSet3,
                      Type = "celltype" ){

##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  #### Define group by different phenotype ####
  source("FUN_Find_Markers.R")
  # scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
  
  # for (i in 1:(ncol(list_files.df)-1)) {
  #   scRNA.SeuObj[[paste0("celltype.",colnames(list_files.df)[1+i])]] <- paste(Idents(scRNA.SeuObj), as.matrix(scRNA.SeuObj[[colnames(list_files.df)[1+i]]]), sep = "_")
  #
  # }
  
  scRNA.SeuObj[[paste0(Type,".",classSet2)]] <- paste(Idents(scRNA.SeuObj),
                                                         as.matrix(scRNA.SeuObj[[classSet2]]), sep = "_")
  scRNA.SeuObj[[paste0(Type,".",classSet2,".",classSet3)]] <- paste(Idents(scRNA.SeuObj),
                                                                       as.matrix(scRNA.SeuObj[[classSet2]]),
                                                                       as.matrix(scRNA.SeuObj[[classSet3]]), sep = "_")
  Idents(scRNA.SeuObj) <- paste0(Type,".",classSet2,".",classSet3)
  
  
  DefaultAssay(scRNA.SeuObj) <- "RNA"
  
  
  classSet2.set <- list_files.df[[classSet2]] %>% unique()
  classSet3.set <- list_files.df[[classSet3]] %>% unique()
  CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[[Type]]))
  # CellType.list <- CellType.list[-9]
  ####-------------- Find Marker gene in Male --------------####
  Sep_Cla3_FMar.Path <- paste0(sampletype,"_",projectName,"_Separate_",classSet3.set[1],"_FindMarkers")
  dir.create(paste0(Path,"/",Sep_Cla3_FMar.Path))
  
  # About 15 mins
  CCMarker_Male.lt <- list()
  for(i in c(1:length(CellType.list))){
    try({
      CCMarker_Male.lt[[i]] <- Find_Markers(scRNA.SeuObj,
                                            paste0(CellType.list[i],"_",classSet2.set[1],"_",classSet3.set[1]),
                                            paste0(CellType.list[i],"_",classSet2.set[2],"_",classSet3.set[1]),
                                            CellType.list[i],
                                            Path = Path,
                                            ResultFolder = paste0(Sep_Cla3_FMar.Path),
                                            ProjectTitle = projectName)
      # names(CCMarker_Male.lt)[[i]] <- paste0("CCMarker_Male.lt.",CellType.list[i])
      names(CCMarker_Male.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)
  
  CCMarker_Male.lt <- CCMarker_Male.lt[!unlist(lapply(CCMarker_Male.lt,is.null))]
  
  
  ## Generate pdf and tif file for Male VolcanoPlot
  Sep_Cla3_Volcano.Path <- paste0(sampletype,"_",projectName,"_Separate_",classSet3.set[1],"_VolcanoPlot")
  dir.create(paste0(Path,"/",Sep_Cla3_Volcano.Path ))
  
  pdf(file = paste0(Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,".pdf"),
      width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_Male.lt[[i]][[paste0(projectName, "Marker.S")]],
                        CCMarker_Male.lt[[i]][[paste0(projectName, "Marker.S_Pos_List")]],
                        CCMarker_Male.lt[[i]][[paste0(projectName, "Marker.S_Neg_List")]], ShowGeneNum = 6)+
              ggtitle(paste0(projectName,"_",classSet3.set[1],"_",CellType.list[i]))
      )
    })
  }
  dev.off() # graphics.off()
  rm(i)
  
  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,"_",CellType.list[i],".tif"),
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Male.lt[[i]][[paste0(projectName, "Marker.S")]],
                        CCMarker_Male.lt[[i]][[paste0(projectName, "Marker.S_Pos_List")]],
                        CCMarker_Male.lt[[i]][[paste0(projectName, "Marker.S_Neg_List")]])+
              ggtitle(paste0(projectName,"_",classSet3.set[1],"_",CellType.list[i]))
      )
      
      graphics.off()
    })
  }
  rm(i,Sep_Cla3_FMar.Path)
  
  ####-------------- Find Marker gene in Female --------------####
  Sep_Cla3_FMar.Path <- paste0(sampletype,"_",projectName,"_Separate_",classSet3.set[2],"_FindMarkers")
  dir.create(paste0(Path, "/", Sep_Cla3_FMar.Path))
  
  # About 15 mins
  CCMarker_Female.lt <- list()
  
  for(i in c(1:length(CellType.list))){
    try({
      CCMarker_Female.lt[[i]] <- Find_Markers(scRNA.SeuObj,
                                              paste0(CellType.list[i],"_",classSet2.set[1],"_",classSet3.set[2]),
                                              paste0(CellType.list[i],"_",classSet2.set[2],"_",classSet3.set[2]),
                                              CellType.list[i],
                                              Path = Path,
                                              ResultFolder = paste0(Sep_Cla3_FMar.Path),
                                              ProjectTitle = projectName)
      # names(CCMarker_Female.lt)[[i]] <- paste0("CCMarker_Female.lt.",CellType.list[i])
      names(CCMarker_Female.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)
  
  CCMarker_Female.lt <- CCMarker_Female.lt[!unlist(lapply(CCMarker_Female.lt,is.null))]
  
  
  ## Generate pdf and tif file for Female VolcanoPlot
  Sep_Cla3_Volcano.Path <- paste0(sampletype,"_",projectName,"_Separate_",classSet3.set[2],"_VolcanoPlot")
  dir.create(paste0(Path,"/",Sep_Cla3_Volcano.Path ))
  
  # dir.create(paste0(Path,"/PBMC_SSA_Female_VolcanoPlot/"))
  
  pdf(file = paste0(Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,".pdf"),
      width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_Female.lt[[i]][[paste0(projectName, "Marker.S")]],
                        CCMarker_Female.lt[[i]][[paste0(projectName, "Marker.S_Pos_List")]],
                        CCMarker_Female.lt[[i]][[paste0(projectName, "Marker.S_Neg_List")]], ShowGeneNum = 6)+
              ggtitle(paste0(projectName,"_",classSet3.set[2],"_",CellType.list[i]))
      )
    })
  }
  dev.off() # graphics.off()
  rm(i)
  
  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,"_",CellType.list[i],".tif"),
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Female.lt[[i]][[paste0(projectName, "Marker.S")]],
                        CCMarker_Female.lt[[i]][[paste0(projectName, "Marker.S_Pos_List")]],
                        CCMarker_Female.lt[[i]][[paste0(projectName, "Marker.S_Neg_List")]])+
              ggtitle(paste0(projectName,"_",classSet3.set[2],"_",CellType.list[i]))
      )
      
      graphics.off()
    })
  }
rm(i,Sep_Cla3_FMar.Path)

return()

}