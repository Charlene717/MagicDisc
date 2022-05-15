BeautifyVennDiag <- function(CCMarker_Male.lt, CCMarker_Female.lt, CellType.list,
                             Path = PathBiomarkers, projectName = ProjectName,
                             sampletype = Sampletype){

  ##-------------- Venn Pos --------------##
  source("FUN_Venn.R")
  # pdf(file = paste0(Path,"/PBMC_Female_VolcanoPlot.pdf"),width = 7, height = 7 )
  Sep_Cla3_Venn.Path <- paste0(sampletype,"_",projectName,"_Separate_","_VennDiagrame")
  dir.create(paste0(Path,"/",Sep_Cla3_Venn.Path))
  Venn_CCMarker_Pos <- list()
  for(i in c(1:length(CellType.list))){
    try({
      Venn_CCMarker_Pos[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][[paste0(projectName, "Marker.S_Pos_List")]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][[paste0(projectName, "Marker.S_Pos_List")]],
                                               CellType.list[i],"Pos","#9d0208","#f08080",SampleType = sampletype,
                                               PathName = paste0(Path,"/",Sep_Cla3_Venn.Path),
                                               ClassSet3_1 = ClassSet3.set[1], ClassSet3_2 =ClassSet3.set[2])
      names(Venn_CCMarker_Pos)[[i]] <- paste0("Venn_",projectName,"Marker.",CellType.list[i],"_Pos")
    })
  }
  rm(i)

  ##-------------- Venn Neg --------------##
  Venn_CCMarker_Neg <- list()
  for(i in c(1:length(CellType.list))){
    try({
      Venn_CCMarker_Neg[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][[paste0(projectName, "Marker.S_Neg_List")]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][[paste0(projectName, "Marker.S_Neg_List")]],
                                               CellType.list[i],"Neg","#00296b","#1368aa",SampleType = sampletype,
                                               PathName = paste0(Path,"/",Sep_Cla3_Venn.Path))

      names(Venn_CCMarker_Neg)[[i]] <- paste0("Venn_",projectName,"Marker.",CellType.list[i],"_Neg")
    })
  }
  rm(i,Sep_Cla3_Venn.Path)

  Venn_CCMarke.lt <- c(Venn_CCMarker_Pos, Venn_CCMarker_Neg)
  names(Venn_CCMarke.lt)[[1]] <- paste0("CCMarker_Pos")
  names(Venn_CCMarke.lt)[[2]] <- paste0("CCMarker_Neg")
  return(Venn_CCMarke.lt)

}
