BeautifyVennDiag <- function(CCMarker_Male.lt, CCMarker_Female.lt, CellType.list,list_files.df,
                             Path = PathBiomarkers, projectName = ProjectName,
                             sampletype = Sampletype, classSet3 = ClassSet3){


  source("FUN_Venn.R")
  # pdf(file = paste0(Path,"/PBMC_Female_VolcanoPlot.pdf"),width = 7, height = 7 )
  ## Create new folder
  Sep_Cla3_Venn.Path <- paste0(Path,"/",sampletype,"_",projectName,"_Separate","_VennDiagrame")
  if (!dir.exists(Sep_Cla3_Venn.Path)){
    dir.create(Sep_Cla3_Venn.Path)
  }

  ##-------------- Venn Pos --------------##

  ClassSet3.set <- list_files.df[[classSet3]] %>% unique()
  Venn_CCMarker_Pos <- list()
  for(i in c(1:length(CellType.list))){
    try({
      Venn_CCMarker_Pos[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][[paste0(projectName, "Marker.S_Pos_List")]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][[paste0(projectName, "Marker.S_Pos_List")]],
                                               CellType.list[i],"Pos","#9d0208","#f08080",SampleType = sampletype,
                                               PathName = Sep_Cla3_Venn.Path,
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
                                               PathName = Sep_Cla3_Venn.Path,
                                               ClassSet3_1 = ClassSet3.set[1], ClassSet3_2 =ClassSet3.set[2])

      names(Venn_CCMarker_Neg)[[i]] <- paste0("Venn_",projectName,"Marker.",CellType.list[i],"_Neg")
    })
  }
  rm(i,Sep_Cla3_Venn.Path)

  Venn_CCMarke.lt <- list(Venn_CCMarker_Pos, Venn_CCMarker_Neg)
  names(Venn_CCMarke.lt)[[1]] <- paste0("CCMarker_Pos")
  names(Venn_CCMarke.lt)[[2]] <- paste0("CCMarker_Neg")
  return(Venn_CCMarke.lt)

}
