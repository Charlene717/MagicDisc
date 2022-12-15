MetaSummary = function(SeuObj_BFQC.list = scRNA_SeuObj.list,  scRNA_BFQC.SeuObj = scRNA_Ori.SeuObj,
                       SeuObj_AFTQC.list = scRNA_SeuObj_QC.list, scRNA_AFTQC.SeuObj = scRNA.SeuObj,
                       Path = Save.Path,
                       Folder = "B01_CellCount", projectName = ProjectName){

  # Create new folder
  NewPath <- paste0(Path,"/",Folder)
  if (!dir.exists(NewPath)){
    dir.create(NewPath)
  }


  Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
  colnames(Meta.df) <- c("Folder","Cell_Num","Gene_Num")

  ## Before QC
  for (i in 1:length(SeuObj_BFQC.list)) {
    Meta.df[i,1] <- names(SeuObj_BFQC.list[i])
    Meta.df[i,2] <- ncol(SeuObj_BFQC.list[[i]]@assays[["RNA"]]@counts)
    Meta.df[i,3] <- nrow(SeuObj_BFQC.list[[i]]@assays[["RNA"]]@counts)

  }

  # Summary to Meta table
  Meta.df[i+1,1] <- c("Summary")
  Meta.df[i+1,2] <- ncol(scRNA_BFQC.SeuObj@assays[["RNA"]]@counts)
  Meta.df[i+1,3] <- nrow(scRNA_BFQC.SeuObj@assays[["RNA"]]@counts)

  ## After QC
  for (j in 1:length(SeuObj_AFTQC.list)) {
    Meta.df[i+j+1,1] <- paste0(names(SeuObj_AFTQC.list[j]),".QC")
    Meta.df[i+j+1,2] <- ncol(SeuObj_AFTQC.list[[j]]@assays[["RNA"]]@counts)
    Meta.df[i+j+1,3] <- nrow(SeuObj_AFTQC.list[[j]]@assays[["RNA"]]@counts)

  }

  # Summary to Meta table
  Meta.df[i+j+2,1] <- c("Summary.QC")
  Meta.df[i+j+2,2] <- ncol(scRNA_AFTQC.SeuObj@assays[["RNA"]]@counts)
  Meta.df[i+j+2,3] <- nrow(scRNA_AFTQC.SeuObj@assays[["RNA"]]@counts)


  rm(i,j)

  Meta.df <- left_join(Meta.df,list_files.df, by = "Folder")
  Meta.df[is.na(Meta.df)] <- ""

  ## Export data
  write.table( Meta.df ,
               file = paste0(NewPath,"/",projectName,"_CellCount_Meta.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )


  return(Meta.df)
}
