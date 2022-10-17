FUN_ReadscRNA = function( InputFolder,Folder = Folder,
                      Path = "/monocle/outs/filtered_gene_bc_matrices/mm10",
                      list_files.df, Mode="10x" ,projectName=ProjectName) # Mode=c("10x","Exp")
  {

    if(Mode=="10x"){
      ## Read 10x files
      scRNA_SeuObj.list <- list()
      for(i in 1:nrow(list_files.df)){
        Folder <- list_files.df$Folder[i]
        Data.dgCMatrix <- Read10X(data.dir = paste0(InputFolder,"/", Folder,Path))
        Data.SeuObj <- CreateSeuratObject(counts = Data.dgCMatrix,
                                          project = projectName,
                                          min.cells = 3, min.features = 200)

        for (j in 1:(ncol(list_files.df)-1)) {
          Data.SeuObj@meta.data[[colnames(list_files.df)[j+1]]] <- rep(list_files.df[i,j+1],
                                                                       times=nrow(Data.SeuObj@meta.data))
        }

        scRNA_SeuObj.list[[i]] <- Data.SeuObj
        names(scRNA_SeuObj.list)[[i]] <- list_files.df$Folder[i]
      }
      rm(i,j,Folder,Data.dgCMatrix,Data.SeuObj)

    }else if(Mode=="H5AD") {
      scRNA_SeuObj.list <- list()
      for(i in 1:nrow(list_files.df)){
        # Install SeuratDisk # https://github.com/mojaveazure/seurat-disk
        if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes")}
        remotes::install_github("mojaveazure/seurat-disk")
        library(SeuratDisk)

        Folder <- list_files.df$Folder[i]
        ## Convert h5ad to h5seurat
        Convert(paste0(InputFolder,"/", Folder, "/scRNA.h5ad"), paste0(InputFolder,"/", Folder, "/scRNA.h5seurat"),
                assay = "RNA",) # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

        Data.SeuObj <- LoadH5Seurat(paste0(InputFolder,"/", Folder, "/scRNA.h5seurat")) # This .d5seurat object can then be read in manually

        for (j in 1:(ncol(list_files.df)-1)) {
          Data.SeuObj@meta.data[[colnames(list_files.df)[j+1]]] <- rep(list_files.df[i,j+1],
                                                                       times=nrow(Data.SeuObj@meta.data))
        }

        scRNA_SeuObj.list[[i]] <- Data.SeuObj
        names(scRNA_SeuObj.list)[[i]] <- list_files.df$Folder[i]
      }


    }else{

      ## Expression matrix
      ## GSE103322 HNSC
      scRNA_SeuObj.list <- list()
      for(i in 1:nrow(list_files.df)){
        Folder <- list_files.df$Folder[i]
        GeneExp.df <- read.table(paste0(InputFolder,"/", Folder, "/ExpMat.tsv"),
                                     header=T, row.names = 1, sep="\t")
        Data.SeuObj <- CreateSeuratObject(counts = GeneExp.df,
                                          project = projectName,
                                          min.cells = 3, min.features = 200)

        for (j in 1:(ncol(list_files.df)-1)) {
          Data.SeuObj@meta.data[[colnames(list_files.df)[j+1]]] <- rep(list_files.df[i,j+1],
                                                                       times=nrow(Data.SeuObj@meta.data))
        }

        scRNA_SeuObj.list[[i]] <- Data.SeuObj
        names(scRNA_SeuObj.list)[[i]] <- list_files.df$Folder[i]

        scAnno.df <- read.table(paste0(InputFolder,"/", Folder, "/CellMeta.tsv"),
                                     header=T, row.names = 1, sep="\t")
        # scRNA.SeuObj@meta.data[["sample"]] <- GeneExp.df[5,] %>% as.character()

          for (k in 1:nrow(scAnno.df)) {
            scRNA_SeuObj.list[[i]]@meta.data[[row.names(scAnno.df)[k]]] <- c(scAnno.df[k,] %>% as.character())
          }

      }
      rm(i,j,k,Folder,GeneExp.df,Data.SeuObj)

    }
    return(scRNA_SeuObj.list)
  }
