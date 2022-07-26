##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

#### Load the required libraries ####
  ## Check whether the installation of those packages is required from basic
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("Seurat")) install.packages("Seurat")

  ## Load Packages
  library(tidyverse)
  library(Seurat)

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if(!require("SeuratDisk", quietly = TRUE)) BiocManager::install("SeuratDisk")

## Load Packages
  library(SeuratDisk)

#### Load data ####
  # ## Old version for h5ad ##
  # #### Converse h5ad to Seurat ####
  # remotes::install_github("mojaveazure/seurat-disk")
  # library(SeuratDisk)
  #
  # # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
  # Convert("StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad", "PRJCA001063.h5seurat")
  #
  # # This .d5seurat object can then be read in manually
  # scRNA.SeuObj <- LoadH5Seurat("PRJCA001063.h5seurat")

  load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-05_CTAnno_singleR_RefPRJCA001063_PDAC.RData")
  ## Clean up the object
  rm(list=setdiff(ls(), c("scRNA.SeuObj")))

  Meta.df <- scRNA.SeuObj@meta.data %>% as.data.frame()

##### Current path and new folder setting*  #####
  ProjectName = "SeuratSmall"
  Version = paste0(Sys.Date(),"_",ProjectName,"_PADC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){dir.create(Save.Path)}

#### Data preprocessing ####
  # ## Try sample
  # scRNA.SeuObj@meta.data[["orig.ident"]] <- sample(c(1,2),nrow(scRNA.SeuObj@meta.data),replace = TRUE)
  # ## Reorder the col
  # scRNA.SeuObj@meta.data <- scRNA.SeuObj@meta.data %>%
  #                           relocate(orig.ident, .before = colnames(scRNA.SeuObj@meta.data)[1])%>%
  #                           relocate(n_counts, .before = colnames(scRNA.SeuObj@meta.data)[1]) %>%
  #                           relocate(n_genes, .before = colnames(scRNA.SeuObj@meta.data)[1]) %>%
  #                           rename(nFeature_RNA = n_genes, nCount_RNA = n_counts)
  scRNA.SeuObj@meta.data[["Cell_ID"]] <- rownames(scRNA.SeuObj@meta.data)
  scRNA.SeuObj@meta.data <- scRNA.SeuObj@meta.data %>%
                            relocate(Cell_ID, .before = colnames(scRNA.SeuObj@meta.data)[1])

#### Sampling ####
  # ## Test Sampling
  # Celltype_S100.set <- sample(scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique(), 100, replace = TRUE, prob = NULL)
  # scRNA_S50.SeuObj <- scRNA.SeuObj[,sample(1:100,50, replace = FALSE, prob = NULL)]
  # scRNA_S100.SeuObj <- scRNA.SeuObj[,sample(1:100,100, replace = FALSE, prob = NULL)]

  # ## Merging more than two seurat objects
  # ## Ref: https://github.com/satijalab/seurat/issues/706
  # ## Ref: https://mojaveazure.github.io/seurat-object/reference/Seurat-methods.html
  #
  # scRNA_Small.SeuObj_Merge <- merge(scRNA_S50.SeuObj, scRNA_S100.SeuObj) # merge(x,y,add.cell.ids = c(x@project.name,y@project.name))

  ## Extract SubSet
  # scRNA_Fib.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% c("Fibroblast cell")]

  # ## Old version
  # Set_SSize <- 50
  # CellType.set <- scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique()
  # for (i in 1:length(CellType.set)) {
  #
  #   scRNA_Small_Temp.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% CellType.set[i]]
  #   if(nrow(scRNA_Small_Temp.SeuObj@meta.data) < Set_SSize){
  #     scRNA_Small_Temp.SeuObj <- scRNA_Small_Temp.SeuObj[,sample(1:nrow((scRNA_Small_Temp.SeuObj@meta.data)),Set_SSize, replace = TRUE, prob = NULL)]
  #   }else{
  #     scRNA_Small_Temp.SeuObj <- scRNA_Small_Temp.SeuObj[,sample(1:nrow((scRNA_Small_Temp.SeuObj@meta.data)),Set_SSize, replace = FALSE, prob = NULL)]
  #   }
  #
  #   try({
  #     if(i==1){
  #       scRNA_Small.SeuObj <- scRNA_Small_Temp.SeuObj
  #     }else{
  #       scRNA_Small.SeuObj <- merge(scRNA_Small.SeuObj, scRNA_Small_Temp.SeuObj)
  #     }
  #   })
  #
  # }
  # rm(i,scRNA_Small_Temp.SeuObj)


  ## New version (Keep all data)
  Set_SSize <- 50
  Set_seed <- 1

  CellType.set <- scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique()
  for (i in 1:length(CellType.set)) {

    scRNA_Small_Temp.set <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% CellType.set[i]]$Cell_ID
    if(length(scRNA_Small_Temp.set) < Set_SSize){
      set.seed(Set_seed)
      scRNA_Small_Temp.set <- sample(scRNA_Small_Temp.set, Set_SSize, replace = TRUE, prob = NULL)
    }else{
      set.seed(Set_seed)
      scRNA_Small_Temp.set <- sample(scRNA_Small_Temp.set, Set_SSize, replace = FALSE, prob = NULL)

    }

    try({
      if(i==1){
        scRNA_Small.set <- scRNA_Small_Temp.set
      }else{
        scRNA_Small.set <- c(scRNA_Small.set,scRNA_Small_Temp.set)
      }
    })

  }
  scRNA_Small.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_ID"]] %in% scRNA_Small.set]
  rm(i,scRNA_Small.set,  scRNA_Small_Temp.set)

## Plot
  DimPlot(scRNA_Small.SeuObj, reduction = "umap",group.by = "seurat_clusters")
  DimPlot(scRNA_Small.SeuObj, reduction = "umap",group.by = "Cell_type")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")

  scRNA_Small_S1.SeuObj <- scRNA_Small.SeuObj
  # scRNA_Small_S2.SeuObj <- scRNA_Small.SeuObj
  # scRNA_Small_S3.SeuObj <- scRNA_Small.SeuObj


#### Save the RData ####
  # rm(list=setdiff(ls(), c("scRNA_Small.SeuObj","Version","Save.Path",
  # "scRNA_Small_S1.SeuObj", "scRNA_Small_S2.SeuObj", scRNA_Small_S3.SeuObj)))
  save.image(paste0(Save.Path,"/",Version,"_SmallData.RData"))

  scRNA.SeuObj <- scRNA_Small.SeuObj
  save.image(paste0(Save.Path,"/",Version,"_SmallDataS.RData"))
