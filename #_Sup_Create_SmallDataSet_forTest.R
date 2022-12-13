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

##### Current path and new folder setting*  #####
  ProjectName = "SeuratSmall"
  Version = paste0(Sys.Date(),"_",ProjectName,"_PADC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){dir.create(Save.Path)}

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

#### Data preprocessing ####
  # ## Try sample
  # scRNA.SeuObj@meta.data[["orig.ident"]] <- sample(c(1,2),nrow(scRNA.SeuObj@meta.data),replace = TRUE)
  # ## Reorder the col
  # scRNA.SeuObj@meta.data <- scRNA.SeuObj@meta.data %>%
  #                           relocate(orig.ident, .before = colnames(scRNA.SeuObj@meta.data)[1])%>%
  #                           relocate(n_counts, .before = colnames(scRNA.SeuObj@meta.data)[1]) %>%
  #                           relocate(n_genes, .before = colnames(scRNA.SeuObj@meta.data)[1]) %>%
  #                           rename(nFeature_RNA = n_genes, nCount_RNA = n_counts)

#### Sampling ####
  # ## Test Sampling
  # Celltype_S100.set <- sample(scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique(), 100, replace = TRUE, prob = NULL)
  # scRNA_S50.SeuObj <- scRNA.SeuObj[,sample(1:100,50, replace = FALSE, prob = NULL)]
  # scRNA_S100.SeuObj <- scRNA.SeuObj[,sample(1:100,100, replace = FALSE, prob = NULL)]

  # ## Merging more than two seurat objects
  # ## Ref: https://github.com/satijalab/seurat/issues/706
  # ## Ref: https://mojaveazure.github.io/seurat-object/reference/Seurat-methods.html
  #
  # scRNA.SeuObj_Small_Merge <- merge(scRNA_S50.SeuObj, scRNA_S100.SeuObj) # merge(x,y,add.cell.ids = c(x@project.name,y@project.name))

  ## Extract SubSet
  # scRNA_Fib.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% c("Fibroblast cell")]


  CellType.set <- scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique()
  for (i in 1:length(CellType.set)) {
    if(i==1){
      scRNA.SeuObj_Small <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% CellType.set[i]]
      scRNA.SeuObj_Small <- scRNA.SeuObj_Small[,sample(1:nrow((scRNA.SeuObj_Small@meta.data)),50, replace = FALSE, prob = NULL)]
    }else{
      scRNA.SeuObj_Small_Temp <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% CellType.set[i]]
      scRNA.SeuObj_Small_Temp <- scRNA.SeuObj_Small_Temp[,sample(1:nrow((scRNA.SeuObj_Small_Temp@meta.data)),50, replace = FALSE, prob = NULL)]
      scRNA.SeuObj_Small <- merge(scRNA.SeuObj_Small,scRNA.SeuObj_Small_Temp)
    }
  }
  rm(i,scRNA.SeuObj_Small_Temp)



#### Save the RData ####
  rm(list=setdiff(ls(), c("scRNA.SeuObj_Small","Version","Save.Path")))
  save.image(paste0(Save.Path,"/",Version,"_Seurat-Small-Data-Set.RData"))
