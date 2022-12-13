##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

#### Load the required libraries ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","Seurat")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("SeuratDisk")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)

  ##### Current path and new folder setting*  #####
  ProjectName = "SeuratSmall"
  Version = paste0(Sys.Date(),"_",ProjectName,"_PADC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

#### Load data ####
  #### Converse h5ad to Seurat ####
  remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)

  # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
  Convert("StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad", "PRJCA001063.h5seurat")

  # This .d5seurat object can then be read in manually
  seuratObject <- LoadH5Seurat("PRJCA001063.h5seurat")
  seuratObject_Ori <- seuratObject
#### Data preprocessing ####
  seuratObject@meta.data[["orig.ident"]] <- sample(c(1,2),nrow(seuratObject@meta.data),replace = TRUE)
  seuratObject@meta.data <- seuratObject@meta.data %>%
                            relocate(orig.ident, .before = colnames(seuratObject@meta.data)[1])%>%
                            relocate(n_counts, .before = colnames(seuratObject@meta.data)[1]) %>%
                            relocate(n_genes, .before = colnames(seuratObject@meta.data)[1]) %>%


                            rename(nFeature_RNA = n_genes, nCount_RNA = n_counts)



#### Sampling ####
  # TTT <- sample(seuratObject@meta.data[["Cell_type"]] %>% unique(), 100, replace = TRUE, prob = NULL)
  # FibroblastTTT <- seuratObject[,seuratObject@meta.data[["Cell_type"]] %in% c("Fibroblast cell")]
  #
  # seuratObject_Small <- seuratObject[,sample(1:100,50, replace = FALSE, prob = NULL)]
  # seuratObject_Small2 <- seuratObject[,sample(1:100,50, replace = FALSE, prob = NULL)]
  #
  # ## Merging more than two seurat objects
  # ## Ref: https://github.com/satijalab/seurat/issues/706
  # ## Ref: https://mojaveazure.github.io/seurat-object/reference/Seurat-methods.html
  #
  # seuratObject_Small_Merge <- merge(seuratObject_Small,seuratObject_Small2) # merge(x,y,add.cell.ids = c(x@project.name,y@project.name))

  CellType.set <- seuratObject@meta.data[["Cell_type"]] %>% unique()
  for (i in 1:length(CellType.set)) {
    if(i==1){
      seuratObject_Small <- seuratObject[,seuratObject@meta.data[["Cell_type"]] %in% CellType.set[i]]
      seuratObject_Small <- seuratObject_Small[,sample(1:nrow((seuratObject_Small@meta.data)),50, replace = FALSE, prob = NULL)]
    }else{
      seuratObject_Small_Temp <- seuratObject[,seuratObject@meta.data[["Cell_type"]] %in% CellType.set[i]]
      seuratObject_Small_Temp <- seuratObject_Small_Temp[,sample(1:nrow((seuratObject_Small_Temp@meta.data)),50, replace = FALSE, prob = NULL)]
      seuratObject_Small <- merge(seuratObject_Small,seuratObject_Small_Temp)
    }
  }
  rm(i,seuratObject_Small_Temp)



#### Save the RData ####
  rm(list=setdiff(ls(), c("seuratObject_Small","Version","Save.Path")))
  save.image(paste0(Save.Path,"/",Version,"_Seurat-Small-Data-Set.RData"))
