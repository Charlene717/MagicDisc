##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)


  #
  # i=160
  # j=0.3
  # k=25
  # seuratObject@reductions[["umap"]]@cell.embeddings <- seuratObject@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]]
  #
  #
  # FeaturePlot(seuratObject, features = c("TOP2A"))
  # FeaturePlot(seuratObject, features = c("ALDH1A1","ALDH1A2","ALDH1A3","NR5A2","KRT19","SOX9"))
  #
  #
  #
  #
  # FeaturePlot(seuratObject, features = c("HNF1B","HNF6","FOXA2","HNF4A","HEX","GATA4","GATA6"))
  # FeaturePlot(seuratObject, features = c("MNX1","PTF1A","PDX1","HNF1B","HNF6","FOXA2","HNF4",
  #                                        "HEX","GATA4","GATA6"))
  #
  # FeaturePlot(seuratObject, features = c("SOX9","NKX6-1","NKX2-2","MNX1","PTF1A","PDX1",
  #                                        "HNF1B","HNF6","FOXA2","HNF4","HEX","GATA4","GATA6"))
  #
  # FeaturePlot(seuratObject, features = c("HES1","PROX1","SOX9","NKX6-1","NKX2-2","PTF1A","PDX1",
  #                                        "HNF1B","HNF6","FOXA2","HNF4A","GATA4","GATA6"))
  #
  # FeaturePlot(seuratObject, features = c("PTF1A","CPA1","MYC","NR5A2","HNF1B","MNX1","SOX9"))
  # FeaturePlot(seuratObject, features = c("HNF1B","HNF6","FOXA2","NR5A2","GLIS3","PROX1"))








##### Data prepocessing #####
  load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  seuratObject_1 <- scRNA.SeuObj
  seuratObject_1@meta.data[["DataSetID"]] <- rep("PRJCA001063")
  rm(list=setdiff(ls(), "seuratObject_1"))

  load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-25_PDAC_GSE131886_SC/04_Perform_an_integrated_analysis.RData")
  seuratObject_2 <- scRNA.SeuObj
  DefaultAssay(seuratObject_2) <- "RNA"
  rm(list=setdiff(ls(), str_subset(objects(), pattern = "seuratObject")))


  load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-23_PDAC_GSE154778_SC/04_Perform_an_integrated_analysis.RData")
  seuratObject_3 <- scRNA.SeuObj
  DefaultAssay(seuratObject_3) <- "RNA"

  scRNA_SeuObj.list <- list(PRJCA001063 = seuratObject_1,
                            GSE131886 = seuratObject_2,
                            GSE154778 = seuratObject_3)

  rm(list=setdiff(ls(), "scRNA_SeuObj.list"))


##### Current path and new folder setting* #####
  ProjectName = "TrajAna_PCA"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }


##### Load Packages #####
  #### Basic installation ####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                   "stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("fgsea","AnnotationHub","ensembldb",
                   "SeuratDisk","monocle",
                   "SingleR","scRNAseq","celldex","scran")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)

  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  library(monocle)
  devtools::install_github("cole-trapnell-lab/garnett")
  devtools::install_github('cole-trapnell-lab/monocle3')
  devtools::install_github("LTLA/SingleR")

  library(monocle3)
  library(garnett)
  # library(SingleR)




##### Function setting #####
  ## Call function
  source("FUN_Cal_Mit.R")
  source("FUN_CombineSeuObj.R")
  source("FUN_Beautify_ggplot.R")


##### CombineSeuObj #####
  scRNA.SeuObj <- CombineSeuObj(scRNA_SeuObj.list)
  rm(scRNA_SeuObj.list)

  # Run the standard workflow for visualization and clustering
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 1000, verbose = FALSE)
  # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:1000)
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Ori.RData"))


  ##### Plot #####
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap")
  DimPlot(scRNA.SeuObj, reduction = "umap",label = T)
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")

##### Export figures #####
## Export TIFF
for (i in seq(80,400,80)) {
  for (j in seq(0.1,0.9,0.2)) {
    for (k in seq(20,800,60)) {
      try({
        set.seed(1)
        scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:i,n.neighbors = k, min.dist= j)
        scRNA.SeuObj@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
        # scRNA.SeuObj@reductions[["umap"]]@cell.embeddings <- scRNA.SeuObj@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]]

        Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
        p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
          ggtitle(paste0("CellType","  PCA:",i,"  NNe:",k,"  MD:",j)) +
          theme(plot.title = element_text(hjust = 0.5,vjust = 0))
        tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_CellType.tiff"),
             width = 28, height = 20, units = "cm", res = 200)
        print(p)
        graphics.off()

        Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster"]]
        p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
          ggtitle(paste0("ReCluster","  PCA:",i,"  NNe:",k,"  MD:",j)) +
          theme(plot.title = element_text(hjust = 0.5,vjust = 0))
        tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_ReCluster.tiff"),
             width = 35, height = 20, units = "cm", res = 200)
        print(p)
        graphics.off()

        Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["DataSetID"]]
        p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
          ggtitle(paste0("ReCluster","  PCA:",i,"  NNe:",k,"  MD:",j)) +
          theme(plot.title = element_text(hjust = 0.5,vjust = 0))
        tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_DataSetID.tiff"),
             width = 35, height = 20, units = "cm", res = 200)
        print(p)
        graphics.off()

        p <-  FeaturePlot(scRNA.SeuObj, features = c("TOP2A")) %>% BeautifyggPlot(LegPos = c(1.02, 0.15)) +
          ggtitle(paste0("TOP2A","  PCA:",i,"  NNe:",k,"  MD:",j)) +
          theme(plot.title = element_text(hjust = 0.5,vjust = 0))
        tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_TOP2A.tiff"),
             width = 28, height = 20, units = "cm", res = 200)
        print(p)
        graphics.off()
      })
    }
  }
}



# scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:100,n.neighbors = 20, min.dist=0.05)
# # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:100,n.neighbors = 1000, min.dist=0.1)
# # scRNA.SeuObj@meta.data[["UMAP_NNe1000"]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
# # scRNA.SeuObj@meta.data <- scRNA.SeuObj@meta.data[,!colnames(scRNA.SeuObj@meta.data)=="UMAP_NNe1000"]
# scRNA.SeuObj@meta.data[["UMAP_NNe20_MD03"]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
#
# DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster"]]
DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))

  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_TryCondition.RData"))




###########################################################################################
library(tidyverse)
Genes <- c("ALDH1A1","ALDH1A2","ALDH1A3","NR5A2","KRT19","SOX9",
           "HNF1B","HNF6","FOXA2","HNF4A","HEX","GATA4","GATA6",
           "MNX1","PTF1A","PDX1","HNF1B","HNF6","FOXA2","HNF4", "HEX","GATA4","GATA6",
           "PTF1A","CPA1","MYC","NR5A2","HNF1B","MNX1","SOX9",

           "SOX9","NKX6-1","NKX2-2","MNX1","PTF1A","PDX1",
           "HNF1B","HNF6","FOXA2","HNF4","HEX","GATA4","GATA6",
           "HES1","PROX1","SOX9","NKX6-1","NKX2-2","PTF1A","PDX1","HNF1B","HNF6",
           "FOXA2","HNF4A","GATA4","GATA6") %>% unique()

library(Seurat)
DoHeatmap(scRNA.SeuObj, features = Genes) + NoLegend()

HM <- DoHeatmap(scRNA.SeuObj, features = Genes, # group.by = Type,
                size = 4, angle = 90) +
  scale_fill_gradient2(low = "#5283ff",
                       mid = "white",
                       high = "#ff5c5c") +
  theme(axis.text.y = element_text(size  = 20)) +
  theme(legend.position = "bottom" )
print(HM)

## Ref: https://github.com/satijalab/seurat/issues/1369
scRNA.SeuObj <- ScaleData(object = scRNA.SeuObj, features = rownames(scRNA.SeuObj))
