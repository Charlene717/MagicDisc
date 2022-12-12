##### Load Packages #####
source("FUN_Package_InstLoad.R")
# FUN_Basic.set <- c("tidyverse","Seurat","monocle","ggplot2","ggpmisc","broom",
#                    "stringr","magrittr","dplyr",
#                    "patchwork","reticulate","anndata")
# FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb",
#                          "basilisk","zellkonverter","SeuratDisk",
#                          "SingleR","scRNAseq","celldex","scran")
# ## Set the desired organism
# # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
# # c(organism,"fgsea")

FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
rm(FUN_Basic.set, FUN_BiocManager.set)

#### GitHub installation ####
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("cole-trapnell-lab/garnett")
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("sqjin/CellChat")
# devtools::install_github("LTLA/SingleR")

library(monocle3)
library(garnett)
library(CellChat)
# library(SingleR)

##### Function setting #####
## Prepossession
source("FUN_ReadscRNA.R")
source("FUN_Cal_Mit.R")
source("FUN_scRNAQC.R")
source("FUN_HSsymbol2MMsymbol.R")
source("FUN_CombineSeuObj.R")
source("FUN_DRCluster.R")
source("FUN_UMAP_CellTypeMarker.R")
source("FUN_Export_All_DRPlot.R")

## Summarize
source("FUN_MetaSummary.R")
source("FUN_AnnoSummary.R")
source("FUN_Export_CellCount.R")

## Find Markers
source("FUN_Find_Markers.R")
source("FUN_VolcanoPlot.R")
source("FUN_Venn.R")
source("FUN_BioMarker1Index.R")
source("FUN_BioMarker2Index.R")

## Enrichment analysis
source("FUN_GSEA_LargeGeneSet.R")
source("FUN_GSEA_ggplot.R")
source("FUN_GSEA_MultiCell.R")

## Cell-Cell interaction
source("FUN_CellChatOne.R")

## inferCNV
source("FUN_inferCNV.R")

## Beautify Plot
source("FUN_Beautify_ggplot.R")
source("FUN_Beautify_UMAP.R")
source("FUN_Beautify_Heatmap_Seurat.R")
source("FUN_BeautifyVennDiag.R")
source("FUN_BeautifyDotPlot.R")

## Verification
source("Fun_Draw_ConfuMax.R")
