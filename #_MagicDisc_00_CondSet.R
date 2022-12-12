##### Setting ######
##### Load Packages #####
FUN_Basic.set <- c("tidyverse","Seurat","monocle","ggplot2","ggpmisc","broom",
                   "stringr","magrittr","dplyr", "patchwork","reticulate","anndata")
FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb",
                         "basilisk","zellkonverter","SeuratDisk",
                         "SingleR","scRNAseq","celldex","scran")

#### Current path and new folder setting* ####
ProjectName = "CC"
Sampletype = "PBMC"
Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
Save.Path = paste0(getwd(),"/",Version)

## Create new folder
if (!dir.exists(Save.Path)){dir.create(Save.Path)}

##### Import information sheet setting* #####
## Folder name
InputFolder = "Input_files_10x"

## Metadata spreadsheet
InputAnno = "PBMC_Ano.csv"

## GSEA Genesets file
InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"

## InferCNV Ref file
InputInferCNV = "mm10_genomic_mapinfo_one.tsv"

##### Parameter setting* #####
ClassSet1 = "Sample"
ClassSet2 = "Cachexia"
ClassSet3 = "Sex"

DataMode = "10x"
Species = "Mouse" # Species = c("Mouse","Human")
