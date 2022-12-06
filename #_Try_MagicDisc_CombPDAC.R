##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)
# options(stringsAsFactors = FALSE)
# Sys.setlocale(category = "LC_ALL", locale = "UTF-8")

##### Load datasets  #####
## Load all
load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-05_CTAnno_singleR_RefPRJCA001063_PDAC.RData")
## Clean up the object
rm(list=setdiff(ls(), c("scRNA.SeuObj")))

scRNA.SeuObj$celltype <- scRNA.SeuObj$singleR_classic_PredbyscRNA

##### Setting ######
  #### Current path and new folder setting* ####
  ProjectName = "Combine"
  Sampletype = "PDAC"
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

  # ## InferCNV Ref file
  # InputInferCNV = "mm10_genomic_mapinfo_one.tsv"

  ##### Parameter setting* #####
  ClassSet1 = "DataSetID"
  ClassSet2 = "Type"
  ClassSet3 = "Sex"

  DataMode = "10x"
  Species = "Human" # Species = c("Mouse","Human")

  ## Extract the original Meta term
  Ori_Meta.set <- colnames(scRNA.SeuObj@meta.data)

# ##### Export the log file (Start) #####
#   ## Create new folder for log file
#   if (!dir.exists("LogFiles")){dir.create("LogFiles")}
#
#   my_log <- file(paste0("LogFiles/",Version,"_log.txt")) # File name of output log
#
#   sink(my_log, append = TRUE, type = "output") # Writing console output to log file
#   sink(my_log, append = TRUE, type = "message")
#
#   cat(readChar(rstudioapi::getSourceEditorContext()$path, # Writing currently opened R script to file
#                file.info(rstudioapi::getSourceEditorContext()$path)$size))


##### Load Packages #####
source("FUN_Package_InstLoad.R")
FUN_Basic.set <- c("tidyverse","Seurat","monocle","ggplot2","ggpmisc","broom",
                   "stringr","magrittr","dplyr",
                   "patchwork","reticulate","anndata")
FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb",
                         "basilisk","zellkonverter","SeuratDisk",
                         "SingleR","scRNAseq","celldex","scran")
## Set the desired organism
# organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
# c(organism,"fgsea")

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


##### Load datasets  #####
  # ## Annotation table
  # list_files.df <- read.csv(paste0(InputFolder,"/",InputAnno))
  # Feature.set <- colnames(list_files.df)[-1]
  #
  # ## Read 10x files
  # source("FUN_ReadscRNA.R")
  # scRNA_SeuObj.list <- FUN_ReadscRNA(InputFolder,Folder = Folder,
  #                                Path =  "/monocle/outs/filtered_gene_bc_matrices/mm10",
  #                                list_files.df, Mode = DataMode, projectName = ProjectName)

##### 07 Count Cell number  #####
  source("FUN_AnnoSummary.R")
  source("FUN_Export_CellCount.R")

  ## Create new folder
  PathCellCount <- paste0(Save.Path,"/","B01_CellCount")
  if (!dir.exists(PathCellCount)){
    dir.create(PathCellCount)
  }

  ## Annotation Summary Table
  AnnoSummary.lt <- AnnoSummary(scRNA.SeuObj,  list_files.df, Ori_Meta.set,
                                ClassSet = ClassSet1, ClassSet2 = ClassSet2)

  ## ExportCellCount
  ExportCellCount(AnnoSummary.lt,Path = PathCellCount)

  #### Save RData ####
  save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))

##### 08 Find CCmarker in different Cell type ########
  ## Create new folder
  PathBiomarkers <- paste0(Save.Path,"/","B02_Biomarkers")
  if (!dir.exists(PathBiomarkers)){
    dir.create(PathBiomarkers)
  }

  #### Define group by different phenotype ####
  source("FUN_Find_Markers.R")
  source("FUN_BioMarker2Index.R")


  # scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)

  ### Define group by different phenotype ###
  Type = "celltype"
  scRNA.SeuObj[[paste0(Type,".",ClassSet2)]] <- paste(Idents(scRNA.SeuObj),
                                                      as.matrix(scRNA.SeuObj[[ClassSet2]]), sep = "_")
  scRNA.SeuObj[[paste0(Type,".",ClassSet2,".",ClassSet3)]] <- paste(Idents(scRNA.SeuObj),
                                                                    as.matrix(scRNA.SeuObj[[ClassSet2]]),
                                                                    as.matrix(scRNA.SeuObj[[ClassSet3]]), sep = "_")

  CellType.set <- scRNA.SeuObj@meta.data[["celltype"]] %>% unique()
  ##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  Idents(scRNA.SeuObj) <- paste0(Type,".",ClassSet2)

  DefaultAssay(scRNA.SeuObj) <- "RNA"
  CCMarker.lt <- BioMarker1Index(scRNA.SeuObj, Path = PathBiomarkers, projectName = ProjectName,
                                 sampletype = Sampletype, cellType.set = CellType.set, classSet2 = ClassSet2,
                                 Type = paste0("celltype.",ClassSet2) )

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_1_Find__",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Pooled).RData"))

  if(ClassSet3 != ""){

##### 08_2 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  Idents(scRNA.SeuObj) <- paste0(Type,".",ClassSet2,".",ClassSet3)
  DefaultAssay(scRNA.SeuObj) <- "RNA"

  ## Find CCmarker in different Cell type
  CCMarker2Index.lt <-  BioMarker2Index(scRNA.SeuObj, Path = PathBiomarkers, projectName = ProjectName,
                                        sampletype = Sampletype, cellType.set = CellType.set, classSet2 = ClassSet2,
                                        classSet3 = ClassSet3,Type = "celltype")

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_2_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Separate).RData"))


##### 08_3 Find CCmarker in different Cell type and VennDiagrame (SSA_IntersectCT) ########
  ##-------------- Intersect_CellType --------------##
  CCMarker_Male.lt <- CCMarker2Index.lt[[1]]
  CCMarker_Female.lt <- CCMarker2Index.lt[[2]]
  intersect_CellType <- intersect(names(CCMarker_Male.lt),names(CCMarker_Female.lt))

  CCMarker_Male.lt <- CCMarker_Male.lt[names(CCMarker_Male.lt) %in% intersect_CellType]
  CCMarker_Female.lt <- CCMarker_Female.lt[names(CCMarker_Female.lt) %in% intersect_CellType]

  CellType.list <- names(CCMarker_Male.lt)

  Venn_CCMarke.lt <- BeautifyVennDiag(CCMarker_Male.lt, CCMarker_Female.lt, CellType.list,list_files.df,
                                      classSet3 = ClassSet3)
  }
  #### Save RData ####
  save.image(paste0(Save.Path,"/08_3_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_Venn.RData"))



################## (Pending) CCmarker matrix (Heatmap) ##################
################## (Pending) CCmarker matrix LogFC (Heatmap) ##################

    #####------------------------------------------------------------------------------------------------------------#####

##### 09_0 GSEA Analysis (Geneset Prepare) #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_GSEA_ggplot.R")

  ## Load the GSEA Dataset
  load("GSEA_Analysis_Geneset.RData")

  # # Geneset from GSEA
  # # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  # Pathway.all <- read.delim2(paste0(getwd(),"/",InputGSEA),
  #                            col.names = 1:max(count.fields(paste0(getwd(),"/",InputGSEA))),
  #                            header = F,sep = "\t")


##### 09_1 GSEA Analysis (SPA) #####
  ## Create folder
  GSEA_SPA.lt <- FUN_GSEA_MultiCell(CCMarker.lt, CellType.list,Path = Save.Path, sampletype = Sampletype,
                                    projectName = ProjectName,NES_TH = 1.5, Padj_TH = 0.01)

  ##### save.image #####
  save.image(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

  if(ClassSet3 != ""){
##### 09_2 GSEA Analysis (SSA_MAle) #####
  ClassSet3.set <- list_files.df[,ClassSet3] %>% unique()
  GSEA_SSA_Male.lt <- FUN_GSEA_MultiCell(CCMarker_Male.lt, CellType.list,Path = Save.Path, sampletype = Sampletype,
                                    projectName = ProjectName, Type = ClassSet3.set[1],NES_TH = 1.5, Padj_TH = 0.01)

  ##### save.image #####
  save.image(paste0(Save.Path,"/09_2_GSEA_Analysis_(SSA_Male).RData"))

##### 09_2 GSEA Analysis (SSA_Female) #####
  GSEA_SSA_Female.lt <- FUN_GSEA_MultiCell(CCMarker_Female.lt, CellType.list,Path = Save.Path, sampletype = Sampletype,
                                         projectName = ProjectName, Type = ClassSet3.set[2],NES_TH = 1.5, Padj_TH = 0.01)


  ##### save.image #####
  save.image(paste0(Save.Path,"/09_2_GSEA_Analysis_(SSA_Female).RData"))
  }

##### 010 Cell-cell interaction #####
  ## ECM-Receptor
  CellChatOne(scRNA.SeuObj,
              signalingtype = "ECM-Receptor", projectName = "ECM",
              save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
              groupby = "celltype",species = Species
              ) ->   CellChat_ECM.lt

  ## Cell-Cell Contact
  CellChatOne(scRNA.SeuObj,
              signalingtype = "Cell-Cell Contact", projectName = "CC",
              save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
              groupby = "celltype",species = Species
              ) -> CellChat_CC.lt

  ## Secreted Signaling
  CellChatOne(scRNA.SeuObj,
              signalingtype = "Secreted Signaling", projectName = "Secret",
              save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
              groupby = "celltype",species = Species
              ) -> CellChat_Secret.lt

  ##### save.image #####
  save.image(paste0(Save.Path,"/010_Cell_Cell_Interaction.RData"))

##### GO/Metascape #####

##### Deconvolution #####

##### Clinical analysis #####


# ##### 013 inferCNV #####
#   ## Create new folder
#   PathinferCNV <- paste0(Save.Path,"/","D01_inferCNV")
#   if (!dir.exists(PathinferCNV)){
#     dir.create(PathinferCNV)
#   }
#
#   infercnv_obj <- inferCNV(scRNA.SeuObj, AnnoSet = "celltype",
#                            SpeciSet = Species,
#                            Path = PathinferCNV,
#                            RefSet = c("T","B"),
#                            CreateInfercnvObject.lt = list(chr_exclude = c("chrM")))
#   ##### save.image #####
#   save.image(paste0(Save.Path,"/013_inferCNV.RData"))


##### Beautify Figs #####
## Test

# ##### Export the log file (End) #####
#   closeAllConnections() # Close connection to log file

