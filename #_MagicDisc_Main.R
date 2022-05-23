## Ref: https://statisticsglobe.com/r-save-all-console-input-output-to-file
## Ref: https://blog.gtwang.org/r/r-data-input-and-output/

##### Export the log file (Start) #####
  my_log <- file("MagicDisc_log.txt") # File name of output log

  sink(my_log, append = TRUE, type = "output") # Writing console output to log file
  sink(my_log, append = TRUE, type = "message")

  cat(readChar(rstudioapi::getSourceEditorContext()$path, # Writing currently opened R script to file
               file.info(rstudioapi::getSourceEditorContext()$path)$size))


##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc","broom",
                   "stringr","magrittr","dplyr",
                   "CellChat","patchwork","reticulate","anndata")
  ## Check whether the installation of those packages is required from basic
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
  Package.set <- c("fgsea","AnnotationHub","ensembldb",
                   "basilisk","zellkonverter","SeuratDisk")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)


##### Function setting #####
  ## Call function
  source("FUN_ReadscRNA.R")
  source("FUN_Cal_Mit.R")
  source("FUN_scRNAQC.R")
  source("FUN_MetaSummary.R")
  source("FUN_AnnoSummary.R")

  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_Venn.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Beautify_UMAP.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")

  source("Fun_Draw_ConfuMax.R")

  source("FUN_CombineSeuObj.R")
  source("FUN_DRCluster.R")
  source("FUN_Export_CellCount.R")
  source("FUN_Beautify_Heatmap_Seurat.R")
  source("FUN_UMAP_CellTypeMarker.R")
  source("FUN_Export_All_DRPlot.R")
  source("FUN_BeautifyDotPlot.R")
  source("FUN_BioMarker2Index.R")
  source("FUN_BeautifyVennDiag.R")
  source("FUN_BioMarker1Index.R")
  source("FUN_CellChatOne.R")
  source("FUN_GSEA_MultiCell.R")

##### Current path and new folder setting* #####
  ProjectName = "CC"
  Sampletype = "PBMC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_","CC_PBMC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }


  ## Import information
  InputFolder = "Input_files_10x"
  InputAnno = "PBMC_Ano.csv"

  InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"

##### Parameter setting* #####
  ClassSet1 = "Sample"
  ClassSet2 = "Cachexia"
  ClassSet3 = "Sex"
  DataMode="10x"

##### Load datasets  #####
  ## Annotation table
  list_files.df <- read.csv(paste0(InputFolder,"/",InputAnno))
  Feature.set <- colnames(list_files.df)[-1]

  ## Read 10x files
  scRNA_SeuObj.list <- ReadscRNA(InputFolder,Folder = Folder,
                                 Path =  "/monocle/outs/filtered_gene_bc_matrices/mm10",
                                 list_files.df, Mode = DataMode)

##### 01 Combine different datasets before QC #####
  source("FUN_CombineSeuObj.R")
  ## Combine SeuObjs from list before QC
  # (About 30 min for 20000 cells)
  scRNA.SeuObj <- CombineSeuObj(scRNA_SeuObj.list)

  ## Extract the original Meta term
  Ori_Meta.set <- colnames(scRNA.SeuObj@meta.data)

  #### Save RData ####
  save.image(paste0(Save.Path,"/01_Combine_different_datasets_before_QC.RData"))


##### 02 Quality Control #####
  source("FUN_scRNAQC.R")
  ## Create new folder
  PathQC <- paste0(Save.Path,"/","A01_QC")
  if (!dir.exists(PathQC)){
    dir.create(PathQC)
  }

  ## QC for all samples
  scRNA.SeuObj_QCTry <- scRNAQC(scRNA.SeuObj, Path = PathQC ,FileName = paste0(ProjectName,"_QCTry"))

  ## QC for each sample for the new integration
  #Test# scRNA_Ori.SeuObj.list <- SplitObject(scRNA_Ori.SeuObj, split.by = "ID")
    scRNA_SeuObj_QC.list <- list()
  for (i in 1:length(scRNA_SeuObj.list)) {

    Name <- names(scRNA_SeuObj.list)[[i]]
    scRNA_SeuObj_QC.list[[i]] <- scRNAQC(scRNA_SeuObj.list[[i]], Path = PathQC ,
                                         FileName = paste0(ProjectName,"_", Name,"_QC"))
    names(scRNA_SeuObj_QC.list)[[i]] <- Name

    }
  rm(i,Name)

  scRNA_Ori.SeuObj <- scRNA.SeuObj # Save the original obj
  rm(scRNA.SeuObj,scRNA.SeuObj_QCTry)

  #### Save RData ####
  save.image(paste0(Save.Path,"/02_Quality_Control.RData"))

##### 03 Combine different data sets after QC #####
  ## Combine SeuObjs from list after QC
  # (About 30 min for 20000 cells)
  scRNA.SeuObj <- CombineSeuObj(scRNA_SeuObj_QC.list)

  ## Check QC
  scRNAQC(scRNA.SeuObj,AddMitInf = "No",CheckOnly="Yes", Path = PathQC ,FileName = paste0(ProjectName,"_QCTry"))

  #### Save RData ####
  save.image(paste0(Save.Path,"/03_Combine_different_data_sets_after_QC.RData"))

##### 04 Perform an integrated analysis #####
  source("FUN_DRCluster.R")
  ## Create new folder
  PathCluster <- paste0(Save.Path,"/","A02_Cluster")
  if (!dir.exists(PathCluster)){
    dir.create(PathCluster)
  }

  scRNA.SeuObj <- DRCluster(scRNA_SeuObj.list, seed=1, PCAdims = 30,
                            Path = PathCluster, projectName= ProjectName,
                            MetaSet = Ori_Meta.set)

  ##### Meta Table  #####
  Meta.df <- MetaSummary(scRNA_SeuObj.list, scRNA.SeuObj,
                         scRNA_SeuObj_QC.list,scRNA_Ori.SeuObj)

  #### Save RData ####
  save.image(paste0(Save.Path,"/04_Perform_an_integrated_analysis.RData"))

################## (Pending) Cell Cycle Regression ##################


################## (Pending) Auto Cell type annotation ##################



##### 05 Identify conserved cluster markers  #####
  ## Create new folder
  PathCellType <- paste0(Save.Path,"/","A03_CellTypeAno")
  if (!dir.exists(PathCellType)){
    dir.create(PathCellType)
  }

  ## Identify conserved cluster markers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  set.seed(1) # Fix the seed
  PBMC.markers <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  scRNA.SeuObj <- Beautify_Heatmap_Seurat(scRNA.SeuObj, PBMC.markers, topN = 7, Path = PathCluster,
                                          Type = "seurat_clusters",
                                          projectName = ProjectName)


  # --------------- Check specific tissue marker --------------- #
  #(Pending)
  UMAP_CellTypeMarker()

  #### Save RData ####
  save.image(paste0(Save.Path,"/05_Identify_conserved_cell_type_markers.RData"))

##### 06 Cell type annotation  #####
  # scRNA.SeuObj.copy <- scRNA.SeuObj

  ## CD4+T: CD4+T Cell; CD8+T: CD8+T Cell; T: T Cell; B: B Cell; Mac: Macrophages;
  ## Neu: Neutrophils; NK: NK Cell; Mast: Mast Cell; Ery: Erythrocytes;
  ## Thr: Thrombocytes
  scRNA.SeuObj <- RenameIdents(scRNA.SeuObj, `0` = "CD4+T", `1` = "B", `2` = "Mac3",
                                `3` = "Neu", `4` = "CD8+T", `5` = "CD8+T", `6` = "Mac2", `7` = "CD4+T",
                                `8` = "NK", `9` = "Neu",`10` = "Mast1", `11` = "T", `12` = "Ery", `13` = "Mac1",
                                `14` = "B", `15` = "B", `16` = "Mast2", `17` = "Mac0", `18` = "Neu")

  Cell_Type_Order.set <- c("T", "CD4+T", "CD8+T", "B" , "Mac0", "Mac1", "Mac2", "Mac3", "Mast1", "Mast2", "NK", "Neu", "Ery")

  scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
  # Idents(scRNA.SeuObj) <- "celltype"

  ## Heatmap
  Heatmap_Color.lt <- list(low="#5283ff",mid ="white", high ="#ff5c5c")
  scRNA.SeuObj <- Beautify_Heatmap_Seurat(scRNA.SeuObj, PBMC.markers, topN = 7, Path = PathCellType,
                                          Type = "celltype", HMColor.lt = Heatmap_Color.lt,
                                          projectName = ProjectName)

  ## Export All DRPlot(UMAP,tSNE)
  Export_All_DRPlot(scRNA.SeuObj)

  ## Summary
  markers.to.plot <- c("Cd3d","Cd3e", "Cd4","Cd8a", "Csf1r", "Lyz2","Chil3","Il1b", "S100a9","Nkg7",
                       "Gzmb", "Cd79a", "Ms4a1","Clu","Hbb-bs","Ppbp")
  ## DotPlot
  #(Pending)
  BeautifyDotPlot(scRNA.SeuObj, Path = PathCellType, projectName = ProjectName,
                  Features = markers.to.plot)


  #### Save RData ####
  save.image(paste0(Save.Path,"/06_Cell_type_annotation.RData"))

##### 07 Count Cell number  #####
  source("FUN_AnnoSummary.R")
  source("FUN_Export_CellCount.R")

  ## Annotation Summary Table
  AnnoSummary.lt <- AnnoSummary(scRNA.SeuObj,  list_files.df, Ori_Meta.set,
                                ClassSet = ClassSet1, ClassSet2 = ClassSet2)

  ## ExportCellCount
  ExportCellCount(AnnoSummary.lt)

  #### Save RData ####
  save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))

##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  #### Define group by different phenotype ####
  source("FUN_Find_Markers.R")
  source("FUN_BioMarker2Index.R")

  ## Create new folder
  PathBiomarkers <- paste0(Save.Path,"/","B02_Biomarkers")
  if (!dir.exists(PathBiomarkers)){
    dir.create(PathBiomarkers)
  }

  # scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)

  ### Define group by different phenotype ###
  Type = "celltype"
  scRNA.SeuObj[[paste0(Type,".",ClassSet2)]] <- paste(Idents(scRNA.SeuObj),
                                                      as.matrix(scRNA.SeuObj[[ClassSet2]]), sep = "_")
  scRNA.SeuObj[[paste0(Type,".",ClassSet2,".",ClassSet3)]] <- paste(Idents(scRNA.SeuObj),
                                                                    as.matrix(scRNA.SeuObj[[ClassSet2]]),
                                                                    as.matrix(scRNA.SeuObj[[ClassSet3]]), sep = "_")
  Idents(scRNA.SeuObj) <- paste0(Type,".",ClassSet2,".",ClassSet3)
  DefaultAssay(scRNA.SeuObj) <- "RNA"

  ## Find CCmarker in different Cell type
  CCMarker2Index.lt <-  BioMarker2Index(scRNA.SeuObj, Path = PathBiomarkers, projectName = ProjectName,
                                  classSet2 = ClassSet2, classSet3 = ClassSet3,Type = "celltype")

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_1_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Separate).RData"))


##### 08_2 Find CCmarker in different Cell type and VennDiagrame (SSA_IntersectCT) ########
  ##-------------- Intersect_CellType --------------##
  CCMarker_Male.lt <- CCMarker2Index.lt[[1]]
  CCMarker_Female.lt <- CCMarker2Index.lt[[2]]
  intersect_CellType <- intersect(names(CCMarker_Male.lt),names(CCMarker_Female.lt))

  CCMarker_Male.lt <- CCMarker_Male.lt[names(CCMarker_Male.lt) %in% intersect_CellType]
  CCMarker_Female.lt <- CCMarker_Female.lt[names(CCMarker_Female.lt) %in% intersect_CellType]

  CellType.list <- names(CCMarker_Male.lt)

  Venn_CCMarke.lt <- BeautifyVennDiag(CCMarker_Male.lt, CCMarker_Female.lt, CellType.list,list_files.df,
                                      classSet3 = ClassSet3)

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_2_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_Venn.RData"))


##### 08_3 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")

  Idents(scRNA.SeuObj) <- paste0(Type,".",ClassSet2)
  DefaultAssay(scRNA.SeuObj) <- "RNA"
  CCMarker.lt <- BioMarker1Index(scRNA.SeuObj, Path = PathBiomarkers, projectName = ProjectName,
                                 sampletype = Sampletype, cellType.list = CellType.list, classSet2 = ClassSet2,
                                 Type = paste0("celltype.",ClassSet2) )

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_3_Find__",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Pooled).RData"))

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

##### 09_2 GSEA Analysis (SSA_MAle) #####
  GSEA_SSA_Male.lt <- FUN_GSEA_MultiCell(CCMarker_Male.lt, CellType.list,Path = Save.Path, sampletype = Sampletype,
                                    projectName = ProjectName,NES_TH = 1.5, Padj_TH = 0.01)

  ##### save.image #####
  save.image(paste0(Save.Path,"/09_2_GSEA_Analysis_(SSA_Male).RData"))

##### 09_3 GSEA Analysis (SSA_MAle) #####
  GSEA_SSA_Female.lt <- FUN_GSEA_MultiCell(CCMarker_Female.lt, CellType.list,Path = Save.Path, sampletype = Sampletype,
                                         projectName = ProjectName,NES_TH = 1.5, Padj_TH = 0.01)

  ##### save.image #####
  save.image(paste0(Save.Path,"/09_2_GSEA_Analysis_(SSA_Female).RData"))

##### Cell-cell interaction #####
  ## ECM-Receptor
  CellChatOne(scRNA.SeuObj,
              signalingtype = "ECM-Receptor", projectName = "ECM",
              save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
              groupby = "celltype",species = "mouse"
              ) ->   CellChat_ECM.lt

  ## Cell-Cell Contact
  CellChatOne(scRNA.SeuObj,
              signalingtype = "Cell-Cell Contact", projectName = "CC",
              save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
              groupby = "celltype",species = "mouse"
              ) -> CellChat_CC.lt

  ## Secreted Signaling
  CellChatOne(scRNA.SeuObj,
              signalingtype = "Secreted Signaling", projectName = "Secret",
              save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
              groupby = "celltype",species = "mouse"
  ) -> CellChat_Secret.lt


##### GO/Metascape #####



##### inferCNV #####

##### Deconvolution #####

##### Clinical analysis #####

##### Beautify Figs #####


##### Export the log file (End) #####
  closeAllConnections() # Close connection to log file

