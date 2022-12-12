## Ref: https://satijalab.org/seurat/
## Ref: https://cole-trapnell-lab.github.io/monocle3/

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)
  # options(stringsAsFactors = FALSE)
  # Sys.setlocale(category = "LC_ALL", locale = "UTF-8")

##### Setting ######
  source("#_MagicDisc_00_CondSet.R")

##**************************** Export the log file (Start) ****************************##
  ## Ref: https://statisticsglobe.com/r-save-all-console-input-output-to-file
  ## Ref: https://blog.gtwang.org/r/r-data-input-and-output/
  ## Create new folder for log file
  if (!dir.exists("LogFiles")){dir.create("LogFiles")}

  my_log <- file(paste0("LogFiles/",Version,"_log.txt")) # File name of output log

  sink(my_log, append = TRUE, type = "output") # Writing console output to log file
  sink(my_log, append = TRUE, type = "message")

  cat(readChar(rstudioapi::getSourceEditorContext()$path, # Writing currently opened R script to file
               file.info(rstudioapi::getSourceEditorContext()$path)$size))
##**************************** Export the log file (Start) ****************************##

##### Load Packages #####
  FUN_Basic.set <- c("tidyverse","Seurat","monocle","ggplot2","ggpmisc","broom",
                     "stringr","magrittr","dplyr", "patchwork","reticulate","anndata")
  FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb",
                           "basilisk","zellkonverter","SeuratDisk",
                           "SingleR","scRNAseq","celldex","scran")
  source("#_MagicDisc_00_PKG_FUN.R")

##### Load datasets  #####
  ## Annotation table
  list_files.df <- read.csv(paste0(InputFolder,"/",InputAnno))
  Feature.set <- colnames(list_files.df)[-1]

  ## Read 10x files
  source("FUN_ReadscRNA.R")
  scRNA_SeuObj.list <- FUN_ReadscRNA(InputFolder, Folder = Folder,
                                     Path =  "/monocle/outs/filtered_gene_bc_matrices/mm10", ## Hm: hg19 ; Mm: mm10
                                     list_files.df, Mode = DataMode, projectName = ProjectName)


  ##****************************************************************************##
  ################## (Pending) Cell Cycle Regression ##################
  ##****************************************************************************##

##### 01 Quality Control #####

##### 01 Combine different datasets before QC #####
  source("FUN_Cal_Mit.R")
  source("FUN_CombineSeuObj.R")

  ## Combine SeuObjs from list before QC   # (About 30 min for 20000 cells)
  scRNA.SeuObj <- CombineSeuObj(scRNA_SeuObj.list)

  ## Extract the original Meta term
  Ori_Meta.set <- colnames(scRNA.SeuObj@meta.data)

  #### Save RData ####
  save.image(paste0(Save.Path,"/01_Combine_different_datasets_before_QC.RData"))


##### 02 Quality Control #####
  source("FUN_scRNAQC.R")
  ## Create new folder
  PathQC <- paste0(Save.Path,"/","A01_QC")
  if (!dir.exists(PathQC)){dir.create(PathQC)}

  ## QC for all samples
  scRNA.SeuObj_QCTry <- scRNAQC(scRNA.SeuObj, Path = PathQC ,SpeciSet = Species,
                                FileName = paste0(ProjectName))

  ## QC for each sample for the new integration
  #Test# scRNA_Ori.SeuObj.list <- SplitObject(scRNA_Ori.SeuObj, split.by = "ID")
  scRNA_SeuObj_QC.list <- list()
  for (i in 1:length(scRNA_SeuObj.list)) {

    Name <- names(scRNA_SeuObj.list)[[i]]
    scRNA_SeuObj_QC.list[[i]] <- scRNAQC(scRNA_SeuObj.list[[i]], Path = PathQC ,SpeciSet = Species,
                                         FileName = paste0(ProjectName,"_",Name))
    names(scRNA_SeuObj_QC.list)[[i]] <- Name

    }
  rm(i,Name)

  scRNA_Ori.SeuObj <- scRNA.SeuObj # Save the original obj
  rm(scRNA.SeuObj,scRNA.SeuObj_QCTry)

  #### Save RData ####
  save.image(paste0(Save.Path,"/02_Quality_Control.RData"))

##### 03 Combine different data sets after QC #####
  ## Combine SeuObjs from list after QC  # (About 30 min for 20000 cells)
  scRNA.SeuObj <- CombineSeuObj(scRNA_SeuObj_QC.list)

  ## Check QC
  scRNAQC(scRNA.SeuObj,AddMitInf = "No",CheckOnly="Yes", Path = PathQC ,
          SpeciSet = Species, FileName = paste0(ProjectName,"_Check"))

  #### Save RData ####
  save.image(paste0(Save.Path,"/03_Combine_different_data_sets_after_QC.RData"))

##### 04 Perform an integrated analysis #####
  source("FUN_DRCluster.R")
  ## Create new folder
  PathCluster <- paste0(Save.Path,"/","A02_Cluster")
  if (!dir.exists(PathCluster)){dir.create(PathCluster)}

  scRNA.SeuObj <- DRCluster(scRNA.SeuObj, scRNA_SeuObj.list, seed=1, PCAdims = 30,
                            Path = PathCluster, projectName= ProjectName,
                            MetaSet = Ori_Meta.set)

  ##### Meta Table  #####
  Meta.df <- MetaSummary(scRNA_SeuObj.list, scRNA.SeuObj,
                         scRNA_SeuObj_QC.list, scRNA_Ori.SeuObj)

  #### Save RData ####
  save.image(paste0(Save.Path,"/04_Perform_an_integrated_analysis.RData"))

  ##****************************************************************************##
  ################## (Pending) ROUGE ##################
  ##****************************************************************************##

##### 05 Identify conserved cluster markers  #####
  ## Create new folder
  PathCellType <- paste0(Save.Path,"/","A03_CellTypeAno")
  if (!dir.exists(PathCellType)){dir.create(PathCellType)}

  ## Identify conserved cluster markers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  set.seed(1) # Fix the seed
  PBMC.markers.df <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(PBMC.markers.df, file = paste0(PathCluster,"/CC_ClusterMarker_AllGene.txt"),
              quote = F,sep = "\t",row.names = F)

  scRNA.SeuObj <- Beautify_Heatmap_Seurat(scRNA.SeuObj, PBMC.markers.df, topN = 7, Path = PathCluster,
                                          Type = "seurat_clusters",
                                          projectName = ProjectName)


  # --------------- Check specific tissue marker --------------- #
  #(Pending)
  UMAP_CellTypeMarker()

  #### Save RData ####
  save.image(paste0(Save.Path,"/05_Identify_conserved_cell_type_markers.RData"))

##### 06 Cell type annotation*  #####

  ##### 06 Cell type annotation: Traditional method (Manual annotation by cell type marker)  #####
  source("#_MagicDisc_06A_CellType_Anno_TM.R")

  #### Save RData ####
  save.image(paste0(Save.Path,"/06A_Cell_type_annotation_TM.RData"))


  ##### 06 Cell type annotation: Auto Cell type annotation  #####
  source("#_MagicDisc_06A_CellType_Anno_Auto.R")

  #### Save RData ####
  save.image(paste0(Save.Path,"/06A_Cell_type_annotation_Auto.RData"))

  ##### Verification (CellCheck) #####
  source("#_MagicDisc_06C_CellType_Anno_VER.R")

  # #### Save RData ####
  # save.image(paste0(Save.Path,"/06C_Cell_type_annotation_VER.RData"))

##### 07 Count Cell number  #####
  source("FUN_AnnoSummary.R")
  source("FUN_Export_CellCount.R")

  ## Create new folder
  PathCellCount <- paste0(Save.Path,"/","B01_CellCount")
  if (!dir.exists(PathCellCount)){dir.create(PathCellCount)}

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
  if (!dir.exists(PathBiomarkers)){dir.create(PathBiomarkers)}

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


  ##****************************************************************************##
  ################## (Pending) Trajectory inference ##################
  ##****************************************************************************##


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
  source("#_MagicDisc_DSA_04_RNADeconvolution.R")

  ##### save.image #####
  save.image(paste0(Save.Path,"/DSA_04_RNADeconvolution.RData"))

##### Clinical analysis #####


##### 013 inferCNV #####
  ## Create new folder
  PathinferCNV <- paste0(Save.Path,"/","D01_inferCNV")
  if (!dir.exists(PathinferCNV)){dir.create(PathinferCNV)}

  infercnv_obj <- inferCNV(scRNA.SeuObj, AnnoSet = "celltype",
                           SpeciSet = Species,
                           Path = PathinferCNV,
                           RefSet = c("T","B"),
                           CreateInfercnvObject.lt = list(chr_exclude = c("chrM")))
  ##### save.image #####
  save.image(paste0(Save.Path,"/013_inferCNV.RData"))


##### Beautify Figs #####
## Test


##**************************** Export the log file (End) ****************************##
  closeAllConnections() # Close connection to log file
##**************************** Export the log file (End) ****************************##

