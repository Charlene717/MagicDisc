## Ref: https://statisticsglobe.com/r-save-all-console-input-output-to-file
## Ref: https://blog.gtwang.org/r/r-data-input-and-output/

#************************************************************************************************************************#
##### To-Do List ######
  # - [ ] Ensembl gene name conversion
  # - [ ] Clean up the code

  # - [ ] Gene correlation
  # - [ ] Delete unrecognized cell type

  # - [ ] Find Biomarker
  # - [ ] GSEA
  # - [ ] GO
  # - [ ] Cell-cell interaction
  # - [ ] inferCNV
  # - [ ] Deconvolution
  # - [ ] Clinical analysis
  # - [ ] Beautify Figs

#************************************************************************************************************************#

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)
# options(stringsAsFactors = FALSE)
# Sys.setlocale(category = "LC_ALL", locale = "UTF-8")

##### Setting ######
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

##### Export the log file (Start) #####
  ## Create new folder for log file
  if (!dir.exists("LogFiles")){dir.create("LogFiles")}

  my_log <- file(paste0("LogFiles/",Version,"_log.txt")) # File name of output log

  sink(my_log, append = TRUE, type = "output") # Writing console output to log file
  sink(my_log, append = TRUE, type = "message")

  cat(readChar(rstudioapi::getSourceEditorContext()$path, # Writing currently opened R script to file
               file.info(rstudioapi::getSourceEditorContext()$path)$size))


##### Load Packages #####
FUN_Basic.set <- c("tidyverse","Seurat","monocle","ggplot2","ggpmisc","broom",
                   "stringr","magrittr","dplyr",
                   "patchwork","reticulate","anndata")
FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb",
                         "basilisk","zellkonverter","SeuratDisk",
                         "SingleR","scRNAseq","celldex","scran")
source("#_MagicDisc_PKG_FUN.R")

##### Load datasets  #####
  ## Annotation table
  list_files.df <- read.csv(paste0(InputFolder,"/",InputAnno))
  Feature.set <- colnames(list_files.df)[-1]

  ## Read 10x files
  source("FUN_ReadscRNA.R")
  scRNA_SeuObj.list <- FUN_ReadscRNA(InputFolder,Folder = Folder,
                                     Path =  "/monocle/outs/filtered_gene_bc_matrices/mm10",
                                     list_files.df, Mode = DataMode, projectName = ProjectName)

##### 01 Combine different datasets before QC #####
  source("FUN_Cal_Mit.R")
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
  ## Combine SeuObjs from list after QC
  # (About 30 min for 20000 cells)
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
  if (!dir.exists(PathCluster)){
    dir.create(PathCluster)
  }

  scRNA.SeuObj <- DRCluster(scRNA.SeuObj, scRNA_SeuObj.list, seed=1, PCAdims = 30,
                            Path = PathCluster, projectName= ProjectName,
                            MetaSet = Ori_Meta.set)

  ##### Meta Table  #####
  Meta.df <- MetaSummary(scRNA_SeuObj.list, scRNA.SeuObj,
                         scRNA_SeuObj_QC.list,scRNA_Ori.SeuObj)

  #### Save RData ####
  save.image(paste0(Save.Path,"/04_Perform_an_integrated_analysis.RData"))

################## (Pending) Cell Cycle Regression ##################

##### 05 Identify conserved cluster markers  #####
  ## Create new folder
  PathCellType <- paste0(Save.Path,"/","A03_CellTypeAno")
  if (!dir.exists(PathCellType)){
    dir.create(PathCellType)
  }

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
  scRNA.SeuObj <- Beautify_Heatmap_Seurat(scRNA.SeuObj, PBMC.markers.df, topN = 7, Path = PathCellType,
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

  ## Create cell type markers dataframe
  # DefaultAssay(scRNA.SeuObj_Small) <- "RNA"
  scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
  CellType.markers.df <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

  write.table(CellType.markers.df, file = paste0(PathCellType,"/CC_CelltypeMarker_AllGene.txt"),
              quote = F,sep = "\t",row.names = F)

  #### Save RData ####
  save.image(paste0(Save.Path,"/06_Cell_type_annotation.RData"))


################## (Pending) Auto Cell type annotation ##################
##### 06 Cell type annotation  #####
  ##### scSorter #####
    library(scSorter)
    ## Ref: https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

    #### Example ####
    load(url('https://github.com/hyguo2/scSorter/blob/master/inst/extdata/TMpancreas.RData?raw=true'))
    anno_ori <- anno
    anno <- anno[!anno$Type %in% c("Pancreatic_Acinar_cells","Pancreatic_PP_cells"),]
    # anno$Marker <- toupper(anno$Marker)

    #### Test ####
    ## Create anno.df
    CellType.markers.df %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> CTTop.markers
    #DoHeatmap(scRNA.SeuObj, features = CTTop.markers$gene) + NoLegend()
    anno.df <- data.frame(Type = CTTop.markers$cluster,
                          Marker = CTTop.markers$gene,
                          Weight = CTTop.markers$avg_log2FC)

    ## Create small sample for test
    scRNA.SeuObj_Small <- scRNA.SeuObj[,scRNA.SeuObj$cells %in% sample(scRNA.SeuObj$cells,1000)]
    scSorter.obj <- scSorter(scRNA.SeuObj_Small@assays[["RNA"]]@counts %>% as.data.frame(),anno.df)

    scRNA.SeuObj_Small$scSorterPred <- scSorter.obj[["Pred_Type"]]
    DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="scSorterPred" ,label = TRUE, pt.size = 0.5) + NoLegend()
    DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()

  ##### SingleR #####
    ## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
    ## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html
    library(SingleR)

    #### Example ####
    library(scRNAseq)
    hESCs <- LaMannoBrainData('human-es')
    hESCs <- hESCs[,1:100]

    library(celldex)
    hpca.se <- HumanPrimaryCellAtlasData()
    hpca.se

    singler.results.Demo <- SingleR(scRNA.cds_Small, hpca.se, labels = hpca.se$label.main)

    #### Test ####
    Temp <- as.SingleCellExperiment(scRNA.SeuObj_Small)
    # Ref <- HumanPrimaryCellAtlasData()
    Ref <- ImmGenData()
    singler.results <- SingleR(Temp, Ref, labels = Ref$label.main)
    scRNA.SeuObj_Small[["SingleR.labels"]] <- singler.results$labels
    DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="SingleR.labels" ,label = TRUE, pt.size = 0.5) + NoLegend()
    DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()

    ## Create reference
    # ## Ref: https://github.com/dviraran/SingleR/issues/136
    # ## Ref: http://bioconductor.org/books/3.15/OSCA.basic/cell-type-annotation.html
    # ### Ref: https://rdrr.io/github/dviraran/SingleR/f/vignettes/SingleR_create.Rmd
    # library(scran)
    #
    # library(scRNAseq)
    # sceM <- MuraroPancreasData()
    # out <- pairwiseBinom(counts(sceM), sceM$label, direction="up")
    # markers <- getTopMarkers(out$statistics, out$pairs, n=10)
    sce.muraro


  ##### CelliD #####
  ## Ref: https://github.com/RausellLab/CelliD
  ## CelliD Vignette Ref: https://bioconductor.org/packages/release/bioc/vignettes/CelliD/inst/doc/BioconductorVignette.html

    ## install
    if(!require("tidyverse")) install.packages("tidyverse")
    if(!require("ggpubr")) install.packages("ggpubr")
    BiocManager::install("CelliD")

    ## Load libraries
    library(CelliD)
    library(tidyverse) # general purpose library for data handling
    library(ggpubr) #library for plotting

    ## Obtaining pancreatic cell-type gene signatures
    panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

    # restricting the analysis to pancreas specific gene signatues
    panglao_pancreas <- panglao %>% dplyr::filter(organ == "Pancreas")

    library(Hmisc)
    capitalize(tolower("TOP2A"))
    capitalize(tolower("TOP2A"))
    panglao_pancreas$`official gene symbol` <- capitalize(tolower(panglao_pancreas$`official gene symbol`))

    # restricting to human specific genes
    panglao_pancreas <- panglao_pancreas %>%  dplyr::filter(str_detect(species,"Hs"))

    # converting dataframes into a list of vectors, which is the format needed as input for CelliD
    panglao_pancreas <- panglao_pancreas %>%
      group_by(`cell type`) %>%
      summarise(geneset = list(`official gene symbol`))
    pancreas_gs <- setNames(panglao_pancreas$geneset, panglao_pancreas$`cell type`)

    ## Obtaining gene signatures for all cell types in the Panglao database
    #filter to get human specific genes
    panglao_all <- panglao %>%  dplyr::filter(str_detect(species,"Hs"))

    # convert dataframes to a list of named vectors which is the format for CelliD input
    panglao_all <- panglao_all %>%
      group_by(`cell type`) %>%
      summarise(geneset = list(`official gene symbol`))
    all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)

    #remove very short signatures
    all_gs <- all_gs[sapply(all_gs, length) >= 10]

    ## Assessing per-cell gene signature enrichments against pre-established marker lists
    # Performing per-cell hypergeometric tests against the gene signature collection
    scRNA.SeuObj_Small <- RunMCA(scRNA.SeuObj_Small)
    DimPlotMC(scRNA.SeuObj_Small, reduction = "mca", group.by = "celltype", features = c("CTRL", "INS", "MYZAP", "CDH11"), as.text = TRUE) + ggtitle("MCA with some key gene markers")

    HGT_pancreas_gs <- RunCellHGT(scRNA.SeuObj_Small, pathways = pancreas_gs, dims = 1:50, n.features = 200)
    # For each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
    pancreas_gs_prediction <- rownames(HGT_pancreas_gs)[apply(HGT_pancreas_gs, 2, which.max)]

    # For each cell, evaluate if the lowest p-value is significant
    pancreas_gs_prediction_signif <- ifelse(apply(HGT_pancreas_gs, 2, max)>2, yes = pancreas_gs_prediction, "unassigned")

    # Save cell type predictions as metadata within the Seurat object
    scRNA.SeuObj_Small$pancreas_gs_prediction <- pancreas_gs_prediction_signif


    DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="pancreas_gs_prediction" ,label = TRUE, pt.size = 0.5) + NoLegend()
    DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()

  ##### Garnett #####
    library(monocle3)
    library(garnett)

    ## Convert Seurat object to Monocle3
    ## Ref: https://github.com/satijalab/seurat-wrappers
    remotes::install_github('satijalab/seurat-wrappers')
    library(SeuratWrappers)
    scRNA.cds_Small <- as.cell_data_set(scRNA.SeuObj_Small)

    ## Ref: https://github.com/satijalab/seurat-wrappers/issues/54
    ## Calculate size factors using built-in function in monocle3
    scRNA.cds_Small <- estimate_size_factors(scRNA.cds_Small)

    ## Add gene names into CDS
    scRNA.cds_Small@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(scRNA.SeuObj_Small[["RNA"]])

    ## Automated annotation with Garnett
    assigned_type_marker_test_res <- top_markers(scRNA.cds_Small,
                                                 group_cells_by="celltype",
                                                 #reference_cells=1000,
                                                 cores=4)

  ##### Verification (CellCheck) #####
    #### Install ####
    ## Check whether the installation of those packages is required
    Package.set <- c("tidyverse","caret","cvms","DescTools","devtools")
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        install.packages(Package.set[i])
      }
    }
    ## Load Packages
    # library(Seurat)
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)

    ## install CellCheck
    # Install the CellCheck package
    detach("package:CellCheck", unload = TRUE)
    devtools::install_github("Charlene717/CellCheck")
    # Load CellCheck
    library(CellCheck)

    #### Run CellCheck ####
    ## Create check dataframe
    scSorter.df <- data.frame(Actual = scRNA.SeuObj_Small@meta.data[["celltype"]],
                              Predict = scRNA.SeuObj_Small@meta.data[["scSorterPred"]])
    scSorter_Anno.df <- data.frame(TestID = "Predict",
                                   Tool = "scSorter",
                                   Type = "PDAC",
                                   PARM = "1")
    ## For one prediction
    DisMultCM.lt <- list(Actual = "Actual", Predict = "Predict", Type = "Type", Type2 = "PDAC" )
    cm_DisMult.lt <- CellCheck_DisMult(scSorter.df, scSorter_Anno.df, Mode = "One", DisMultCM.lt,
                                       Save.Path = Save.Path, ProjectName = ProjectName)
    ## For multiple prediction
    Sum_DisMult.df <- CellCheck_DisMult(scSorter.df, scSorter_Anno.df,
                                        Mode = "Multiple",DisMultCM.lt=DisMultCM.lt,
                                        Save.Path = Save.Path, ProjectName = ProjectName)

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


##### 013 inferCNV #####
  ## Create new folder
  PathinferCNV <- paste0(Save.Path,"/","D01_inferCNV")
  if (!dir.exists(PathinferCNV)){
    dir.create(PathinferCNV)
  }

  infercnv_obj <- inferCNV(scRNA.SeuObj, AnnoSet = "celltype",
                           SpeciSet = Species,
                           Path = PathinferCNV,
                           RefSet = c("T","B"),
                           CreateInfercnvObject.lt = list(chr_exclude = c("chrM")))
  ##### save.image #####
  save.image(paste0(Save.Path,"/013_inferCNV.RData"))


##### Beautify Figs #####
## Test

##### Export the log file (End) #####
  closeAllConnections() # Close connection to log file

