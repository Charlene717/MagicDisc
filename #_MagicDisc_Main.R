## Ref: https://statisticsglobe.com/r-save-all-console-input-output-to-file
## Ref: https://blog.gtwang.org/r/r-data-input-and-output/

# setwd("../") # Set the path at the same location of Demo_CellTypeAnno.R
my_log <- file("CellCheck_log.txt") # File name of output log

sink(my_log, append = TRUE, type = "output") # Writing console output to log file
sink(my_log, append = TRUE, type = "message")

cat(readChar(rstudioapi::getSourceEditorContext()$path, # Writing currently opened R script to file
             file.info(rstudioapi::getSourceEditorContext()$path)$size))



##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####

  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc","broom",
                   "stringr","magrittr","dplyr",
                   "CellChat","patchwork","reticulate","anndata")
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
  source("FUN_Export_All_DRPlot")
  source("FUN_BeautifyDotPlot.R")
  source("FUN_BioMarker2Index.R")
  source("FUN_BeautifyVennDiag.R")

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
  scRNA.SeuObj <- Beautify_Heatmap_Seurat(scRNA.SeuObj, PBMC.markers, topN = 7, Path = PathCluster,
                                          Type = "celltype", HMColor.lt = Heatmap_Color.lt,
                                          projectName = ProjectName)

  ## Export All DRPlot(UMAP,tSNE)
  Export_All_DRPlot(scRNA.SeuObj)

  ## DotPlot
  #(Pending)
  BeautifyDotPlot(scRNA.SeuObj, Path = PathCellType, projectName = ProjectName)


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

  scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)

  ## Find CCmarker in different Cell type
  CCMarker.lt <-  BioMarker2Index(scRNA.SeuObj, Path = PathBiomarkers, projectName = ProjectName,
                                  classSet2 = ClassSet2, classSet3 = ClassSet3,Type = "celltype")

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_1_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Separate).RData"))


##### 08_2 Find CCmarker in different Cell type and VennDiagrame (SSA_IntersectCT) ########
  ##-------------- Intersect_CellType --------------##
  CCMarker_Male.lt <- CCMarker.lt[1]
  CCMarker_Female.lt <- CCMarker.lt[2]
  intersect_CellType <- intersect(names(CCMarker_Male.lt),names(CCMarker_Female.lt))

  CCMarker_Male.lt <- CCMarker_Male.lt[names(CCMarker_Male.lt) %in% intersect_CellType]
  CCMarker_Female.lt <- CCMarker_Female.lt[names(CCMarker_Female.lt) %in% intersect_CellType]

  CellType.list <- names(CCMarker_Male.lt)

  Venn_CCMarke.lt <- BeautifyVennDiag(CCMarker_Male.lt, CCMarker_Female.lt, CellType.list)

  #### Save RData ####
  save.image(paste0(Save.Path,"/08_2_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_Venn.RData"))


##### 08_3 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")

  CCMarker.lt <- BioMarker1Index(scRNA.SeuObj, Path = PathBiomarkers, projectName = ProjectName,
                                  sampletype = Sampletype, cellType.list = CellType.list, classSet2 = ClassSet2,
                                  Type = paste0("celltype.",classSet2) )

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

  # Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  Pathway.all <- read.delim2(paste0(getwd(),"/",InputGSEA),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",InputGSEA))),
                             header = F,sep = "\t")

  ##### Converting the Human gene name to Mouse gene name #####
    # #  Need to be optimized
    # # (Method1) bind the different length of column (Cannot use rbind)
    # # (Method2) Save the data as list first and than use do.call() to unlist to have dataframe
    #
    # ## (Ori method)
    # Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*2))
    # for (i in 1:nrow(Pathway.all)) {
    #   #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    #   PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()
    #   colnames(PathwayN)="Temp"
    #   PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
    #   Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
    # }
    #
    # Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
    # colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
    #
    # rm(PathwayN)
    #
    # # assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
    # # assign(colnames(Pathway.all)[i],Pathway.all[,i])

    ## (Method1)
    # Refer # https://stackoverflow.com/questions/3699405/how-to-cbind-or-rbind-different-lengths-vectors-without-repeating-the-elements-o
    # How to cbind or rbind different lengths vectors without repeating the elements of the shorter vectors?
    ## Modify by Charlene: Can use in MultRow
    bind_diff <- function(x, y){
      if(ncol(x) > ncol(y)){
        len_diff <- ncol(x) - ncol(y)
        y <- data.frame(y, rep(NA, len_diff) %>% t() %>% as.data.frame())
        colnames(x) <- seq(1:ncol(x))
        colnames(y) <- seq(1:ncol(y))
      }else if(ncol(x) < ncol(y)){
        len_diff <- ncol(y) - ncol(x)
        x <- data.frame(x, rep(NA, len_diff) %>% t() %>% as.data.frame())
        colnames(x) <- seq(1:ncol(x))
        colnames(y) <- seq(1:ncol(y))
      }
      rbind(x, y)
    }

    ## Converting
    for (i in 1:nrow(Pathway.all)) {
      PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()  %>% as.data.frame()
      colnames(PathwayN)="Temp"
      PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
      PathwayN <- PathwayN[!PathwayN$MM.symbol  == 0,]
      PathwayN <- PathwayN[!is.na(PathwayN$MM.symbol),]
      if(i==1){
        Pathway.all.MM <- unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame()
      }else{
        Pathway.all.MM <- bind_diff(Pathway.all.MM,unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame())
        # Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
      }
    }

    Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
    colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
    #Pathway.all.MM[Pathway.all.MM==0] <-NA

    rm(PathwayN)

save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))

##### 09_1 GSEA Analysis (SPA) #####
  ## Create folder
  Sep_Cla3_GSEA.Path <- paste0(Sampletype,"_",ProjectName,"_GSEA")
  dir.create(paste0(Save.Path,"/",Sep_Cla3_GSEA.Path))


  GSEA_Large <- list()
  GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large.df.TOP <- GSEA_Large.df



  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,".pdf"),width = 15, height = 7 )

  for(i in 1:length(CellType.list)){

    gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))

    GSEA_Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)

    fgseaRes <- GSEA_Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    pathwaysH <- GSEA_Large.Output[["Pathway.all.list"]]

    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

    topPathways <- GSEA_Large.Output[["topPathways"]]

    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway],
                  ranks,
                  fgseaRes,
                  gseaParam = 0.5) + title( paste0(Sampletype,".",ProjectName,".",CellType.list[i]), adj = 0, line =3)

    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+
                           labs(title= paste0(ampletype,".",ProjectName,".",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+
                           labs(title= paste0(ampletype,".",ProjectName,".",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1

    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large[[i]] <- Sum
    names(GSEA_Large)[[i]] <- paste0(CellType.list[i])

    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large.df <- rbind(GSEA_Large.df,fgseaRes2 )

    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP, topPathways2)

    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

  }

  dev.off()

  ## GSEA_Large.Sum.TOP ##
  GSEA_Large.Sum.TOP <- rbind(GSEA_Large.df.TOP)
  GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large.Sum.TOP, file=paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_Pathway_LargeTOP_SPA.txt"),sep="\t",
              row.names=F, quote = FALSE)

  ##### Bubble plot #####
  library(ggplot2)
  library(scales)
  GSEA_Color.lt = list(high = "#ef476f",mid = "white",low = "#0077b6")

  GSEA_Large.Sum.TOP$PhenoType <- factor(GSEA_Large.Sum.TOP$PhenoType,
                                         levels = Cell_Type_Order.set)

  GSEA_ggplot_SPA.lt <- GSEA_ggplot(GSEA_Large.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
  GSEA_Large.Sum.TOP.S <- GSEA_ggplot_SPA.lt[["GSEA_TOP.df"]]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]

  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_Bubble_SPA.pdf"),width = 17, height = 12 )
    GSEA_ggplot_SPA.lt[["BBPlot_Ori"]]
    GSEA_ggplot_SPA.lt[["BBPlot"]]
    GSEA_ggplot_SPA.lt[["BBPlot2"]]
    GSEA_ggplot_SPA.lt[["BBPlotB1"]]
    GSEA_ggplot_SPA.lt[["BBPlotB1"]]
  dev.off()


  ##### Extract SubType #####
  ## T Cell
  # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
  GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]

  BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
    geom_point() +
    scale_size_area(max_size = 7)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  BBPlot_T

  BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                           XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  # BBPlot_TB <- BBPlot_TB +theme(axis.title.y=element_blank(),
  #                  axis.text.y=element_blank(),
  #                  axis.ticks.y=element_blank())
  BBPlot_TB

  BBPlot_TB1 <- BBPlot_TB %>%
    insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
  BBPlot_TB1


  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_Bubble_SPA_SubType_T.pdf"),width = 17, height = 7 )
  BBPlot_TB
  BBPlot_TB1
  dev.off()


  ## Mac
  GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]

  BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
    geom_point() +
    scale_size_area(max_size = 5)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  BBPlot_Mac

  BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                               XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)

  BBPlot_MacB1 <- BBPlot_MacB %>%
    insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
  BBPlot_MacB1

  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_Bubble_SPA_SubType_Mac.pdf"),width = 17, height = 20 )
    BBPlot_MacB
    BBPlot_MacB1
    dev.off()

  rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
     df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

  ##### save.image #####
  save.image(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

##### 09_2 GSEA Analysis (SSA_MAle) #####
  GSEA_Large_Male <- list()
  GSEA_Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large_Male.df.TOP <- GSEA_Large_Male.df


  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_SSA_Male.pdf"),width = 15, height = 7 )

  for(i in 1:length(CellType.list)){

    gseaDat <- CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))


    #GSEA_Large_Male.Output <- FUN_GSEA_Large_MaleGeneSet(ranks,Pathway.all,10)
    GSEA_Large_Male.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)

    fgseaRes <- GSEA_Large_Male.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    pathwaysH <- GSEA_Large_Male.Output[["Pathway.all.list"]]

    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

    topPathways <- GSEA_Large_Male.Output[["topPathways"]]

    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway],
                  ranks,
                  fgseaRes,
                  gseaParam = 0.5) + title( paste0(Sampletype,".",ProjectName,".",CellType.list[i]), adj = 0, line =3)

    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0(Sampletype,".",ProjectName,".",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+
                           labs(title= paste0(Sampletype,".",ProjectName,".",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1

    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large_Male[[i]] <- Sum
    names(GSEA_Large_Male)[[i]] <- paste0(CellType.list[i])

    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large_Male.df <- rbind(GSEA_Large_Male.df,fgseaRes2 )

    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large_Male.df.TOP <- rbind(GSEA_Large_Male.df.TOP,topPathways2 )

    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

  }

  dev.off()

  ## GSEA_Large_Male.Sum.TOP ##
  GSEA_Large_Male.Sum.TOP <- rbind(GSEA_Large_Male.df.TOP)
  GSEA_Large_Male.Sum.TOP <- GSEA_Large_Male.Sum.TOP[,!colnames(GSEA_Large_Male.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large_Male.Sum.TOP, file=paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"Pathway_LargeTOP_SSA_Male.txt"),sep="\t",
              row.names=F, quote = FALSE)



  ##### Bubble plot #####
  library(ggplot2)
  library(scales)

  GSEA_Large_Male.Sum.TOP$PhenoType <- factor(GSEA_Large_Male.Sum.TOP$PhenoType,
                                              levels = Cell_Type_Order.set)

  GSEA_ggplot_SSA_Male.lt <- GSEA_ggplot(GSEA_Large_Male.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
  GSEA_Large_Male.Sum.TOP.S <- GSEA_ggplot_SSA_Male.lt[["GSEA_TOP.df"]]

  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_Bubble_SSA_Male.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SSA_Male.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot2"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
  dev.off()


  ##### Extract SubType #####

  ## T Cell
  # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
  GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]

  BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
    geom_point() +
    scale_size_area(max_size = 7)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  BBPlot_T

  BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                           XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  BBPlot_TB

  BBPlot_TB1 <- BBPlot_TB %>%
    insert_left(GSEA_ggplot_SSA_Male.lt[["Y_Order"]],width = 0.2)
  BBPlot_TB1


  pdf(file = paste0(Save.Path, "/",Sep_Cla3_GSEA.Path,"/",Sep_Cla3_GSEA.Path,"_Bubble_SSA_Male_SubType_T.pdf"),width = 17, height = 7 )
  BBPlot_TB
  BBPlot_TB1
  dev.off()



##### GO/Metascape #####



##### inferCNV #####

##### Deconvolution #####

##### Beautify Figs #####

closeAllConnections() # Close connection to log file

