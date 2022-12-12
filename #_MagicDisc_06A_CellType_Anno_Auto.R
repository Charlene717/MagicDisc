##### 06 Cell type annotation*  #####
##### 06 Cell type annotation: Auto Cell type annotation  #####
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

# #### Save RData ####
# save.image(paste0(Save.Path,"/06A_Cell_type_annotation_Auto.RData"))

