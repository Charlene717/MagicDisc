# https://satijalab.org/seurat/articles/integration_introduction.html
# lupus
# Stim: (IFN)-??
#############
  rm(list = ls()) # Clean variable

  memory.limit(150000)


##### Setup the Seurat objects #####
  library(Seurat)
  library(SeuratData) # devtools::install_github('satijalab/seurat-data')
  library(patchwork)

  # install dataset
  InstallData("ifnb")

  # load dataset
  LoadData("ifnb")


  # split the dataset into a list of two seurat objects (stim and CTRL)
  ifnb.list <- SplitObject(ifnb, split.by = "stim")

  # ## Export Seurat Object in 10X format
  # # https://www.jianshu.com/p/771536c4eb5a
  # ifnb.CTRL.updated = UpdateSeuratObject(object = ifnb.list[["CTRL"]])
  # ifnb.STIM.updated = UpdateSeuratObject(object = ifnb.list[["STIM"]])
  #
  # # https://github.com/satijalab/seurat/issues/884
  # library(DropletUtils)
  # write10xCounts(ifnb.CTRL.updated@assays[["RNA"]]@counts , path = paste0(setwd(getwd()),"/CTRL"))
  # write10xCounts(ifnb.STIM.updated@assays[["RNA"]]@counts , path = paste0(setwd(getwd()),"/STIM"))
  #

  # normalize and identify variable features for each dataset independently
  set.seed(1) # Fix the seed
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  # select features that are repeatedly variable across datasets for integration
  set.seed(1) # Fix the seed
  features <- SelectIntegrationFeatures(object.list = ifnb.list)

##### Perform integration #####
  #(Seurat) https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8
  #(MNN) https://www.nature.com/articles/nbt.4091
  set.seed(1) # Fix the seed
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

  # this command creates an 'integrated' data assay
  set.seed(1) # Fix the seed
  immune.combined <- IntegrateData(anchorset = immune.anchors)


#### Perform an integrated analysis #####
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(immune.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  set.seed(1) # Fix the seed
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  set.seed(1) # Fix the seed
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  set.seed(1) # Fix the seed
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)


  # Visualization
  p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
  p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
  p1 + p2

  #
  DimPlot(immune.combined, reduction = "umap", split.by = "stim")


##### Identify conserved cell type markers #####
  # For performing differential expression after integration, we switch back to the original
  # data
  library("BiocManager")
  library("multtest")
  library("metap")
  DefaultAssay(immune.combined) <- "RNA"
  nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
  head(nk.markers)

  FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                            "CCL2", "PPBP"), min.cutoff = "q9")

  immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                  `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                  `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
  DimPlot(immune.combined, label = TRUE)

  ##
  nk.markers2 <- FindConservedMarkers(immune.combined, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
  head(nk.markers2)


  Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
                                                                        "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
                                                                        "CD4 Naive T", "CD4 Memory T"))
  markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                       "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                       "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
  DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
    RotatedAxis()


##### Identify differential expressed genes across conditions #####
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  t.cells <- subset(immune.combined, idents = "CD4 Naive T")
  Idents(t.cells) <- "stim"
  avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
  avg.t.cells$gene <- rownames(avg.t.cells)

  cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
  Idents(cd14.mono) <- "stim"
  avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
  avg.cd14.mono$gene <- rownames(avg.cd14.mono)

  genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
  p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
  p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
  p1 + p2


  immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
  immune.combined$celltype <- Idents(immune.combined)
  Idents(immune.combined) <- "celltype.stim"
  b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
  head(b.interferon.response, n = 15)


  FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3,
              cols = c("grey", "red"))


  plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",
                   pt.size = 0, combine = FALSE)
  wrap_plots(plots = plots, ncol = 1)

##### Performing integration on datasets normalized with SCTransform #####
  LoadData("ifnb")
  ifnb.list <- SplitObject(ifnb, split.by = "stim")
  ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
  ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                           anchor.features = features)
  immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

  immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
  immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)


  p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
  p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
                repel = TRUE)
  p1 + p2

##### Export Seurat Object in 10X format #####
  ## Ref: https://github.com/satijalab/seurat/issues/884
  # write10xCounts(x = seurat.object@assays$RNA@counts, path = output.path)
  library(DropletUtils)
  # write10xCounts(x = immune.combined.sct@assays$RNA@counts, path = "Seurat_Int")
  write10xCounts(x = ifnb.list[["CTRL"]]@assays$RNA@counts, path = "Seurat_Int/CTRL")
  write10xCounts(x = ifnb.list[["STIM"]]@assays$RNA@counts, path = "Seurat_Int/STIM")
