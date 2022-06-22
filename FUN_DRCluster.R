DRCluster <- function(scRNA_SeuObj.list, seed=1, PCAdims = 30,
                      UMAP_NNEighbors = 30, UMAP_MinDist = 0.3,
                      Path = Save.Path,
                      projectName= ProjectName,
                      MetaSet = Ori_Meta.set
                      ){

    # Run the standard workflow for visualization and clustering

    # # # !!
    # set.seed(seed) # Fix the seed
    # all.genes <- rownames(scRNA.SeuObj)
    # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, features = all.genes)

    ## Issues: re_clustering in seurat v3
    ## https://github.com/satijalab/seurat/issues/1528
    # DefaultAssay(scRNA.SeuObj) <- "RNA"

    # specify that we will perform downstream analysis on the corrected data note that the
    # original unmodified data still resides in the 'RNA' assay

    if(length(scRNA_SeuObj.list)==1){
      DefaultAssay(scRNA.SeuObj) <- "RNA"
      scRNA.SeuObj <- FindVariableFeatures(object = scRNA.SeuObj)
    }else{
      DefaultAssay(scRNA.SeuObj) <- "integrated"
    }

    ##?
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
    AnnoNames.set <- colnames(scRNA.SeuObj@meta.data)
    # ## Run if use filter
    # set.seed(seed) # Fix the seed
    # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)


    ### RunPCA
    # set.seed(seed) # Fix the seed
    # scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 30, verbose = FALSE)
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- RunPCA(scRNA.SeuObj,npcs = PCAdims, features = VariableFeatures(object = scRNA.SeuObj))

    print(scRNA.SeuObj[["pca"]], dims = 1:5, nfeatures = 5)

    pdf(
      file = paste0(Path,"/",projectName,"_PCA.pdf"),
      width = 10,  height = 8
    )
    VizDimLoadings(scRNA.SeuObj, dims = 1:2, reduction = "pca")
    DimPlot(scRNA.SeuObj, reduction = "pca")
    DimHeatmap(scRNA.SeuObj, dims = 1, cells = 500, balanced = TRUE)
    DimHeatmap(scRNA.SeuObj, dims = 1:15, cells = 500, balanced = TRUE)
    DimHeatmap(scRNA.SeuObj, dims = 16:30, cells = 500, balanced = TRUE)

    # # Determine the 'dimensionality' of the dataset
    # # NOTE: This process can take a long time for big datasets, comment out for expediency. More
    # # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
    # # computation time
    # scRNA.SeuObj <- JackStraw(scRNA.SeuObj, num.replicate = 100)
    # scRNA.SeuObj <- ScoreJackStraw(scRNA.SeuObj, dims = 1:20)
    # JackStrawPlot(scRNA.SeuObj, dims = 1:20)
    ElbowPlot(scRNA.SeuObj, ndims = PCAdims)
    dev.off()

    ElbowPlot(scRNA.SeuObj, ndims = PCAdims)

    ## Issues: RunUMAP causes R exit
    ## https://github.com/satijalab/seurat/issues/2259
    # The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    # To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    # This message will be shown once per session
    #### UMAP
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:PCAdims ,
                            n.neighbors = UMAP_NNEighbors, min.dist= UMAP_MinDist)

    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:PCAdims)
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

    #### tSNE
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- RunTSNE(scRNA.SeuObj, reduction = "pca", dims = 1:PCAdims)
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:PCAdims)
    set.seed(seed) # Fix the seed
    scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)


    ## Visualization
    DimPlot(scRNA.SeuObj, reduction = "umap", group.by = colnames(list_files.df)[3] ) %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)

    pdf(
      file = paste0(Path,"/",projectName,"_nlDR_Cluster.pdf"),
      width = 12,  height = 8
    )

    DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>%
      BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)


    for (i in 1:(length(MetaSet)-3)) {
      print(DimPlot(scRNA.SeuObj, reduction = "umap", group.by = AnnoNames.set[i+3]) %>%
              BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+
              theme(plot.title = element_text(vjust = 0.85)))
      print(DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2, split.by = AnnoNames.set[i+3], label = TRUE, label.size = 4) %>%
              BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                             SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))

      print(DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = AnnoNames.set[i+3]) %>%
              BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+
              theme(plot.title = element_text(vjust = 0.85)))
      print(DimPlot(scRNA.SeuObj, reduction = "tsne", ncol = 2, split.by = AnnoNames.set[i+3], label = TRUE, label.size = 4) %>%
              BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                             SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))

    }
    rm(i)

    dev.off()    # graphics.off()



    return(scRNA.SeuObj)
}
