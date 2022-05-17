## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

CellChatOne <- function(seuratObject,
                        signalingtype = "ECM-Receptor", projectName = "ECM",
                        save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
                        groupby = "celltype", species = "human" # species = c("human","mouse")
                        ){

#### Load the required libraries ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","CellChat","patchwork","reticulate","anndata","Seurat","NMF","ggalluvial")
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
  Package.set <- c("basilisk","zellkonverter","SeuratDisk")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)

##### Current path and new folder setting*  #####
  SignalingType = signalingtype # Secreted Signaling, ECM-Receptor, Cell-Cell Contact
  ProjectName = projectName # Secret, ECM, CC
  Save_Path = save.path

  ## Create new folder
  if (!dir.exists(Save_Path)){
    dir.create(Save_Path)
  }


##### Part I: Data input & processing and initialization of CellChat object #####
  ## Extract the CellChat input files from a Seurat V3 object
  # Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
  data.input <- GetAssayData(seuratObject, assay = "RNA", slot = "data") # normalized data matrix
  meta <- seuratObject@meta.data # create a dataframe of the cell labels

  #### Create a CellChat object ####
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = groupby)

  #### Set the ligand-receptor interaction database ####
  if(species == "mouse"){
    CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data

  }else{
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

  }

  showDatabaseCategory(CellChatDB)

  pdf(file = paste0(Save_Path,"/",ProjectName,"_CellChatDB.pdf"),
      width = 7,  height = 7
  )
  p <- showDatabaseCategory(CellChatDB)
  print(p)
  dev.off()
  rm(p)

  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)

  # use a subset of CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  # CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use ECM-Receptor Signaling
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Cell-Cell Contact Signaling
  CellChatDB.use <- subsetDB(CellChatDB, search = SignalingType) # use Secreted Signaling

  # use all CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- CellChatDB # simply use the default CellChatDB

  # set the used database in the object
  cellchat@DB <- CellChatDB.use

  #### Preprocessing the expression data for cell-cell communication analysis ####
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4) # do parallel

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.human)

##### Part II: Inference of cell-cell communication network #####
  #### Compute the communication probability and infer cellular communication network ####
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  #### Infer the cell-cell communication at a signaling pathway level ####
  cellchat <- computeCommunProbPathway(cellchat)

  #### Calculate the aggregated cell-cell communication network ####
  cellchat <- aggregateNet(cellchat)

  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


  pdf(file = paste0(Save_Path,"/",ProjectName,"_interaction_strength.pdf"),
      width = 7,  height = 7
  )
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

  dev.off()


  pdf(file = paste0(Save_Path,"/",ProjectName,"_interaction_strength_Sep.pdf"),
      width = 17,  height = 15
  )
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

  mat <- cellchat@net$weight
  par(mfrow = c(3,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  rm(i)

##### Part III: Visualization of cell-cell communication network #####
  ## Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
  PathPN <- paste0(Save_Path,"/Pathway_Network")
  ## Create new folder
  if (!dir.exists(PathPN)){
    dir.create(PathPN)
  }

  pathway.set <- cellchat@netP[["pathways"]]

  #### Plot Summary ####
  pdf(file = paste0(PathPN,"/",ProjectName,"_pathway_network.pdf"),
      width = 7,  height = 7
  )

  for (i in 1:length(pathway.set)) {
    pathways.show <- pathway.set[i] # pathways.show <- c("CXCL")
    # Hierarchy plot
    # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
    vertex.receiver = seq(1,4) # a numeric vector.
    netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
    # # Circle plot
    # par(mfrow=c(1,1))
    # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

    # Heatmap
    par(mfrow=c(1,1))
    Heatmap <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    print(Heatmap)
    #> Do heatmap based on a single object


    # Chord diagram
    par(mfrow=c(1,1))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

    ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
    #par(mfrow=c(1,1))
    p <- netAnalysis_contribution(cellchat, signaling = pathways.show)
    print(p)

  }
  dev.off()
  rm(i,p,Heatmap)

  #### Plot all ####
  for (i in 1:length(pathway.set)) {
    pathways.show <- pathway.set[i] # pathways.show <- c("CXCL")

    pdf(file = paste0(PathPN,"/",ProjectName,"_pathway_network_",pathways.show,".pdf"),
        width = 7,  height = 7
    )


    # Hierarchy plot
    # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
    vertex.receiver = seq(1,4) # a numeric vector.
    netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
    # # Circle plot
    # par(mfrow=c(1,1))
    # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

    # Heatmap
    par(mfrow=c(1,1))
    Heatmap <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    print(Heatmap)
    #> Do heatmap based on a single object


    # Chord diagram
    par(mfrow=c(1,1))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

    # # Chord diagram
    # group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
    # names(group.cellType) <- levels(cellchat@idents)
    # netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
    #
    # #> Plot the aggregated cell-cell communication network at the signaling pathway level
    # #> Note: The first link end is drawn out of sector 'Inflam. FIB'.


    ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
    #par(mfrow=c(1,1))
    p <- netAnalysis_contribution(cellchat, signaling = pathways.show)
    print(p)

    pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
    for (j in 1:nrow(pairLR.CXCL)) {
      LR.show <- pairLR.CXCL[j,] # show one ligand-receptor pair
      # # Hierarchy plot
      # vertex.receiver = seq(1,4) # a numeric vector
      # netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
      #> [[1]]
      # Circle plot
      netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
      # Chord diagram
      netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

    }
    rm(j)

    dev.off()
  }
  rm(i,p,Heatmap)

  #### Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ####
  PathLR <- paste0(Save_Path,"/LRPair")
  ## Create new folder
  if (!dir.exists(PathLR)){
    dir.create(PathLR)
  }

  #### Bubble plot ####
  ## Sum
  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_Bubble_Sum.pdf"),
      width = 15,  height = 20
  )
  netVisual_bubble(cellchat,  remove.isolate = FALSE)
  dev.off()


  ## All
  try({
  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_Bubble_All.pdf"),
      width = 5,  height = 12
  )
  for (i in 1:ncol(mat)) {
    try({

      # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
      P <- netVisual_bubble(cellchat, sources.use = i, remove.isolate = FALSE)
      print(P)
      #> Comparing communications on a single object

    })
  }
  dev.off()
  })
  rm(i,p)

  #### Chord diagram ####
  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_ChordDiagram.pdf"),
      width = 12,  height = 12
  )
  for (i in 1:ncol(mat)) {
    try({

      # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
      P1 <-netVisual_chord_gene(cellchat, sources.use = i,  lab.cex = 0.5,legend.pos.y = 30)
      print(P1)
      P2 <-netVisual_chord_gene(cellchat, targets.use =i , lab.cex = 0.5,legend.pos.y = 30)
      print(P2)
      #> Comparing communications on a single object

    })
  }
  dev.off()
  rm(i,P1,P2)

  #### Plot the signaling gene expression distribution using violin/dot plot ####
  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_Violin.pdf"),
      width = 12,  height = 8
  )
  for (i in 1:length(pathway.set)) {
    pathways.show <- pathway.set[i] # pathways.show <- c("CXCL")

    P <- plotGeneExpression(cellchat, signaling = pathways.show)
    print(P)
  }
  dev.off()
  rm(i,P)

##### Part IV: Systems analysis of cell-cell communication network #####
  ##### Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling #####
  #### Compute and visualize the network centrality scores ####
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_NWCentScores_Heatmap.pdf"),
      width = 12,  height = 8
  )
  for (i in 1:length(pathway.set)) {
    pathways.show <- pathway.set[i] # pathways.show <- c("CXCL")

    P <-  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
    print(P)
  }
  dev.off()
  rm(i,P)



  #### Visualize the dominant senders (sources) and receivers (targets) in a 2D space ####
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat) + ggtitle("All signaling pathways")+
    theme(plot.title = element_text(color="black", size=14, face="bold"))
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest
  gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show) +
    ggtitle(paste0(pathways.show," signaling pathway network"))+
    theme(plot.title = element_text(color="black", size=14, face="bold"))
  #> Signaling role analysis on the cell-cell communication network from user's input
  gg1 + gg2

  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_SourcesTargets_2Dspace.pdf"),
      width = 7,  height = 7
  )
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat) + ggtitle("All signaling pathways")+
    theme(plot.title = element_text(color="black", size=14, face="bold"))
  gg1

  for (i in 1:length(pathway.set)) {
    pathways.show <- pathway.set[i] # pathways.show <- c("CXCL")
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    # Signaling role analysis on the cell-cell communication networks of interest
    gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show) +
      ggtitle(paste0(pathways.show," signaling pathway network"))+
      theme(plot.title = element_text(color="black", size=14, face="bold"))
    print(gg2)
    rm(gg2)
  }
  dev.off()
  rm(i,gg1)


  #### Identify signals contributing most to outgoing or incoming signaling of certain cell groups ####
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht1 + ht2
  # Signaling role analysis on the cell-cell communication networks of interest
  ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show)
  ht

  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_mostOutIn_Heatmap.pdf"),
      width = 12,  height = 8
  )
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht1 + ht2

  for (i in 1:length(pathway.set)) {

    pathways.show <- pathway.set[i] # pathways.show <- c("CXCL")

    # Signaling role analysis on the cell-cell communication networks of interest
    try({
      ht3 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show, pattern = "outgoing")
      print(ht3)
    })
    try({
      ht4 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show, pattern = "incoming")
      print(ht4)
    })
    rm(ht3,ht4)
  }
  dev.off()
  rm(i,ht1,ht2)

  ##### Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together #####
  ##### Identify and visualize outgoing communication pattern of secreting cells #####
  library(NMF)
  library(ggalluvial)
  P.outgoing <- selectK(cellchat, pattern = "outgoing")
  P.incoming <- selectK(cellchat, pattern = "incoming")

  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_GlobalPatterns_outgoing.pdf"),
      width = 12,  height = 8
  )

  nPatterns = 4
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  # P.nHeatmap.outgoing <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  # P.nHeatmap.outgoing
  # river plot
  netAnalysis_river(cellchat, pattern = "outgoing")
  #> Please make sure you have load `library(ggalluvial)` when running this function

  # dot plot
  netAnalysis_dot(cellchat, pattern = "outgoing")

  P.outgoing
  #graphics.off()
  dev.off()
  #### Identify and visualize incoming communication pattern of target cells ####

  pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_GlobalPatterns_incoming.pdf"),
      width = 12,  height = 8
  )


  nPatterns = 4
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

  # river plot
  netAnalysis_river(cellchat, pattern = "incoming")
  #> Please make sure you have load `library(ggalluvial)` when running this function

  # dot plot
  netAnalysis_dot(cellchat, pattern = "incoming")

  P.incoming
  #graphics.off()
  dev.off()

  # ##### Manifold and classification learning analysis of signaling networks #####
  #   #### Identify signaling groups based on their functional similarity ####
  #     cellchat <- computeNetSimilarity(cellchat, type = "functional")
  #     cellchat <- netEmbedding(cellchat, type = "functional")
  #     #> Manifold learning of the signaling networks for a single dataset
  #     cellchat <- netClustering(cellchat, type = "functional")
  #     #> Classification learning of the signaling networks for a single dataset
  #     # Visualization in 2D-space
  #     netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  #
  #   #### Identify signaling groups based on structure similarity  ####
  #     cellchat <- computeNetSimilarity(cellchat, type = "structural")
  #     cellchat <- netEmbedding(cellchat, type = "structural")
  #     #> Manifold learning of the signaling networks for a single dataset
  #     cellchat <- netClustering(cellchat, type = "structural")
  #     #> Classification learning of the signaling networks for a single dataset
  #     # Visualization in 2D-space
  #     netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  #
  #     netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

##### Part V: Save the CellChat object #####
  saveRDS(cellchat, file = paste0(Save_Path,"cellchat_humanSkin_LS.rds"))

    ##### Save CellChatDataBase #####
    PathDB <- paste0(Save_Path,"/DataBase")
    ## Create new folder
    if (!dir.exists(PathDB)){
      dir.create(PathDB)
    }

    ## Export all database
    DB_Interact_All.df <- data.frame(Term = row.names(CellChatDB[["interaction"]]), CellChatDB[["interaction"]])
    write.table(DB_Interact_All.df,
                file=paste0(PathDB,"/",ProjectName,"_DBAll_Interact.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Complex_All.df <- data.frame(Term = row.names(CellChatDB[["complex"]]), CellChatDB[["complex"]])
    write.table(DB_Complex_All.df,
                file=paste0(PathDB,"/",ProjectName,"_DBAll_Complex.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Cofactor_All.df <- data.frame(Term = row.names(CellChatDB[["cofactor"]]), CellChatDB[["cofactor"]])
    write.table(DB_Cofactor_All.df,
                file=paste0(PathDB,"/",ProjectName,"_DBAll_Cofactor.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_GeneInfo_All.df <- data.frame(Term = row.names(CellChatDB[["geneInfo"]]), CellChatDB[["geneInfo"]])
    write.table(DB_GeneInfo_All.df,
                file=paste0(PathDB,"/",ProjectName,"_DBAll_GeneInfo.tsv"),sep="\t",
                row.names=F, quote = FALSE)

    ## Export used database
    DB_Interact.df <- data.frame(Term = row.names(CellChatDB.use[["interaction"]]), CellChatDB.use[["interaction"]])
    write.table(DB_Interact.df,
                file=paste0(PathDB,"/",ProjectName,"_DB_Interact.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Complex.df <- data.frame(Term = row.names(CellChatDB.use[["complex"]]), CellChatDB.use[["complex"]])
    write.table(DB_Complex.df,
                file=paste0(PathDB,"/",ProjectName,"_DB_Complex.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Cofactor.df <- data.frame(Term = row.names(CellChatDB.use[["cofactor"]]),CellChatDB.use[["cofactor"]])
    write.table(DB_Cofactor.df,
                file=paste0(PathDB,"/",ProjectName,"_DB_Cofactor.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_GeneInfo.df <- data.frame(Term = row.names(CellChatDB.use[["geneInfo"]]),CellChatDB.use[["geneInfo"]])
    write.table(DB_GeneInfo.df,
                file=paste0(PathDB,"/",ProjectName,"_DB_GeneInfo.tsv"),sep="\t",
                row.names=F, quote = FALSE)

    #### Catch the significant path ####
    DB_Interact_Sig.df <- DB_Interact.df[DB_Interact.df$pathway_name %in% pathway.set,]
    DB_Interact_Sig.df$pathway_name %>% unique()

    write.table(DB_Interact.df,
                file=paste0(PathDB,"/",ProjectName,"_DBSig_Interact.tsv"),sep="\t",
                row.names=F, quote = FALSE)


    #### Automatically save the plots of the all inferred network for quick exploration ####
    # # Access all the signaling pathways showing significant communications
    # pathways.show.all <- cellchat@netP$pathways
    # # check the order of cell identity to set suitable vertex.receiver
    # levels(cellchat@idents)
    # vertex.receiver = seq(1,4)
    # for (i in 1:length(pathways.show.all)) {
    #   # Visualize communication network associated with both signaling pathway and individual L-R pairs
    #   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
    #   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    #   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
    #   ggsave(filename=paste0(Version,"/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
    # }

    #### Save the RData ####
    save.image(paste0(Version,"_CellChat_Example_PRJCA001063.RData"))

    CellChat.lt <- list(CellChatObj = cellchat,
                        DataBase_Use = CellChatDB.use,
                        DB_Interact_Sig = DB_Interact_Sig.df
    )

    return(CellChat.lt)
}
