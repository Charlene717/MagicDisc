FUN_GSEA_MultiCell <- function(CCMarker.lt, CellType.list,
                               Sampletype = sampletype,
                               ProjectName = projectName
                               ) {

  GSEAFilname <- paste0(sampletype,"_",projectName,"_GSEA")
  PathGSEA <- paste0(paste0(Save.Path,"/",GSEAFilname))
  ## Create new folder
  if (!dir.exists(PathGSEA)){
    dir.create(PathGSEA)
  }


  GSEA_Large <- list()
  GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large.df.TOP <- GSEA_Large.df



  pdf(file = paste0(PathGSEA,"/",GSEAFilname,".pdf"),width = 15, height = 7 )

  for(i in 1:length(CellType.list)){

    gseaDat <- CCMarker.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
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
                  gseaParam = 0.5) + title( paste0(sampletype,".",projectName,".",CellType.list[i]), adj = 0, line =3)

    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+
      labs(title= paste0(sampletype,".",projectName,".",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+
      labs(title= paste0(sampletype,".",projectName,".",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
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
  write.table(GSEA_Large.Sum.TOP, file=paste0(PathGSEA,"/",GSEAFilname,"_Pathway_LargeTOP_SPA.txt"),sep="\t",
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

  pdf(file = paste0(PathGSEA,"/",GSEAFilname,"_Bubble_SPA.pdf"),width = 17, height = 12 )
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


  pdf(file = paste0(PathGSEA,"/",GSEAFilname,"_Bubble_SPA_SubType_T.pdf"),width = 17, height = 7 )
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

  pdf(file = paste0(PathGSEA,"/",GSEAFilname,"_Bubble_SPA_SubType_Mac.pdf"),width = 17, height = 20 )
  BBPlot_MacB
  BBPlot_MacB1
  dev.off()

  rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
     df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)
  return(OUTPUT)
}
