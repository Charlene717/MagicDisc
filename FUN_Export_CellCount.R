ExportCellCount <- function(AnnoSummary.lt,Path = Save.Path,
                            Folder = "B01_CellCount", projectName = ProjectName) {
  # Create new folder
  NewPath <- paste0(Path,"/",Folder)
  dir.create(NewPath)

  # Extract annotation dataframe
  Anno.df <- AnnoSummary.lt[["Anno.df"]]

  #### Export table ####
  write.table( AnnoSummary.lt[["Freq_All_Cla.lt"]][["Freq_All_Cla.df"]] ,
               file = paste0(NewPath,"/",projectName,"_CellCount_CT_Cla.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )


  #### LinePlot ####
  # https://ithelp.ithome.com.tw/articles/10186047
  # Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
  #                                    levels = sort(unique(as.character(Freq_All.df$Cell_Type))))

  Freq_All.df <- AnnoSummary.lt[["Freq_All.df"]]
  Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
                                  levels = Cell_Type_Order.set)

  ## Plot by State
  CellNum_P1 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Number,
                                        colour = Pheno_Type, group = Pheno_Type)) +
    geom_line(linetype = "dashed",size=1.5) +
    geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  CellNum_P1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=15,  YtextSize=15, xangle = 90,
                                LegTextSize = 15) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P1
  CellNum_P1

  CellNum_P2 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Percent,
                                        colour = Pheno_Type, group = Pheno_Type)) +
    geom_line(linetype = "dashed",size=1.5) +
    geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  CellNum_P2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=18, xangle = 90,
                                LegTextSize = 15) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P2
  CellNum_P2

  ## Combine sex
  Freq_All_Cla.df <- AnnoSummary.lt[["Freq_All_Cla.lt"]][["Freq_All_Cla.df"]]

  # Freq_All_Cla.df$Cell_Type <- factor(Freq_All_Cla.df$Cell_Type,
  #                                    levels = sort(unique(as.character(Freq_All_Cla.df$Cell_Type))))
  Freq_All_Cla.df$Cell_Type <- factor(Freq_All_Cla.df$Cell_Type,
                                      levels = Cell_Type_Order.set)

  CellNum_P3 <- ggplot(Freq_All_Cla.df, aes(x = factor(Cell_Type), y = Number,
                                            colour = Pheno_Type,
                                            group = Pheno_Type,linetype=Pheno_Type)) +
    geom_line(aes(linetype=Pheno_Type))+
    geom_line(size=1.5) +
    scale_linetype_manual(name="Pheno_Type",
                          values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                          #labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")
    ) +
    scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
    geom_point(shape = 12, size = 4, fill = "white") + theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


  CellNum_P3 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.4, 0.8),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=18,xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P3
  CellNum_P3

  # library(eoffice)
  # topptx(CellNum_P3,paste0(NewPath,"/Temp.pptx"))

  CellNum_P4 <- ggplot(Freq_All_Cla.df, aes(x = factor(Cell_Type), y = Percent,
                                            colour = Pheno_Type,
                                            group = Pheno_Type,linetype=Pheno_Type)) +
    geom_line(aes(linetype=Pheno_Type))+
    geom_line(size=1.5) +
    scale_linetype_manual(name="Pheno_Type",
                          values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                          #labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")
    ) +
    scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
    geom_point(shape = 12, size = 4, fill = "white") +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
    #theme_set(theme_bw())+ # Remove the background
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  CellNum_P4 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.15, 0.82),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P4
  CellNum_P4

  rm(Freq_All_Cla.df)


  #### BarPlot ####
  # https://blog.gtwang.org/r/ggplot2-tutorial-layer-by-layer-plotting/3/
  colnames(Anno.df)[ncol(Anno.df)] <- "Cell_Type"
  Anno.df$Cell_Type <- factor(Anno.df$Cell_Type,
                              levels = sort(unique(as.character(Anno.df$Cell_Type))))

  # sample
  BarPlot1_1 <- ggplot(Anno.df, aes(Cell_Type, fill=Anno.df[,1])) +
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
    labs(fill=colnames(Anno.df)[1]) -> BarPlot1_1
  BarPlot1_1

  #### Export PDF file ####
  pdf(file = paste0(NewPath,"/",projectName,"_CellCount_LinePlot.pdf"),
      width = 7, height = 7 )
  CellNum_P4 %>% print()
  CellNum_P3 %>% print()
  CellNum_P1 %>% print()
  CellNum_P2 %>% print()
  for (i in 1:(ncol(Anno.df)-1)) {
    # sample
    BarPlot1_1 <- ggplot(Anno.df, aes(Cell_Type, fill = Anno.df[,i])) +
      geom_bar(position="dodge")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
      labs(fill=colnames(Anno.df)[i]) -> BarPlot1_1
    print(BarPlot1_1)

    BarPlot1_2 <- ggplot(Anno.df, aes(Cell_Type, fill = Anno.df[,i])) +
      geom_bar(position="fill")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot1_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) +
      labs(fill=colnames(Anno.df)[i])+
      labs(y="Proportion") -> BarPlot1_2
    print(BarPlot1_2)
    rm(BarPlot1_1,BarPlot1_2)

  }
  rm(i)

  dev.off() # graphics.off()

  rm(CellNum_P1,CellNum_P2,CellNum_P3,CellNum_P4,Anno.df)

  # return(OUTPUT)
}
