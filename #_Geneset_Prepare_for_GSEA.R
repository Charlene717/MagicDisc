##### 09_0 GSEA Analysis (Geneset Prepare) #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  library(tidyverse)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_GSEA_ggplot.R")


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


  InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"


  ## Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  Pathway.all <- read.delim2(paste0(getwd(),"/",InputGSEA),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",InputGSEA))),
                             header = F,sep = "\t")

##### GSEA_SaveMMRData.R #####
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


    #### Save RData ####
    SaveGSEA.Path <- Save.Path
    rm(list=setdiff(ls(),c("Pathway.all.MM","SaveGSEA.Path")))
    save.image(paste0(SaveGSEA.Path,"/GSEA_Analysis_Geneset.RData"))
    rm(list=setdiff(ls(),c("Pathway.all.MM")))
    save.image(paste0("GSEA_Analysis_Geneset.RData"))
