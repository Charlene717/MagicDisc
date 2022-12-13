## To-do list ##
# - [x] Basic setting
# - [x] Mouse & Human
#   - [x] Basic setting
#   - [x] genecode file setting
#   - [x] create new version of genecode file
# - [x] Parameter setting
#   - [x] Basic setting
#   - [x] CreateInfercnvObject
#   - [x] infercnv::run
# - [x] RefGroup setting

## Add to-do list in advance ##
# - [ ] Beautify graphics
# - [ ] Plot connections in line or plot Chord Diagram
        # Ref: https://cole-trapnell-lab.github.io/cicero-release/docs_m3/
# - [ ] Extract genes of interest for statistical analysis and visualization
# - [ ] Integrated multi-omics analysis of CNV, patient status and other data with clinical databases such as TCGA (bulk data)

inferCNV <- function(scRNA.SeuObj, AnnoSet = "celltype",
                     SpeciSet = Species,
                     Path = "", infercnvCutOff = 0.1, denoiseSet = TRUE, HMMSet = TRUE,
                     GenecodeSet.list = list(Default = TRUE,
                                        HumanGenecode = paste0(getwd(),"/Input_files/Genecode/gencode.v40.annotation.txt"), # "/Input_files/Genecode/gencode_v19_gene_pos.txt"
                                        MouseGenecode = paste0(getwd(),"/Input_files/Genecode/gencode.vM29.annotation.txt")),
                     RefSet = c("normal"),
                     CreateInfercnvObject.lt = "",
                     inferCNVRun.lt = ""
                     ) {

  memory.limit(150000)

  ##### Load package #####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("infercnv")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


  ##### Data preprocessing #####
  ## Create expression matrix

  # # ## Old version (Without normalizaiton) ##
  # EM.mt <-  scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame() %>%
  #           dplyr::filter(., rowSums(.) > 0, .preserve = F) %>%
  #           as.matrix()

  if(!require("Seurat")) install.packages("Seurat")
  library(Seurat)
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)

  EM.mt <-  GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() %>%
            dplyr::filter(., rowSums(.) > 0, .preserve = F) %>%
            as.matrix()

  # ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  # GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix



  # ## SpeciSet
  # if(SpeciSet == "Mouse"){
  #   ## Need to create mouse gene_order_file
  #   ## Ref: https://zhuanlan.zhihu.com/p/111562837
  #   rownames(EM.mt) <- rownames(EM.mt) %>% toupper()
  # }

  ## Create annotaion matrix
  Anno.mt <- scRNA.SeuObj@meta.data %>%
             dplyr::select(AnnoSet) %>%
             as.matrix()

  ##### inferCNV #####
  ## May need to create new version of gene_order_file by yourself
  ## Ref: https://www.jieandze1314.com/post/cnposts/206/
  ## Ref: https://github.com/broadinstitute/infercnv/blob/master/scripts/gtf_to_position_file.py
  ## Ref(OldVersion): https://data.broadinstitute.org/Trinity/CTAT/cnv/

  if(GenecodeSet.list[["Default"]]== TRUE){
    if(SpeciSet == "Mouse"){
      GenecodePath <- paste0(getwd(),"/Input_files/Genecode/gencode.vM29.annotation.txt")
    }else{
      GenecodePath <- paste0(getwd(),"/Input_files/Genecode/gencode.v40.annotation.txt") # "/Input_files/Genecode/gencode_v19_gene_pos.txt"
    }
  }else{
    if(SpeciSet == "Mouse"){
      GenecodePath <- GenecodePath.list[["MouseGenecode"]]
    }else{
      GenecodePath <- GenecodePath.list[["HumanGenecode"]]
    }
  }


  #### Set reference cell ####
  for (i in 1:length(RefSet)) {
    if(i==1){
      RefGroup <- dplyr::filter(as.data.frame(Anno.mt), grepl(RefSet[i] , Anno.mt)) %>%
        unique() %>% unlist()
    }
    RefGroup_Temp <- dplyr::filter(as.data.frame(Anno.mt), grepl(RefSet[i] , Anno.mt)) %>%
      unique() %>% unlist()
    RefGroup <- c(RefGroup,RefGroup_Temp) %>% unique()
  }
  rm(i,RefGroup_Temp)

  #### create the infercnv object ####
  if(CreateInfercnvObject.lt == ""){
    CreateInfercnvObject.lt = list(delim="\t",max_cells_per_group = NULL,min_max_counts_per_cell = c(100, +Inf),chr_exclude = c("chrX", "chrY", "chrM"))
  }

  formals(CreateInfercnvObject)[names(CreateInfercnvObject.lt)] <- CreateInfercnvObject.lt
  CreateInfercnvObject.lt = CreateInfercnvObject.lt
  formals(CreateInfercnvObject)[names(CreateInfercnvObject.lt)] <- CreateInfercnvObject.lt
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = EM.mt,
                                      annotations_file = Anno.mt,
                                      gene_order_file = GenecodePath,
                                      # gene_order_file = system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                      ref_group_names = RefGroup)


  #### Run inferCNV ####

  inferCNVRun.lt = list(cluster_by_groups=TRUE, plot_steps=FALSE, no_plot=FALSE, resume_mode = FALSE, k_nn = 30)
  formals(run)[names(inferCNVRun.lt)] <- inferCNVRun.lt
  inferCNVRun.lt = inferCNVRun.lt
  formals(run)[names(inferCNVRun.lt)] <- inferCNVRun.lt
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff = infercnvCutOff,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= paste0(Path), # out_dir=tempfile(),  #  out_dir= "output_dir",
                               denoise=denoiseSet,
                               HMM=HMMSet)

  return(infercnv_obj)

}


## inferCNV parameter
## Ref: https://rdrr.io/github/broadinstitute/inferCNV/man/run.html
## Ref: https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
