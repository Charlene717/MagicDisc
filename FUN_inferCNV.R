# rm(list = ls()) #?M???ܼ?


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
  EM.mt <-  scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame() %>%
            dplyr::filter(., rowSums(.) > 0, .preserve = F) %>%
            as.matrix()

  ## Need to create mouse gene_order_file
  ## Ref: https://zhuanlan.zhihu.com/p/111562837
  rownames(EM.mt) <- rownames(EM.mt) %>% toupper()

  ## Create annotaion matrix
  Anno.mt <- scRNA.SeuObj@meta.data %>%
             dplyr::select("celltype") %>%
             as.matrix()

##### inferCNV #####
  ## create the infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = EM.mt,
                                      annotations_file = Anno.mt,
                                      #        delim="\t",
                                      gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                      # ref_group_names=c("AC","nAtD","ND01"))
                                      ref_group_names=c("T","B"))

  ## Run inferCNV
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               #  out_dir= "output_dir",
                               out_dir=tempfile(),
                               cluster_by_groups=TRUE,
                               plot_steps=FALSE,
                               no_plot=FALSE,
                               denoise=TRUE,
                               resume_mode = FALSE,
                               HMM=TRUE)





