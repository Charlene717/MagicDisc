CombineSeuObj <- function(scRNA_SeuObj.list,selection.method = "vst", nfeatures = 2000 ,...) {

  if(length(scRNA_SeuObj.list)==1){
    scRNA.SeuObj <- scRNA_SeuObj.list[[1]]
  }else{
    # normalize and identify variable features for each dataset independently
    set.seed(1) # Fix the seed
    scRNA_SeuObj.list <- lapply(X = scRNA_SeuObj.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = selection.method, nfeatures = nfeatures)
    })

    # select features that are repeatedly variable across datasets for integration
    set.seed(1) # Fix the seed
    features <- SelectIntegrationFeatures(object.list = scRNA_SeuObj.list, nfeatures = nfeatures)

    ## Perform integration
    set.seed(1) # Fix the seed
    scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA_SeuObj.list,
                                            anchor.features = features)
    # this command creates an 'integrated' data assay
    set.seed(1) # Fix the seed
    scRNA.SeuObj <- IntegrateData(anchorset = scRNA.anchors,...)

    set.seed(1) # Fix the seed
    DefaultAssay(scRNA.SeuObj) <- "integrated"

  }

  return(scRNA.SeuObj)
}
