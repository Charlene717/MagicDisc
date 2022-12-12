##### Verification (CellCheck) #####
#### Install ####
## Check whether the installation of those packages is required
Package.set <- c("tidyverse","caret","cvms","DescTools","devtools")
for (i in 1:length(Package.set)) {
  if (!requireNamespace(Package.set[i], quietly = TRUE)){
    install.packages(Package.set[i])
  }
}
## Load Packages
# library(Seurat)
lapply(Package.set, library, character.only = TRUE)
rm(Package.set,i)

## install CellCheck
# Install the CellCheck package
detach("package:CellCheck", unload = TRUE)
devtools::install_github("Charlene717/CellCheck")
# Load CellCheck
library(CellCheck)

#### Run CellCheck ####
## Create check dataframe
scSorter.df <- data.frame(Actual = scRNA.SeuObj_Small@meta.data[["celltype"]],
                          Predict = scRNA.SeuObj_Small@meta.data[["scSorterPred"]])
scSorter_Anno.df <- data.frame(TestID = "Predict",
                               Tool = "scSorter",
                               Type = "PDAC",
                               PARM = "1")
## For one prediction
DisMultCM.lt <- list(Actual = "Actual", Predict = "Predict", Type = "Type", Type2 = "PDAC" )
cm_DisMult.lt <- CellCheck_DisMult(scSorter.df, scSorter_Anno.df, Mode = "One", DisMultCM.lt,
                                   Save.Path = Save.Path, ProjectName = ProjectName)
## For multiple prediction
Sum_DisMult.df <- CellCheck_DisMult(scSorter.df, scSorter_Anno.df,
                                    Mode = "Multiple",DisMultCM.lt=DisMultCM.lt,
                                    Save.Path = Save.Path, ProjectName = ProjectName)
