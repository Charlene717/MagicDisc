##### 06 Cell type annotation*  #####
##### 06 Cell type annotation: Traditional method (Manual annotation by cell type marker)  #####
# scRNA.SeuObj.copy <- scRNA.SeuObj

## CD4+T: CD4+T Cell; CD8+T: CD8+T Cell; T: T Cell; B: B Cell; Mac: Macrophages;
## Neu: Neutrophils; NK: NK Cell; Mast: Mast Cell; Ery: Erythrocytes;
## Thr: Thrombocytes
scRNA.SeuObj <- RenameIdents(scRNA.SeuObj, `0` = "CD4+T", `1` = "B", `2` = "Mac3",
                             `3` = "Neu", `4` = "CD8+T", `5` = "CD8+T", `6` = "Mac2", `7` = "CD4+T",
                             `8` = "NK", `9` = "Neu",`10` = "Mast1", `11` = "T", `12` = "Ery", `13` = "Mac1",
                             `14` = "B", `15` = "B", `16` = "Mast2", `17` = "Mac0", `18` = "Neu")

Cell_Type_Order.set <- c("T", "CD4+T", "CD8+T", "B" , "Mac0", "Mac1", "Mac2", "Mac3", "Mast1", "Mast2", "NK", "Neu", "Ery")

scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
# Idents(scRNA.SeuObj) <- "celltype"

## Heatmap
Heatmap_Color.lt <- list(low="#5283ff",mid ="white", high ="#ff5c5c")
scRNA.SeuObj <- Beautify_Heatmap_Seurat(scRNA.SeuObj, PBMC.markers.df, topN = 7, Path = PathCellType,
                                        Type = "celltype", HMColor.lt = Heatmap_Color.lt,
                                        projectName = ProjectName)

## Export All DRPlot(UMAP,tSNE)
Export_All_DRPlot(scRNA.SeuObj)

## Summary
markers.to.plot <- c("Cd3d","Cd3e", "Cd4","Cd8a", "Csf1r", "Lyz2","Chil3","Il1b", "S100a9","Nkg7",
                     "Gzmb", "Cd79a", "Ms4a1","Clu","Hbb-bs","Ppbp")
## DotPlot
#(Pending)
BeautifyDotPlot(scRNA.SeuObj, Path = PathCellType, projectName = ProjectName,
                Features = markers.to.plot)

## Create cell type markers dataframe
# DefaultAssay(scRNA.SeuObj_Small) <- "RNA"
scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
CellType.markers.df <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

write.table(CellType.markers.df, file = paste0(PathCellType,"/CC_CelltypeMarker_AllGene.txt"),
            quote = F,sep = "\t",row.names = F)

# #### Save RData ####
# save.image(paste0(Save.Path,"/06_Cell_type_annotation.RData"))
