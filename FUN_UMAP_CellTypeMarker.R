# --------------- Check specific tissue marker --------------- #
UMAP_CellTypeMarker <- function(){
    pdf(
      file = paste0(PathCellType,"/",ProjectName,"_nlDR_CTMarker.pdf"),
      width = 10,  height = 8
    )
    
      # PMID: 31771616 #!!!!!!
      FeaturePlot(scRNA.SeuObj, features = c("Cd3d", "Cd4", "Cd8a", "Csf1r", "Foxp3", "S100a9"), min.cutoff = "q9",
                  ncol = 3, coord.fixed = 1)
      # T Cell: Cd3d;  CD4+ T Cell: Cd4; CD8+ T Cell: Cd8a; Macrophages: Csf1r; regulatory T cells(Treg): Foxp3; Neutrophils: S100a9
      
      # PMID: 34296197 #!!!!!
      FeaturePlot(scRNA.SeuObj, features = c("Cd3d", "Cd3e", "Lyz1", "Lyz2","Clu","Cd79a","Ms4a1","Nkg7","Gzmb"), min.cutoff = "q9",
                  ncol = 3, coord.fixed = 1)
      # T Cell: Cd3d,Cd3e;  Macrophages: Lyz; Mast Cell: Clu; B Cell: Cd79a,Ms4a1; NK Cell: Nkg7,Gzmb
      
      # http://biocc.hrbmu.edu.cn/CellMarker/
      # Mast cell
      FeaturePlot(scRNA.SeuObj, features = c("Cd117", "Cd25","Cd203c","Slc18a2","Kit","Fcer1a","Cd9"), min.cutoff = "q9", coord.fixed = 1)
      # PMID: 30356731
      FeaturePlot(scRNA.SeuObj, features = c("Cd9"), min.cutoff = "q9", coord.fixed = 1)
      # PMID: 34296197 #!!!!!
      FeaturePlot(scRNA.SeuObj, features = c("Cpa3"), min.cutoff = "q9", coord.fixed = 1)
      # https://www.panglaodb.se/markers.html?cell_type=%27Mast%20cells%27
      
      ## Mac
      # Macrophage-Markers
      # https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
      ## M0
      FeaturePlot(scRNA.SeuObj, features = c("Cd68", "Adgre1","Cd14","Csf1r","Ly6c1",
                                             "Cx3cr1","Fcgr1a","Itgam","Mertk"), min.cutoff = "q9", coord.fixed = 1)
      
      ## M1
      # M1 http://biocc.hrbmu.edu.cn/CellMarker/
      FeaturePlot(scRNA.SeuObj, features = c("Cd16","Cd32","Cd64","Cd68","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1)
      # M1 https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
      FeaturePlot(scRNA.SeuObj, features = c("Marco","Nos2","Tlr2","Cd80","Cd86","Csf2",
                                             "Tnf","Il1b","Il6","Tlr4","Cxcl2","Ifng","Il1r1"), min.cutoff = "q9", coord.fixed = 1)
      FeaturePlot(scRNA.SeuObj, features = c("Il1a","Il1b","Il6","Nos2","Tlr2","Tlr4","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1)
      
      
      ## M2
      # M2 http://biocc.hrbmu.edu.cn/CellMarker/
      FeaturePlot(scRNA.SeuObj, features = c("Chil3","Csf1r","Mrc1","Pparg","Arg1","Cd163","Clec10a","Clec7a",
                                             "Cd206","Cd209","Ccl18","Fizz1"), min.cutoff = "q9", coord.fixed = 1)
      # https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
      FeaturePlot(scRNA.SeuObj, features = c("Cd115", "Cd206", "Pparg", "Arg1", "Cd163", "Cd301",
                                             "Dectin-1", "Pdcd1lg2", "Fizz1"), min.cutoff = "q9", coord.fixed = 1)
      
      FeaturePlot(scRNA.SeuObj, features = c("Chil3"), min.cutoff = "q9", coord.fixed = 1)
      
      
      
      ## Tumor associated macrophage(TAM)
      FeaturePlot(scRNA.SeuObj, features = c("Ccr2","Csf1r","Marco","Pdl2","Cd40","Ccl2","Csf1","Cd16"),
                  min.cutoff = "q9", coord.fixed = 1)
      
      # Erythrocytes
      FeaturePlot(scRNA.SeuObj, features = c("Hbb-bs"), min.cutoff = "q9", coord.fixed = 1)
      # Platelet
      FeaturePlot(scRNA.SeuObj, features = c("Ppbp"), min.cutoff = "q9", coord.fixed = 1)
      
      ## Summary
      markers.to.plot <- c("Cd3d","Cd3e", "Cd4","Cd8a", "Csf1r", "Lyz2","Chil3","Il1b", "S100a9","Nkg7",
                           "Gzmb", "Cd79a", "Ms4a1","Clu","Hbb-bs","Ppbp")
      
      FeaturePlot(scRNA.SeuObj, features = markers.to.plot, min.cutoff = "q9", coord.fixed = 1)
      # T Cell: Cd3d,Cd3e;  CD4+ T Cell: Cd4; CD8+ T Cell: Cd8a; Macrophages: Csf1r,Lyz,Chil3;  Neutrophils: S100a9;
      # NK Cell: Nkg7,Gzmb; B Cell: Cd79a,Ms4a1; Mast Cell: Clu; Erythrocytes: Hbb-bs; Platelet: Ppbp
    
    
    dev.off()
}