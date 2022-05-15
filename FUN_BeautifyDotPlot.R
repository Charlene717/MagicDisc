BeautifyDotPlot <- function(scRNA.SeuObj, Path = PathCellType,
                            projectName = ProjectName, Features = markers.to.plot){

  ## DotPlot
  DotPlot_Color1.set <- c("#de3767", "#de3767", "#4169e1", "#4169e1")
  DotPlot_Color2.set <- c("#5b8e7d","#7b2cbf")
  DotPlot_Color3.set <- c("#de3767", "#4169e1")
  
  pdf(
    file = paste0(Path,"/",projectName,"_DotPlot_CellType",".pdf"),
    width = 10,  height = 8
  )
  
  # https://satijalab.org/seurat/reference/dotplot
  p1 <- DotPlot(scRNA.SeuObj, features = Features, cols = c("lightgrey", "blue"),
                dot.scale = 8) + RotatedAxis()%>%
                BeautifyggPlot(.,LegPos = "bottom",AxisTitleSize=1, TitleSize = 20, xangle =90,
                LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 8, XaThick=1, YaThick=1,XtextSize=12,  YtextSize=12)
  print(p1)
  
  p2 <- DotPlot(scRNA.SeuObj, features = Features, cols = DotPlot_Color1.set,
                dot.scale = 8, split.by = "Sample") + RotatedAxis()
  print(p2)
  
  # https://github.com/satijalab/seurat/issues/1541
  p3 <- DotPlot(scRNA.SeuObj, features = Features, cols = DotPlot_Color2.set,
                dot.scale = 8, split.by = "Cachexia") + RotatedAxis()
  print(p3)
  
  p4 <- DotPlot(scRNA.SeuObj, features = Features, cols = DotPlot_Color3.set,
                dot.scale = 8, split.by = "Sex") + RotatedAxis()
  print(p4)
  
  dev.off()
  
return()

}