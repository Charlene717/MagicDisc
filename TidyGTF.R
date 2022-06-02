
## Ref: https://www.biostars.org/p/272889/
## Ref: https://bioconductor.org/packages/release/bioc/html/rtracklayer.html
GenecodePath <- "./Input_files/Genecode"
GenecodeFilename <- "gencode.vM29.annotation"

gtf <- rtracklayer::import(paste0(GenecodePath,"/",GenecodeFilename,".gtf"))
gtf_df = as.data.frame(gtf) %>%
            dplyr::filter(type=="gene") %>%
            relocate(gene_name,.before = seqnames) %>%
            distinct(.[,1:4]) %>%
            distinct(.[,1], .keep_all = TRUE)
gtf_df <- gtf_df[,-5]


write.table(gtf_df, paste0(GenecodePath,"/",GenecodeFilename,".txt"),
            sep = "\t", row.names = F,quote = FALSE, col.names = F)

# TTT <- read.delim(paste0(getwd(),"/Input_files/Genecode/gencode.vM29.annotation.txt"),header  = F)
# TTT$V1 %>% unique() %>% as.data.frame() -> TTT2
