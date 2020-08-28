smpls <- list.files(".","_R_")

h.smpls <- grep("_C[0-9]+_",smpls,value=T)
d.smpls <- grep("_U[0-9]+_",smpls,value=T)

x.smpls <- d.smpls

require(Seurat)
seu.lst <- vector("list",length(x.smpls))
names(seu.lst) <- x.smpls
gene.names <- vector("list",length(x.smpls))
for(i in 1:length(x.smpls)){
    print(x.smpls[i])
    cnts  <- read.table(x.smpls[i],header=T,sep='\t',check.names=F,stringsAsFactors=F)
    print(cnts[1:2,1:2])
    rownames(cnts) <- cnts[,1]
    cnts <- cnts[,-1]
    cnts <- t(cnts)
    gene.names[[i]] <- rownames(cnts)
    print(dim(cnts))
    seu.lst[[i]] <- CreateSeuratObject(cnts,min.cells=3,min.features=200)
    seu.lst[[i]] <- subset(seu.lst[[i]], subset = nFeature_RNA > 400)
    seu.lst[[i]] <- NormalizeData(seu.lst[[i]])
    seu.lst[[i]] <- FindVariableFeatures(seu.lst[[i]])
}
for(i in 1:(length(x.smpls)-1))
    print(all(gene.names[[i]] %in% gene.names[[i+1]]))

saveRDS(seu.lst,"uc_smpls_seu_list.rds")

anchors <- FindIntegrationAnchors(seu.lst)
saveRDS(anchors,"uc_smpls_anchors.rds")

anchors <- readRDS("uc_smpls_anchors.rds")
seu.itg <- IntegrateData(anchors)

seu.itg <- ScaleData(seu.itg)
seu.itg <- RunPCA(seu.itg,npcs=100)

## seu.itg <- JackStraw(seu.itg, num.replicate=100, dims=100)
## seu.itg <- ScoreJackStraw(seu.itg, dims = 1:100)
## JackStrawPlot(seu.itg, dims = 51:70)
saveRDS(seu.itg,"uc_integrated_Seu.rds")

ElbowPlot(seu.itg, ndims=50)

seu.itg <- RunUMAP(seu.itg, dims = 1:25)
seu.itg <- RunTSNE(seu.itg, dims = 1:25)
seu.itg <- FindNeighbors(seu.itg, dims=1:25)
seu.itg <- FindClusters(seu.itg, resolution=seq(0.4,2.0,0.4))
seu.itg <- FindClusters(seu.itg, resolution=seq(0.1,2.0,0.4))
library(clustree)
clustree(seu.itg)
## For Ctrl sample resolution = 0.9 based on clustree analysis (ctrl_clustree.png)
## For UC sample resolution = 0.5 based on clustree analysis (uc_clustree.png)
saveRDS(seu.itg,"uc_integrated_Seu_clustered.rds")

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Cell-type annotations for Control (manual, canonical clusters)
### =================================================================
seu.itg = readRDS("ctrl_integrated_Seu_clustered.rds")
DefaultAssay(seu.itg) = "RNA"
## Annotate larger groups of cells
new.cluster.names <- vector("character",length(unique(Idents(seu.itg))))
names(new.cluster.names) <- 0:(length(unique(Idents(seu.itg)))-1)

## =======
## T-cells
FeaturePlot(seu.itg,features=c("CD3D","CD3E","CD3G"))
##=> Clusters 0,1,3,6,8,9,11 are T-cells

## CD8+ T-cells (cytotoxic) and CD4+ ones
FeaturePlot(seu.itg,features=c("CD8A","CD8B","CD4"))
##=> Clusters 3, 11 CD8+ T-cells
new.cluster.names[c("3","11")] <- "CD8.T.cell"

## Tregs: FOXP3+
FeaturePlot(seu.itg,features=c("CD8A","CD8B","CD4","FOXP3"))
new.cluster.names["9"] <- "Treg.cell"

##=> Clusters 0,1,6,8 CD4+ T cells
new.cluster.names[c("0","1","6","8")] <- "CD4.T.cell"

## ===================
## NK: CD3-NCAM+FCGR3+
FeaturePlot(seu.itg,features=c("CD3D","CD3E","CD3G","NCAM1","FCGR3A"))
## No

## ===================
## B-cells
## Plasma cells (final diff stage of B): CD79A+SDC1+CD78+(gene name??)IL6R+, CD19-CD20-/MS4A1-
FeaturePlot(seu.itg,features=c("CD79A","CD19","MS4A1","SDC1","IL6R"))
##=> Clusters 2,5,7,10,12,13,15 Plasma
new.cluster.names[c("2","5","7","10","12","13","15")] <- "Plasma.cell"
##=> Clusters 4,14 B cells
new.cluster.names[c("4","14")] <- "B.cell"

## =====================
## Macrophages/Monocytes
FeaturePlot(seu.itg,features=c("CD14","FCGR3A"))
##=> Cluster 16
new.cluster.names["16"] <- "Mac"

## FindMarkers(seu.itg,ident.1=17,ident.2=c(0,1,6,8,9,3,11),only.pos=T, min.pct=0.25,logfc.threshold=0.6)
## Identifies Cluster 17 as DC.cell from PanglaoDB by markers:
## IL4I1,LST1,IL23R,AFF3,PTGDR,CD83,KLF4,LTB
new.cluster.names["17"] <- "DC.cell"

new.cluster.names[c("18","19")] <- "Unknown"

## =================================================
seu.itg = RenameIdents(seu.itg, new.cluster.names)
seu.itg$Cell.types = Idents(seu.itg)

write.table(seu.itg[["RNA"]]@data,"ctrl_R_integrated_dat.txt",quote=F,row.names=T,col.names=T)
anno.tbl <- data.frame(cell.name=names(seu.itg@active.ident),cell.type=seu.itg@active.ident)

write.table(anno.tbl,"ctrl_R_integrated_anno.txt",quote=F,row.names=F,col.names=T)
saveRDS(seu.itg,"ctrl_integrated_Seu_annoed.rds")

###############################
#### BELOW DIDN'T WORK TOO FINE
###############################
markers <- FindAllMarkers(seu.itg, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.6)# 0.6==>~x1.5
markers.sign <- subset(markers, p_val_adj < 0.05)
## Unique markers for each cluster
library(dplyr)
for(cl in unique(markers.sign$cluster))
    cat("Cluster",cl,":",dim(subset(markers.sign, cluster == cl))[1],"markers\n")
markers.sign.2 <- as.data.frame(markers.sign %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC))
write.table(markers.sign,file="ctrl_markers.csv",sep=",",row.names=F,col.names=T,quote=F)

library(scCATCH)
x = scCATCH(markers.sign,"Human",tissue=c("Blood","Plasma","Peripheral blood","Serum","Colon","Colorectum","Large intestine","Gut"))
print(subset(x,select=c("cluster","cell_type")))

## They provide this collection of markers (but markers of which clusters??)
FeaturePlot(seu.itg,features=c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
                               "CD4", "NCAM1", "TBX21", "IFNG",
                               "GATA3", "CXCR5", "BCL6", "PDCD1",
                               "FOXP3", "FCER2", "ZBTB7B", "RUNX3",
                               "CD14", "ITGAM", "CEBPA", "CEBPB",
                               "SPI1", "FLT3", "ITGAX", "ZBTB46",
                               "CEACAM8", "ZBTB10", "IL3RA", "CLEC4C",
                               "NRP1", "NKG7", "KLRC1", "KLRC2",
                               "KLRC3", "KLRB1", "Ly6G6D", "ADGRE1",
                               "CD19", "SDC1", "MS4A1", "CD38",
                               "CD27", "TNFRSF13C", "CD84", "CD86",
                               "CD24", "TRAV", "TRAC", "TRBV", "TRBC",
                               "TRGV", "TRGC", "TRDV", "TRDC", "CR2",
                               "CD22", "IGHV", "IGHG", "IGHA", "IGHD",
                               "IGHM", "IGHE", "IGLV", "IGLL", "IGKV",
                               "IGLC", "IGKC"))
### Cell-type annotations were made using SCSA:
## https://github.com/bioinfo-ibms-pumc/SCSA (Cloned Aug 25, 2020)
## python3 SCSA.py -d whole.db -i ~/Documents/luni/data_analysis/hyper_inflammation/data_sets/UC/ctrl_markers.csv -s seurat -E -g Human -f 1.0 -o ~/Documents/luni/data_analysis/hyper_inflammation/data_sets/UC/ctrl_cell_type_anno.txt -m txt

## -f parameter was varied a bit to see annotations for smaller clusters (smaller FC)
## -f: 0.5, 0.7, 0.9, 1.0, 1.3, 1.5
## In the paper they first classify into four big clusters: T, B, NK, Myeloid

## clu.anno.f05 = read.table("ctrl_cell_type_anno_f0.5.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f07 = read.table("ctrl_cell_type_anno_f0.7.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f09 = read.table("ctrl_cell_type_anno_f0.9.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f1 = read.table("ctrl_cell_type_anno_f1.0.txt",header=T,sep="\t",
##                          check.names=F,stringsAsFactors=F)
## clu.anno.f13 = read.table("ctrl_cell_type_anno_f1.3.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f15 = read.table("ctrl_cell_type_anno_f1.5.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno = list(clu.anno.f05, clu.anno.f07, clu.anno.f09, clu.anno.f1,
##                 clu.anno.f13, clu.anno.f15)
## names(clu.anno) = c("f0.5","f0.7","f0.9","f1.0","f1.3","f1.5")
## saveRDS(clu.anno,"ctrl_cell_type_anno_all_f.rds")

## final.anno <- c()
## for(cl in 0:28){
##     cat("Cluster:",cl,"\n")
##     x <- unique(unlist(lapply(clu.anno,function(x) {
##         y <- subset(x, Cluster == cl)
##         y <- y[y[,2] > 0,]
##         ##y$Z.ratio <- y[,2]/min(y[,2])
##         ##y[which.max(y[,2]),1]
##     })))
##     final.anno <- c(final.anno, paste(x, collapse=","))
## }
## final.anno.tbl <- data.frame(Cluster=0:28,Cell.type=final.anno)

## new.cluster.names = vector("character",29)
## names(new.cluster.names) = seq(0,28,1)
## ## unique cell type will be taken as such
## ## non-unique are further checked for the second significant entries/cell-types
## new.cluster.names[c("2","12","15")] <- "T.cell"
## new.cluster.names[c("5","10","13","14","18","20","21","22")] <- "B.cell"
## new.cluster.names["23"] <- "Mac"

## ## > new.cluster.names
## ##          0          1          2          3          4          5          6 
## ##   "T.cell"   "T.cell"   "T.cell"   "T.cell"  "NK.cell"   "B.cell"   "B.cell" 
## ##          7          8          9         10         11         12         13 
## ##   "B.cell"   "T.cell"     "Treg"   "B.cell"   "B.cell"   "T.cell"   "B.cell" 
## ##         14         15         16         17         18         19         20 
## ##   "B.cell"   "T.cell"   "B.cell"   "B.cell"   "B.cell"  "NK.cell"   "B.cell" 
## ##         21         22         23         24         25         26         27 
## ##   "B.cell"   "B.cell"      "Mac"   "B.cell"   "B.cell" "NK2.cell"  "Unknown" 
## ##         28 
## ##  "Unknown" 

## seu.itg = RenameIdents(seu.itg, new.cluster.names)
## saveRDS(seu.itg,"ctrl_integrated_Seu_annoed.rds")

## Save into InterCom-pipeline format
## write.table(seu.itg[["RNA"]]@data,"ctrl_R_integrated_dat.txt",quote=F,row.names=T,col.names=T)
## anno.tbl <- data.frame(cell.name=rownames(seu.itg@meta.data),cell.type=seu.itg@meta.data$Cell.type)

## write.table(anno.tbl,"ctrl_R_integrated_anno.txt",quote=F,row.names=F,col.names=T)



### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Cell-type annotations for UC (manual, canonical clusters)
### =================================================================
seu.itg = readRDS("uc_integrated_Seu_clustered.rds")
DefaultAssay(seu.itg) = "RNA"
## Annotate larger groups of cells
new.cluster.names <- vector("character",length(unique(Idents(seu.itg))))
names(new.cluster.names) <- 0:(length(unique(Idents(seu.itg)))-1)

## =======
## T-cells
FeaturePlot(seu.itg,features=c("CD3D","CD3E","CD3G"))
##=> Clusters 1,3,5,6,2 are T-cells

## CD8+ T-cells (cytotoxic) and CD4+ ones
FeaturePlot(seu.itg,features=c("CD8A","CD8B","CD4"))
##=> Cluster 2 CD8+ T-cells
new.cluster.names["2"] <- "CD8.T.cell"

## Tregs: FOXP3+
FeaturePlot(seu.itg,features=c("CD8A","CD8B","CD4","FOXP3"))
new.cluster.names["5"] <- "Treg.cell"

##=> Clusters 1,3,6 CD4+ T cells
new.cluster.names[c("1","3","6")] <- "CD4.T.cell"

## ===================
## NK: CD3-NCAM+FCGR3+
FeaturePlot(seu.itg,features=c("CD3D","CD3E","CD3G","NCAM1","FCGR3A"))
## Cluster 10
new.cluster.names["10"] <- "NK.cell"

## ===================
## B-cells
## Plasma cells (final diff stage of B): CD79A+SDC1+CD78+(gene name??)IL6R+, CD19-CD20-/MS4A1-
FeaturePlot(seu.itg,features=c("CD79A","CD19","MS4A1","SDC1","IL6R"))
##=> Clusters 4,9,8 Plasma
new.cluster.names[c("4","9","8")] <- "Plasma.cell"
##=> Clusters 0,7,13,14,12 B cells
new.cluster.names[c("0","7","13","14","12")] <- "B.cell"

## =====================
## Macrophages/Monocytes
FeaturePlot(seu.itg,features=c("CD14","FCGR3A"))
##=> Cluster 11
new.cluster.names["11"] <- "Mac"

## FindMarkers(seu.itg,ident.1=15,ident.2=c(0,7,3,1),only.pos=T, min.pct=0.25,logfc.threshold=0.6)
## Identifies Cluster 15 as DC.cell from PanglaoDB by markers:
## TYROBP,FCER1G,GSN,TRDC,SPINK2,CST3
new.cluster.names["15"] <- "DC.cell"

## FindMarkers(seu.itg,ident.1=16,ident.2=c(0,7,3,1),only.pos=T, min.pct=0.25,logfc.threshold=0.6)
## Identifies Cluster 16 as Basal.cell from PanglaoDB by markers.
## Is it contamination from non-immune cells?
## Markers: CST3,IGFBP7,CAV1,GSN,IGFBP4,PRSS23
new.cluster.names["16"] <- "Unknown"


## =================================================
seu.itg = RenameIdents(seu.itg, new.cluster.names)
seu.itg$Cell.types = Idents(seu.itg)

write.table(seu.itg[["RNA"]]@data,"uc_R_integrated_dat.txt",quote=F,row.names=T,col.names=T)
anno.tbl <- data.frame(cell.name=names(seu.itg@active.ident),cell.type=seu.itg@active.ident)

write.table(anno.tbl,"uc_R_integrated_anno.txt",quote=F,row.names=F,col.names=T)
saveRDS(seu.itg,"uc_integrated_Seu_annoed.rds")

### ++++++++++++++++++++++++++++++++++++
### BELOW DIDN'T WORK
### ++++++++++++++++++++++++++++++++++++
### Cell-type annotations were made using SCSA:
## https://github.com/bioinfo-ibms-pumc/SCSA (Cloned Aug 25, 2020)
## python3 SCSA.py -d whole.db -i ~/Documents/luni/data_analysis/hyper_inflammation/data_sets/UC/uc_markers.csv -s seurat -E -g Human -f 0.5 -o ~/Documents/luni/data_analysis/hyper_inflammation/data_sets/UC/uc_cell_type_anno_f0.5.txt -m txt

## -f parameter was varied a bit to see annotations for smaller clusters (smaller FC)
## -f: 0.5, 0.7, 0.9, 1.0, 1.3, 1.5
## In the paper they first classify into four big clusters: T, B, NK, Myeloid

## clu.anno.f05 = read.table("uc_cell_type_anno_f0.5.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f07 = read.table("uc_cell_type_anno_f0.7.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f09 = read.table("uc_cell_type_anno_f0.9.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f1 = read.table("uc_cell_type_anno_f1.0.txt",header=T,sep="\t",
##                          check.names=F,stringsAsFactors=F)
## clu.anno.f13 = read.table("uc_cell_type_anno_f1.3.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno.f15 = read.table("uc_cell_type_anno_f1.5.txt",header=T,sep="\t",
##                           check.names=F,stringsAsFactors=F)
## clu.anno = list(clu.anno.f05, clu.anno.f07, clu.anno.f09, clu.anno.f1,
##                 clu.anno.f13, clu.anno.f15)
## names(clu.anno) = c("f0.5","f0.7","f0.9","f1.0","f1.3","f1.5")
## saveRDS(clu.anno,"uc_cell_type_anno_all_f.rds")

## final.anno <- c()
## for(cl in 0:30){
##     cat("Cluster:",cl,"\n")
##     x <- unique(unlist(lapply(clu.anno,function(x) {
##         y <- subset(x, Cluster == cl)
##         y <- y[y[,2] > 0,]
##         ##y$Z.ratio <- y[,2]/min(y[,2])
##         y[which.max(y[,2]),1]
##     })))
##     final.anno <- c(final.anno, paste(x, collapse=","))
## }
## final.anno.tbl <- data.frame(Cluster=0:30,Cell.type=final.anno)

## new.cluster.names = vector("character",31)
## names(new.cluster.names) = seq(0,30,1)
## ## unique cell type will be taken as such
## ## non-unique are further checked for the second significant entries/cell-types
## new.cluster.names["4"] <- "Treg"
## new.cluster.names["6"] <- "T.cell"
## new.cluster.names[c("7","21")] <- "NK.cell"
## new.cluster.names[c("8","11","13","16","19","27","28")] <- "B.cell"
## new.cluster.names["25"] <- "Monocyte"
## new.cluster.names["29"] <- "NKT.cell"

## cl <- 7
## lapply(clu.anno,function(x) {
##     y <- subset(x, Cluster == cl)
##     y <- y[y[,2] > 0,]
## })

## ## > new.cluster.names
## ##             0             1             2             3             4 
## ##      "T.cell"      "B.cell"      "T.cell"     "NK.cell"        "Treg" 
## ##             5             6             7             8             9 
## ##      "T.cell"      "T.cell"     "NK.cell"      "B.cell"        "Treg" 
## ##            10            11            12            13            14 
## ##      "B.cell"      "B.cell"     "NK.cell"      "B.cell"     "NK.cell" 
## ##            15            16            17            18            19 
## ##      "B.cell"      "B.cell" "Plasma.cell"     "NK.cell"      "B.cell" 
## ##            20            21            22            23            24 
## ##   "Plasma.DC"     "NK.cell"      "B.cell" "Plasma.cell"      "B.cell" 
## ##            25            26            27            28            29 
## ##    "Monocyte"    "NKT.cell"      "B.cell"      "B.cell"    "NKT.cell" 
## ##            30 
## ##     "Unknown" 


## seu.itg = RenameIdents(seu.itg, new.cluster.names)
## seu.itg$Cell.type = Idents(seu.itg)
## saveRDS(seu.itg,"uc_integrated_Seu_annoed.rds")

## Save into InterCom-pipeline format
## write.table(seu.itg[["RNA"]]@data,"uc_R_integrated_dat.txt",quote=F,row.names=T,col.names=T)
## anno.tbl <- data.frame(cell.name=rownames(seu.itg@meta.data),cell.type=seu.itg@meta.data$Cell.type)

## write.table(anno.tbl,"uc_R_integrated_anno.txt",quote=F,row.names=F,col.names=T)
