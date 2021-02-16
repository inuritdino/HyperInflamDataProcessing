smpls <- list.files(".","*_genes.tsv.gz")
smpls <- unlist(strsplit(smpls,"_CSF_GRCh38_genes.tsv.gz"))
smpls <- smpls[!grepl("_PBMCs_",smpls)]
smpls

for(i in 1:length(smpls)){
    if(dir.exists(smpls[i]))
        unlink(smpls[i],recursive=T)
    dir.create(smpls[i])
    print(grep("matrix",grep(smpls[i],list.files(),value=T),value=T))
    file.copy(grep("genes",grep(smpls[i],list.files(),value=T),value=T),
              paste0(smpls[i],"/features.tsv.gz"))
    file.copy(grep("barcodes",grep(smpls[i],list.files(),value=T),value=T),
              paste0(smpls[i],"/barcodes.tsv.gz"))
    file.copy(grep("matrix",grep(smpls[i],list.files(),value=T),value=T),
              paste0(smpls[i],"/matrix.mtx.gz"))
}


### Seurat analysis
library(Seurat)
combined.ctrl <- list()
combined.ms <- list()
for(i in 1:length(smpls)){
    cnts <- Read10X(smpls[i])
    seu <- CreateSeuratObject(cnts)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=2000)
    if(grepl("_MS",smpls[i]))
        combined.ms <- c(combined.ms,seu)
    else
        combined.ctrl <- c(combined.ctrl,seu)
}
saveRDS(combined.ms,"ms_Seurat_list.rds")
saveRDS(combined.ctrl,"ctrl_Seurat_list.rds")

## =========================================
## Anchors
anchors.ctrl <- FindIntegrationAnchors(combined.ctrl,dims=1:20,anchor.features=2000)
saveRDS(anchors.ctrl,"ctrl_integrate_anchors.rds")
anchors.ms <- FindIntegrationAnchors(combined.ms,dims=1:20,anchor.features=2000)
saveRDS(anchors.ms,"ms_integrate_anchors.rds")

## ============================================
## Common features
all.feats <- list(rownames(combined.ctrl[[1]]),rownames(combined.ctrl[[2]]),
                  rownames(combined.ctrl[[3]]),rownames(combined.ctrl[[4]]),
                  rownames(combined.ctrl[[5]]),rownames(combined.ctrl[[6]]))
## Another approach to use Ligs,Recs,TFs,Signalling molecules which we will use in our analyses
load("~/Documents/luni/data_analysis/HeartGladstone/intercom-v2.0/HUMAN_Background_data_v2.RData",verb=T)
genes <- c(unique(tf.db$Symbol),unique(LR$Ligand),unique(LR$Receptor),
           unique(Background_signaling_interactome$Source),
           unique(Background_signaling_interactome$Target))
genes <- unique(genes)
common.feats <- Reduce(intersect, all.feats, genes)
saveRDS(common.feats,"ctrl_common_feats.rds")
## ---
all.feats <- list(rownames(combined.ms[[1]]),rownames(combined.ms[[2]]),
                  rownames(combined.ms[[3]]),rownames(combined.ms[[4]]),
                  rownames(combined.ms[[5]]),rownames(combined.ms[[6]]))
## Another approach to use Ligs,Recs,TFs,Signalling molecules which we will use in our analyses
load("~/Documents/luni/data_analysis/HeartGladstone/intercom-v2.0/HUMAN_Background_data_v2.RData",verb=T)
genes <- c(unique(tf.db$Symbol),unique(LR$Ligand),unique(LR$Receptor),
           unique(Background_signaling_interactome$Source),
           unique(Background_signaling_interactome$Target))
genes <- unique(genes)
common.feats <- Reduce(intersect, all.feats, genes)
saveRDS(common.feats,"ms_common_feats.rds")


### ????
## Save plots remotely to plot them locally
plot.iris <- function(gg.obj, file){
    pdf(file)
    print(gg.obj)
    dev.off()
}

### ======================================================
### Controls
anchors <- readRDS("ctrl_integrate_anchors.rds")
common.feats <- readRDS("ctrl_common_feats.rds")
ctrl.integrated <- IntegrateData(anchorset = anchors, dims = 1:20, features = common.feats)
saveRDS(ctrl.integrated,"ctrl_Seurat_integrated.rds")

ctrl.integrated <- ScaleData(ctrl.integrated, verbose = FALSE)
ctrl.integrated <- RunPCA(ctrl.integrated, npcs = 30, verbose = FALSE)
ctrl.integrated <- RunUMAP(ctrl.integrated, dims = 1:30)
g <- DimPlot(ctrl.integrated)
ctrl.integrated <- FindNeighbors(ctrl.integrated, dims = 1:20)
ctrl.integrated <- FindClusters(ctrl.integrated, resolution = seq(0.4,1.6,by=0.3))
head(Idents(ctrl.integrated))
saveRDS(ctrl.integrated,"ctrl_Seurat_integrated1.rds")

### ----
### MS
anchors <- readRDS("ms_integrate_anchors.rds")
common.feats <- readRDS("ms_common_feats.rds")
ms.integrated <- IntegrateData(anchorset = anchors, dims = 1:20, features = common.feats)
saveRDS(ms.integrated,"ms_Seurat_integrated.rds")

ms.integrated <- ScaleData(ms.integrated, verbose = FALSE)
ms.integrated <- RunPCA(ms.integrated, npcs = 30, verbose = FALSE)
ms.integrated <- RunUMAP(ms.integrated, dims = 1:30)
g <- DimPlot(ms.integrated)
ms.integrated <- FindNeighbors(ms.integrated, dims = 1:20)
ms.integrated <- FindClusters(ms.integrated, resolution = seq(0.4,1.6,by=0.3))
head(Idents(ms.integrated))
saveRDS(ms.integrated,"ms_Seurat_integrated1.rds")

### ======================================================
## Biomarkers

paper.ids <- list(
    cd4.T.pct=c("CD3E","LCK","TRAC","TRAJ16","IL7R","CD4"),
    cd8a.T.pct=c("CD3E","LCK","TRAC","TRAJ16","CD8B","CCL5"),
    cd8na.T.pct=c("CD3E","LCK","TRAC","TRAJ16","CD8B","CCR7"),
    reg.T.pct=c("CD3E","LCK","TRAC","TRAJ16","FOXP3","CTLA4"),
    gd.T.pct=c("CD3E","LCK","TRAC","TRAJ16","TRDC"),
    nk1.pct=c("GNLY","NKG7","FCGR3A","PRF1"),
    nk2.pct=c("GNLY","NKG7","SELL","XCL1"),
    b1.pct=c("CD74","CD79A","CD37","IGHD"), ## all B MS4A1-,SDC1-
    b2.pct=c("CD74","CD79A","CD27","IGHM"),
    plasma.pct=c("CD74","CD79A","CD38","IGHG","TNFRSF17"),
    mDC1.pct=c("LYZ","WDFY4","XCR1","BATF3"),
    mDC2.pct=c("LYZ","FCER1A","CD1C","CLEC10A"),
    gran.pct=c("LYZ","S100A8","S100A9"),
    mono1.pct=c("FCGR3A"),
    mono2.pct=c("CD14"),
    pDC.pct=c("TCF4","TNFRSF21"),
    megaK.pct=c("GNG11","CLU")
)

get.identity <- function(seu,ids,cluster.col="seurat_clusters"){
    ns <- c()
    for(i in 1:length(ids)){
        idx <- ids[[i]] %in% rownames(seu)
        if(sum(idx) == 0)
            next
        if( sum(idx) != length(ids[[i]]) ){
            cat(names(ids)[i],":",ids[[i]][!idx],"not found\n")
            cat("taking only:", ids[[i]][idx],"\n")
            ##seu[[names(ids)[i]]] <- PercentageFeatureSet(seu, features=ids[[i]][idx])
        }
        seu[[names(ids)[i]]] <- apply(seu[['integrated']]@scale.data[ids[[i]][idx],,drop=FALSE],2,sum)
        ns <- c(ns, names(ids)[i])
    }
    marker.pct <- seu@meta.data[,ns]
    voting <- data.frame(cell.type=colnames(marker.pct)[apply(marker.pct,1,which.max)],
                         cluster=seu@meta.data[,cluster.col])
    tbl <- table(voting)
    heatmap(tbl,scale="column",Rowv=NA,Colv=NA,keep.dendro=F)
    print(tbl)
    new.ids <- rownames(tbl)[apply(tbl,2,which.max)]
    names(new.ids) <- colnames(tbl)
    new.ids <- gsub("\\.pct","",new.ids) ## crucial
    print(new.ids)
    seu <- RenameIdents(seu, new.ids)
    print(head(Idents(seu)))
    return(seu)
}
## Load Seurat object here
seu1 <- get.identity(seu, paper.ids)



