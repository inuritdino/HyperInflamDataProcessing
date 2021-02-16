library(Seurat)
library(ggplot2)

type <- 'umicount' ## 'readcount'

smpl <- "d3"
if(smpl == "c1")
    fn <- "GSM3823939_control.s1.dgecounts.rds"
else if(smpl == "c2")
    fn <- "GSM3823940_control.s2.dgecounts.rds"
else if(smpl == "c3")
    fn <- "GSM3823941_control.s3.dgecounts.rds"
else if(smpl == "d1")
    fn <- "GSM3823942_diabetes.s1.dgecounts.rds"
else if(smpl == "d2")
    fn <- "GSM3823943_diabetes.s2.dgecounts.rds"
else if(smpl == "d3")
    fn <- "GSM3823944_diabetes.s3.dgecounts.rds"

d <- readRDS(fn)
counts <- d[[type]]$inex$downsampling$downsampled_

write.table(rownames(counts),paste0(smpl,"_ensembl_genes.txt"),
            quote=F,sep='\t',row.names=F,col.names=F)
## In the meantime convert ENSMBL ID to Gene Symbols, manual work
gs <- read.table(paste0(smpl,"_genes.txt"),header=F)

## This shows many genes of RF-family,LINC(2)..we can safely remove duplicates
print(gs[which(duplicated(gs[,2])),2])

counts <- counts[-which(duplicated(gs[,2])),]
rownames(counts) <- gs[!duplicated(gs[,2]),2]


## Seurat
seu <- CreateSeuratObject(counts,min.cells=3) ##min.features=500?? params mentioned in Materials&Methods
saveRDS(seu,paste0(smpl,"_Seurat0_readsINEX.rds"))

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0(smpl,"_feature_plot.pdf"),device="pdf")

## For the parameters of the analysis see the supplement of the paper.
## Cleaning (the paper mentions the lower thr of 500 (but for what??
## Is it not confused with min.features=500 above???) and upper thr of
## Inf assuming this here... + MT 5%
seu <- subset(seu, subset = nFeature_RNA > 500 & percent.mt < 5)
## NOTE: for c2 sample above thrs are very strict... see the feature plot

seu <- NormalizeData(seu)

## The paper suggest using num.bin=30, x.low/high.cutoff == mean cutoff in Seurat 3
## y.cutoff == dispersion cutoff [1]
seu <- FindVariableFeatures(seu, num.bin=30, mean.cutoff=c(0.0125,8),
                            dispersion.cutoff=c(1,Inf))

seu <- ScaleData(seu)

## num of PC's 40 (from the paper, supplement)
seu <- RunPCA(seu, feature=VariableFeatures(object = seu), npcs=40)

seu <- FindNeighbors(seu, dims=1:20) # 20 dims from the paper
seu <- FindClusters(seu, resolution = 0.6) # 0.6 from the paper

seu <- RunUMAP(seu, dims = 1:12)
DimPlot(seu, reduction='umap')
ggsave(paste0(smpl,"_umap_plot.pdf"),device="pdf")

saveRDS(seu,paste0(smpl,"_Seurat1_umiINEX.rds"))

### ???????????????????????????????????????????????????????????
### Cluster separate datasets to get an idea
### ???????????????????????????????????????????????????????????
library(dplyr)
library(ggplot2)
library(Seurat)

super.markers <- c("CUBN","LRP2","SLC34A1","SLC5A12","SLC5A2","ALDOB","CFH","SLC12A1","SLC12A3","SLC12A2",
                   "SLC8A1","AQP2","AQP6","SLC26A4","ATP6V0D2","KIT","NPHS1","NPHS2","PECAM1","FLT1","ITGA8",
                   "PDGFRB","PTPRC")

kidney.ids <- list(
    pct.mrk.pct=c("CUBN","LRP2","SLC34A1","SLC5A12","SLC5A2","ALDOB"),
    cfh.mrk.pct=c("LRP2","CFH"),
    loh.mrk.pct=c("SLC12A1"),
    dct.mrk.pct=c("SLC12A3","SLC8A1"),
    dct.ct.mrk.pct=c("SLC12A2","SLC8A1"),
    cd.pc.mrk.pct=c("SLC8A1","AQP2"),
    cd.ica.mrk.pct=c("AQP6","ATP6V0D2","KIT"),
    cd.icb.mrk.pct=c("SLC26A4","ATP6V0D2"),
    podo.mrk.pct=c("NPHS1","NPHS2"),
    endo.mrk.pct=c("PECAM1","FLT1","ITGA8"),
    mes.mrk.pct=c("ITGA8","PDGFRB"),
    leuk.mrk.pct=c("PTPRC") #,"CD247","CSF2RA","MS4A1","PAX5","CD38","SDC1")
)

get.identity <- function(seu,ids){
    ns <- c()
    for(i in 1:length(ids)){
        ns <- c(ns, names(ids)[i])
        seu[[names(ids)[i]]] <- PercentageFeatureSet(seu, features=ids[[i]])
    }
    marker.pct <- seu@meta.data[,ns]
    voting <- data.frame(cell.type=colnames(marker.pct)[apply(marker.pct,1,which.max)],
                         cluster=seu@meta.data[,"seurat_clusters"])
    tbl <- table(voting)
    heatmap(tbl,scale="column",Rowv=NA,Colv=NA,keep.dendro=F)
    print(tbl)
    new.ids <- rownames(tbl)[apply(tbl,2,which.max)]
    names(new.ids) <- colnames(tbl)
    new.ids <- gsub("\\.mrk\\.pct","",new.ids) ## crucial
    print(new.ids)
    seu <- RenameIdents(seu, new.ids)
    print(head(Idents(seu)))
    return(seu)
}

smpl <- "c3"
seu <- readRDS(paste0(smpl,"_Seurat1_umiINEX.rds"))

seu1 <- get.identity(seu, kidney.ids)
DimPlot(seu1, label=T)

DotPlot(seu,features=super.markers)
ggsave(paste0(smpl,"_super_marker_dotplot.pdf"),device="pdf")

## Do we save the seu1 object with defined identities???
## ??



## Manual by observing the DotPlot
new.cluster.ids <- rep("",times=length(levels(Idents(seu))))
new.cluster.ids[14] <- "Endo"
new.cluster.ids[c(3,7)] <- "PCT"


markers <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

idx <- top10$gene %in% super.markers
top10$gene[idx]
top10$cluster[idx]

DoHeatmap(seu, features = top10$gene) + NoLegend()

### ???????????????????????????????????????????????????????????

### ==========================================================
### Integration of control samples into one
### ==========================================================
c1 <- readRDS("c1_Seurat1_umiINEX.rds")
c2 <- readRDS("c2_Seurat1_umiINEX.rds")
c3 <- readRDS("c3_Seurat1_umiINEX.rds")
ctrl.list <- list(ctrl1=c1,ctrl2=c2,ctrl3=c3)
saveRDS(ctrl.list,"c123_Seurat1_umiINEX_list.rds")

ctrl.list <- readRDS("c123_Seurat1_umiINEX_list.rds")
## Totally repeat the paper's analyses is not possible since they used an older version
## of Seurat with CCA analysis for integration. Here I will persue the novel approach,
## described in Seurat's vignettes. They mention 1000 features and 12 dims though...
anchors <- FindIntegrationAnchors(object.list = ctrl.list, dims = 1:12,anchor.features=1000)
saveRDS(anchors,"ctrl_integrate_anchors.rds")

## Features to use (all commont features in the three objects)
all.feats <- list(rownames(ctrl.list[[1]]),rownames(ctrl.list[[2]]),rownames(ctrl.list[[3]]))
##common.feats <- Reduce(intersect, all.feats)
## Another approach to use Ligs,Recs,TFs,Signalling molecules which we will use in our analyses
load("~/Documents/luni/data_analysis/HeartGladstone/intercom-v2.0/HUMAN_Background_data_v2.RData",verb=T)
genes <- c(unique(tf.db$Symbol),unique(LR$Ligand),unique(LR$Receptor),
           unique(Background_signaling_interactome$Source),
           unique(Background_signaling_interactome$Target))
genes <- unique(genes)
common.feats <- Reduce(intersect, all.feats, genes)
saveRDS(common.feats,"ctrl_common_feats.rds")

### Memory problem locally... Have to run it elsewhere...
anchors <- readRDS("ctrl_integrate_anchors.rds")
common.feats <- readRDS("ctrl_common_feats.rds")
ctrl.integrated <- IntegrateData(anchorset = anchors, dims = 1:12, features = common.feats)

### ==========================================================
### Integration of diabetes samples into one
### ==========================================================
d1 <- readRDS("d1_Seurat1_umiINEX.rds")
d2 <- readRDS("d2_Seurat1_umiINEX.rds")
d3 <- readRDS("d3_Seurat1_umiINEX.rds")
diab.list <- list(diab1=d1,diab2=d2,diab3=d3)
saveRDS(diab.list,"d123_Seurat1_umiINEX_list.rds")

diab.list <- readRDS("d123_Seurat1_umiINEX_list.rds")
## Totally repeat the paper's analyses is not possible since they used an older version
## of Seurat with CCA analysis for integration. Here I will persue the novel approach,
## described in Seurat's vignettes. They mention 1000 features and 12 dims though...
anchors <- FindIntegrationAnchors(object.list = diab.list, dims = 1:12, anchor.features=1000)
saveRDS(anchors,"diab_integrate_anchors.rds")

## Features to use (all commont features in the three objects)
all.feats <- list(rownames(diab.list[[1]]),rownames(diab.list[[2]]),rownames(diab.list[[3]]))
##common.feats <- Reduce(intersect, all.feats)
## Another approach to use Ligs,Recs,TFs,Signalling molecules which we will use in our analyses
load("~/Documents/luni/data_analysis/HeartGladstone/intercom-v2.0/HUMAN_Background_data_v2.RData",verb=T)
genes <- c(unique(tf.db$Symbol),unique(LR$Ligand),unique(LR$Receptor),
           unique(Background_signaling_interactome$Source),
           unique(Background_signaling_interactome$Target))
genes <- unique(genes)
common.feats <- Reduce(intersect, all.feats, genes)
saveRDS(common.feats,"diab_common_feats.rds")

### Memory problem locally... Have to run it elsewhere...
anchors <- readRDS("diab_integrate_anchors.rds")
common.feats <- readRDS("diab_common_feats.rds")
diab.integrated <- IntegrateData(anchorset = anchors, dims = 1:12, features = common.feats)
