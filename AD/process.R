smpls <- list.files(".","_matrix.mtx.gz")
smpls <- sapply(strsplit(smpls,"_matrix.mtx.gz"),`[`,1)
## They took initially only blister AD/HC samples and AD-biopsies
## AD-biopsy was used to compare against (for integration???)
## See mmc14.xlsx for annotations
## smpls <- grep("AD[1-8]+|HC[1-5]+",smpls,value=T)

files <- list.files(".","*.gz")
require(Seurat)
seu.lst <- vector("list",length(smpls))
names(seu.lst) <- smpls
for(i in 1:length(smpls)){
    print(smpls[i])
    if(!dir.exists("tmp")){
        dir.create("tmp")
    }
    else{
        unlink("tmp",recursive=TRUE)
        dir.create("tmp")
    }
    file.copy(paste0(smpls[i],"_barcodes.tsv.gz"),"tmp/barcodes.tsv.gz")
    file.copy(paste0(smpls[i],"_features.tsv.gz"),"tmp/features.tsv.gz")
    file.copy(paste0(smpls[i],"_matrix.mtx.gz"),"tmp/matrix.mtx.gz")
    cnts <- Read10X("tmp")
    seu.lst[[i]] <- CreateSeuratObject(cnts,min.cells=3,min.features=200)
    seu.lst[[i]][["percent.mt"]] <- PercentageFeatureSet(seu.lst[[i]], pattern = "^MT-")
    seu.lst[[i]] <- subset(seu.lst[[i]], subset = nFeature_RNA > 500 & percent.mt < 25)# 25?
    seu.lst[[i]] <- NormalizeData(seu.lst[[i]])
    seu.lst[[i]] <- FindVariableFeatures(seu.lst[[i]])
}

saveRDS(seu.lst,"all_smpls_seu_list.rds")

## =======================================
## Integration
## =======================================
require(Seurat)
seu.lst <- readRDS("all_smpls_seu_list.rds")

## HC1 sample cause problem in anchor-finding
seu.lst <- seu.lst[-9] ## remove HC1 with 6694 genes and 79 samples only
anchors <- FindIntegrationAnchors(object.list = seu.lst)
##saveRDS(anchors,"anchors.rds")

## Integrate
anchors <- readRDS("anchors.rds")
seu.itg <- IntegrateData(anchors)
saveRDS(seu.itg,"integrated_Seu.rds")

## Cluster
seu.itg <- ScaleData(seu.itg, verbose = FALSE)
seu.itg <- RunPCA(seu.itg, verbose = FALSE, npcs=150)

seu.itg <- JackStraw(seu.itg, num.replicate=100, dims=150)
seu.itg <- ScoreJackStraw(seu.itg, dims = 1:150)
saveRDS(seu.itg,"integrated_Seu_jackstrawed.rds")

JackStrawPlot(seu.itg, dims = 41:50)

seu.itg <- RunUMAP(seu.itg, dims = 1:30) ## 30??

seu.itg <- FindNeighbors(seu.itg)
seu.itg <- FindClusters(seu.itg, resolution=0.4) ## the only param they mention!

saveRDS(seu.itg,"integrated_Seu_clustered.rds")
