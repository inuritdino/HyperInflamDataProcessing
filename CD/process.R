library(Seurat)
library(sctransform)
library(dplyr)
library(SingleR)
library(ggplot2)

##data: GEO:GSE134809 

##data files(matrix.mtx,genes.tsv,barcodes.tsv) associated-
#-with a given sample number organised in a folder named by the sample number: as required for Read10x()

#'path'<-"path to the folder containing all data"
dirs<-list.dirs(path,recursive = F)
samps<-list.dirs(path,recursive = F,full.names = F)
##"pat_samp_data.txt" - metadata of samples as given on sheet1, supplementary file2(excel workbook)(PMID: 31474370)
pat_samp_data<-read.table("pat_samp_data.txt",header = T,sep="\t")
samps<-samps[which(samps%in%pat_samp_data$Sample_ID)]
read<-dirs[which(samps%in%pat_samp_data$Sample_ID)]
pat_samp_data=pat_samp_data[which(pat_samp_data$Sample_ID%in%samps),]

##Reading raw UMI counts matrices; filtering cells with more than 800 UMIs(as done in the paper)
data_all<-lapply(1:length(read),function(x){
  data<-Read10X(data.dir=read[x])
  id=which(pat_samp_data$Sample_ID%in%samps[x])
  project=paste(samps[x],"_",pat_samp_data$Patient.ID[id],"_",pat_samp_data$status[id])
  data_seu<-CreateSeuratObject(counts=data,min.features = 800,project=project)
})
names(data_all)=samps

##filtering out cells with at least 25% mtmRNA(as done in the paper)
for(i in 1:length(data_all)){
  data_all[[i]][["percent.mt"]] <- PercentageFeatureSet(data_all[[i]], pattern = "^MT-")
  data_all[[i]]<-subset(data_all[[i]],percent.mt<25)
}

##filtering out cells with more than 1% UMIs associated with epithelial cells (gene list given in the paper):
##PLA2G2A, CLCA1, REG4, S100A14, ITLN1, ELF3, PIGR, EPCAM, REG1B, REG1A, REG3A, FABP1, RBP2,
##SST, FABP2, SPINK1, FABP6, AGR2, AGR3, CLDN3, CLDN4, DEFA6, DEFA5, SPINK4, ALDOB, LCN2, MUC2, KRT8, KRT18,
##TSPAN8, OLFM4, GPX2, IFI27, PHGR1, MT1G, CLDN7, KRT19, FXYD3, LGALS4, FCGBP, TFF3, TFF1

genes_epi=as.vector(read.table("epithelial_genes.txt",sep="\n",header = F)[,1])
genes_epi=genes_epi[-which(genes_epi=="MUC2")]
for(i in 1:length(data_all)){
  data_all[[i]][["percent.epiGenes"]] <- PercentageFeatureSet(data_all[[i]], features=genes_epi)
  data_all[[i]]<-subset(data_all[[i]],percent.epiGenes<1)
}
#saveRDS(data_all,file="data_all.rds")

##Preparing data for integration
for(i in 1:length(data_all)){
  data_all[[i]]<-SCTransform(data_all[[i]])
}
data_all.features <- SelectIntegrationFeatures(object.list=data_all, nfeatures = 3000)
data_all<-PrepSCTIntegration(object.list = data_all,anchor.features = data_all.features)
#saveRDS(data_all,"data_all_sct.rds")

##Find anchors and integrate
data_all.anchors<-FindIntegrationAnchors(object.list = data_all
                                         , normalization.method = "SCT", 
                                         anchor.features = data_all.features
                                         ,k.filter=50)
#saveRDS(data_all.anchors,"data_all_anc.rds")

##Data integration
data_all.anchors<-readRDS("data_all_anc.rds")
data_all.int<-IntegrateData(anchorset = data_all.anchors,normalization.method = "SCT")
#save(data_all.int,"data_all.integrated.rds")

##Clustering integrated data
#data_all.int<-readRDS("data_all.integrated.rds")
DefaultAssay(data_all.int)<-"integrated"
data_all.int<-ScaleData(data_all.int)
data_all.int<-RunPCA(data_all.int,npcs=45)
data_all.int<-RunUMAP(data_all.int,reduction = "pca",dims = 1:30)
data_all.int<-FindNeighbors(data_all.int,dims=1:30,reduction="pca")
data_all.int<-FindClusters(data_all,resolution = 1,assay="integrated")

##Identifying cell markers associated with clusters 
integrated.markers<-FindAllMarkers(data_all.int,assay="RNA",only.pos = T)
markers<- integrated.markers %>% group_by(cluster) %>% top_n(n=15,wt=avg_logFC)
write.csv(markers,"markers_top15.csv")
markers<-split(markers,markers$cluster)

##Cell type annotation with SingleR
ref_data2<-BlueprintEncodeData()
c_types_SingleR <- SingleR(test = data_all.int$RNA@data, 
                           ref = ref_data2, 
                           labels = ref_data$label.fine, 
                           method = "cluster", 
                           clusters = Idents(data_all.int))

c_types=c_types_SingleR$first.labels[as.numeric(Idents(data_all.int))]
data_all.int[["cell_types"]]=c_types

sample<-sapply(data_all.int$orig.ident,function(x){
  dat=strsplit(x," ")[[1]][1]
})

patient<-sapply(data_all.int$orig.ident,function(x){
  dat=strsplit(x," ")[[1]][4]
})
type<-sapply(data_all.int$orig.ident,function(x){
  dat=strsplit(x," ")[[1]][6]
})

data_all.int[["sample"]]<-sample
data_all.int[["patient"]]<-patient
data_all.int[["disease_state"]]<-type

saveRDS(data_all.int,"data_all_int_annotated.rds")

##Data visualtion: Umap and dotplots

DimPlot(data_all.int,reduction="umap",group.by = "cell_types",label=T,split.by = "disease_state")

cluster=unique(cbind(data_all.int$c_types,data_all.int$seurat_clusters))

dotplot=list()

for(i in 1:length(markers)){
  dotplot[[i]]=DotPlot(object=data_all.int,assay="RNA"
                       ,group.by="cell_types"
                       , split.by="disease_state"
                       ,features=markers[[i]]$gene,
                       cols=c("red","grey","blue"))+
    RotatedAxis()+
    labs(x=NULL,y=NULL,title=cluster[which(cluster[,2]==i),])
}


