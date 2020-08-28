d <- read.table("counts_filt.csv",sep=",",header=T,check.names=F,stringsAsFactors=F)
ann <- read.table("publication_celldata.csv",sep=",",header=T,check.names=F,stringsAsFactors=F)
rownames(d) <- d[,1]
d <- d[,-1]

## Cell types and subtypes
cluster.names <- data.frame(cl=ann$cluster.name,hl=ann$hlad.cluster.name,fol=ann$follicle.cluster.name,stringsAsFactors=F)
new.cluster.names <- apply(cluster.names,1,function(x){
    make.names(paste(x[!is.na(x)],collapse="."))
})

ann$new.cluster.name <- new.cluster.names

## Divide by the tissue/condition
tt <- "psoriasis"
idx <- ann$tissue == tt
d1 <- d[,ann[idx,1]]
dim(d1)
d1[1:5,1:5]
ann.tbl <- data.frame(cell.name=colnames(d1),cell.type=ann$new.cluster.name[idx])
dim(ann.tbl)
head(ann.tbl)
unique(ann.tbl$cell.type)

write.table(d1,paste0(tt,"_dat.txt"),sep='\t',quote=F,row.names=T,col.names=T)
write.table(ann.tbl,paste0(tt,"_anno.txt"),sep="\t",quote=F,row.names=F,col.names=T)
