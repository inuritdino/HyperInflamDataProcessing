library(Matrix)
meta <- read.table("20180822_PolypAll_cleaned_metadata.txt",header=T)
unique(meta$diagnosis)

d <- read.table("20180822_PolypAll_cleaned_rawdata.txt",header=T)
rownames(d) <- d$X
d <- d[,-1]

dim(subset(meta, diagnosis == "AERD"))
dim(subset(meta, diagnosis == "HC"))

#####
## HC
idx <- meta$diagnosis == "HC"
anno.hc <- data.frame(cell.name = meta$X[idx], cell.type = meta$subset[idx])
dat.hc <- d[,anno.hc$cell.name]
## Make the InterCom input-format
write.table(anno.hc,"HC_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)
write.table(dat.hc,"HC_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

#######
## AERD
idx <- meta$diagnosis == "AERD"
anno.aerd <- data.frame(cell.name = meta$X[idx], cell.type=meta$subset[idx])
dat.aerd <- d[,anno.aerd$cell.name]
## Make the InterCom input-format
write.table(anno.aerd,"AERD_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)
write.table(dat.aerd,"AERD_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

######
## CRS
idx <- meta$diagnosis == "CRS"
anno.crs <- data.frame(cell.name = meta$X[idx], cell.type=meta$subset[idx])
dat.crs <- d[,anno.crs$cell.name]
## Make the InterCom input-format
write.table(anno.crs,"CRS_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)
write.table(dat.crs,"CRS_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

#########
## CRSwNP
idx <- meta$diagnosis == "CRSwNP"
anno.crswnp <- data.frame(cell.name = meta$X[idx], cell.type=meta$subset[idx])
dat.crswnp <- d[,anno.crswnp$cell.name]
## Make the InterCom input-format
write.table(anno.crswnp,"CRSwNP_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)
write.table(dat.crswnp,"CRSwNP_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)
