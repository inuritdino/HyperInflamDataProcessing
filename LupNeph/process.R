info = read.table("41590_2019_386_MOESM3_ESM.csv",header=T,sep=',')
meta = read.table("metadata.csv",header=T,sep=',')
meta1 <- read.table("meta.txt",header=T,sep='\t')
rownames(meta1) <- meta1$NAME # for simpler slicing

dat <- read.table("exprMatrixSleMetro.tsv",header=T,sep='\t')
rownames(dat) <- dat[,1]

dat <- dat[,-1]

## Kidney
idx <-info$Disease == "Healthy"
info$Patient_ID[idx]
idx2 <- unlist(lapply(info$Patient_ID[idx],function(x) which(meta$SUBJECT_ID == x &
                                                             meta$TISSUE == "Kidney")))

meta.healthy  <- meta[idx2,]
healthy.dat <- dat[,as.character(meta.healthy$SAMPLE_ID[meta.healthy$SAMPLE_ID %in% colnames(dat)])]
write.table(healthy.dat,"kidney_healthy_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

anno.healthy <- data.frame(cell.name=colnames(healthy.dat), cell.type=meta1[colnames(healthy.dat),2])
write.table(anno.healthy,"kidney_healthy_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)

#---
idx <- info$Disease == "Lupus"
info$Patient_ID[idx]
idx2 <- unlist(lapply(info$Patient_ID[idx],function(x) which(meta$SUBJECT_ID == x &
                                                             meta$TISSUE == "Kidney")))

meta.lupus <- meta[idx2,]
lupus.dat <- dat[,as.character(meta.lupus$SAMPLE_ID[meta.lupus$SAMPLE_ID %in% colnames(dat)])]
write.table(lupus.dat,"kidney_lupus_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

anno.lupus <- data.frame(cell.name=colnames(lupus.dat), cell.type=meta1[colnames(lupus.dat),2])
write.table(anno.lupus,"kidney_lupus_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)

## Skin
idx <-info$Disease == "Healthy"
info$Patient_ID[idx]
idx2 <- unlist(lapply(info$Patient_ID[idx],function(x) which(meta$SUBJECT_ID == x &
                                                             meta$TISSUE == "Skin")))

meta.healthy  <- meta[idx2,]
healthy.dat <- dat[,as.character(meta.healthy$SAMPLE_ID[meta.healthy$SAMPLE_ID %in% colnames(dat)])]
write.table(healthy.dat,"skin_healthy_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

anno.healthy <- data.frame(cell.name=colnames(healthy.dat), cell.type=meta1[colnames(healthy.dat),2])
write.table(anno.healthy,"skin_healthy_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)

#---
idx <- info$Disease == "Lupus"
info$Patient_ID[idx]
idx2 <- unlist(lapply(info$Patient_ID[idx],function(x) which(meta$SUBJECT_ID == x &
                                                             meta$TISSUE == "Skin")))

meta.lupus <- meta[idx2,]
lupus.dat <- dat[,as.character(meta.lupus$SAMPLE_ID[meta.lupus$SAMPLE_ID %in% colnames(dat)])]
write.table(lupus.dat,"skin_lupus_dat.txt",quote=F,sep='\t',row.names=T,col.names=T)

anno.lupus <- data.frame(cell.name=colnames(lupus.dat), cell.type=meta1[colnames(lupus.dat),2])
write.table(anno.lupus,"skin_lupus_anno.txt",quote=F,sep='\t',row.names=F,col.names=T)
