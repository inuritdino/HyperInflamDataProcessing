## Paper defined cluster names
cluster.map <- data.frame(cluster=1:11,name=c("Schwann1","Schwann2","Melanocyte",
                                              "Endothelial","Pericyte","Fibroblast",
                                              "RPE","B.cell","TorNK.cell","Macrophage",
                                              "Mast.cell"))

## Repeat below for all donors and conditions,
## depending on the donor build differently anno table
name <- "Macula_donor1"

fn <- "GSM4037981_macula_donor_1_expression.tsv"

d <- read.table(fn,header=T)

dat <- d[,4:dim(d)[2]]
d[1:5,1:9]
dat <- t(dat)
colnames(dat) <- d$barcode
dat[1:5,1:9]

## For donors 1-3 Macula+peripheral
anno <- data.frame(cell.name=d$barcode,
                   cell.type=
                       sapply(d$final_cluster_labels,function(x)
                           cluster.map$name[cluster.map$cluster == x]))
## For donors 4-7 Macula+peripheral enriched
## We do not know cell types for all clusters, just retain cluster number for now
anno <- data.frame(cell.name=d$barcode,
                   cell.type=d$final_cluster_labels)

write.table(anno,paste0(name,"_anno.txt"),quote=F,sep='\t',row.names=F,col.names=T)
write.table(dat,paste0(name,"_dat.txt"),quote=F,sep='\t',row.names=T,col.names=T)
