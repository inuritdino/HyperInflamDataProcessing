## Macrophages
dmac <- read.table("merged_plaque_gex-umi-data-mnn-cat_macrophages.txt",header=T,check.names=F,
                   stringsAsFactors=F,sep='\t')
rownames(dmac) <- dmac[,1]
dmac <- dmac[,-1]
## T-cells
dtc <- read.table("merged_plaque_gex-umi-data-mnn-cat_t_cells.txt",header=T,check.names=F,
                   stringsAsFactors=F,sep='\t')
rownames(dtc) <- dtc[,1]
dtc <- dtc[,-1]

## Merge
## Is the order of genes in two tables the same??
## paste0(rownames(dtc),collapse="") == paste0(rownames(dmac),collapse="")
d0 <- cbind(dmac,dtc)

n <- colnames(d0)
ns <- strsplit(n,',\\s(\\w+|\\w+\\s\\w+): ')

cell.name <- sapply(ns,`[`,1)
batch <- sapply(ns,`[`,2) # is it a batch?
condition <- sapply(ns,`[`,3)
cell.type <- sapply(ns,`[`,4)
broad.ct <- sapply(ns,`[`,5)

####
cnd <- "Asymptomatic"
idx <- condition == cnd
print(sum(idx))

d1 <- d0[,idx]
colnames(d1) <- cell.name[idx]
print(dim(d1))
print(d1[1:2,1:5])

anno.tbl <- data.frame(cell.name=cell.name[idx],cell.type=cell.type[idx])
print(dim(anno.tbl))
print(head(anno.tbl))

write.table(d1,paste0(cnd,"_all_dat.txt"),quote=F,sep='\t',row.names=T,col.names=T)
write.table(anno.tbl,paste0(cnd,"_all_anno.txt"),quote=F,sep='\t',row.names=F,col.names=T)

