d <- read.table("GSE141105_Joint_expression_matrix.txt",header=T,check.names=F,stringsAsFactors=F)

cell.annos <- colnames(d)
cell.anno.split <- strsplit(cell.annos,"_")

cell.names <- sapply(cell.anno.split,`[`,1)
tissue <- sapply(cell.anno.split,`[`,2)
condition <- sapply(cell.anno.split,`[`,3)
batch <- apply(sapply(cell.anno.split,`[`,4:5),2,paste,collapse="_")
cell.id <- sapply(cell.anno.split,`[`,6)

### Batch specific
b <- "mouse_6"
cnd <- "Sick"

idx <- batch == b & condition == cnd
print(sum(idx))

d1 <- d[,idx]
colnames(d1) <- cell.id[idx]
print(dim(d1))
print(d1[1:2,1:4])

anno.tbl <- data.frame(cell.name=cell.id[idx],cell.type=cell.names[idx])
print(dim(anno.tbl))
print(head(anno.tbl))
print(unique(anno.tbl$cell.type))

write.table(d1,paste0(b,"_",cnd,"_dat.txt"),quote=F,sep='\t',row.names=T,col.names=T)
write.table(anno.tbl,paste0(b,"_",cnd,"_anno.txt"),quote=F,sep='\t',row.names=F,col.names=T)

### Joined, batch-less
cnd <- "Sick"

idx <- condition == cnd
print(sum(idx))

d1 <- d[,idx]
colnames(d1) <- cell.id[idx]
print(dim(d1))
print(d1[1:2,1:4])

anno.tbl <- data.frame(cell.name=cell.id[idx],cell.type=cell.names[idx])
print(dim(anno.tbl))
print(head(anno.tbl))
print(unique(anno.tbl$cell.type))

write.table(d1,paste0(cnd,"_joined_dat.txt"),quote=F,sep='\t',row.names=T,col.names=T)
write.table(anno.tbl,paste0(cnd,"_joined_anno.txt"),quote=F,sep='\t',row.names=F,col.names=T)

