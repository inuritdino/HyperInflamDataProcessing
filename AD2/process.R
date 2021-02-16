smpls <- list.files(pattern="*.clean.data.txt")

tot.n.cells <- 0
for(s in smpls){
    d <- read.table(s,header=T,sep=',',check.names=F,stringsAsFactors=F)
    rownames(d) <- d[,1]
    d <- d[,-1]
    tot.n.cells <- tot.n.cells + dim(d)[2]
}
