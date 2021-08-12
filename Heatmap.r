arg <- commandArgs(T)

library(pheatmap)


dataFile <- arg[1];
groupFile <- arg[2];
IDFile<-arg[3]
outDir <- arg[4];

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];

data=read.table(dataFile,sep='\t',header=T,row.names=1)
grp=read.table(groupFile,sep='\t',header=T)
rownames(grp)=as.character(grp[,1])
ID=read.table(IDFile,sep='\t',header=T)
rownames(ID)=as.character(ID[,1])
ii=intersect(rownames(grp),colnames(data))
dd=data[,ii]
gg=grp[ii,]
	
d=apply(dd,1,function(x){a=(x-min(x))/(max(x)-min(x))})
d=t(d)
g=as.data.frame(gg[colnames(d),2])
rownames(g)=colnames(d)
colnames(g)='group'
if(nrow(data)<150){
	pdf(paste(outDir,"/",SamID,"heatmap.pdf",sep=""),w=8,h=18)
} else {
	pdf(paste(outDir,"/",SamID,"heatmap.pdf",sep=""),w=8,h=27)
}
pheatmap(d,cluster_cols = FALSE,show_colnames=F,annotation_col=g,labels_row=ID[rownames(d),2])
dev.off()


	

