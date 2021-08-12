arg <- commandArgs(T)

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];


data=read.table(dataFile,sep='\t',header=T,row.names=1)
group=read.table(groupFile,sep='\t',header=T)
rownames(group)=as.character(group[,1])
data=data[,rownames(group)]
re=matrix(NA,nrow=nrow(data),ncol=2)
rownames(re)=rownames(data)
colnames(re)=c('p','p.adj')
gg=unique(group[,2])
	
for(j in 1:nrow(data)){
	j1=group[,2]%in% gg[1]
	j2=group[,2]%in% gg[2]
	d1=as.numeric(data[j,j1])
	d2=as.numeric(data[j,j2])		
	Pvalue <- wilcox.test(d1,d2)$p.value;
	re[j,1]=Pvalue
}
re[,2]=p.adjust(re[,1],method='fdr')
write.table(re,paste(outDir,"/",SamID,".p.txt",sep=""),sep='\t',quote=F,row.names=T,col.names=NA)



