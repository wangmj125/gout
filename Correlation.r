arg <- commandArgs(T)
library('RColorBrewer')

dataFile <- arg[1];
groupFile <- arg[2];
No<-as.numeric(arg[3])
outDir <- arg[4];

data=read.table(dataFile,sep='\t',header=T,row.names=1)

d1=data[,1:No]
d2=data[,(No+1):ncol(data)]

col1=brewer.pal(8,"Set3")[c(1,3:6,8)]
col2=brewer.pal(7,"Set2")[c(1:7)]
col=c(col1,col2)
names(col)=colnames(d2)

group=read.table(groupFile,sep='\t',header=T)
rownames(group)=as.character(group[,1])

gg=unique(as.character(group[,2]))
pch=as.character(group[,2])
pch[which(pch==gg[1])]=3
pch[which(pch==gg[2])]=8
pch=as.numeric(pch)

for(i in 1:No){
	ff=colnames(d1)[i]
	pdf(paste(outDir,"/",ff,".pdf",sep=""),h=10,w=16)
	layout(matrix(1:15,nrow=3,byrow=T))
	xx=as.numeric(d1[,i])
	for(j in 1:ncol(d2)){
		yy=as.numeric(d2[,j])
		y=yy[!is.na(yy)]
		x=xx[!is.na(yy)]
		plot(x=x,y=y,col=NA,ylab=colnames(d2)[j],xlab='',cex.lab=1.6)
		points(x=x,y=y,col=col[j],pch=pch)
		abline(lm(y~x),col=col[j],lwd=2)
		rr=signif(cor(x,y,method='pearson'),2)
		pp=signif(cor.test(x,y,method='pearson')$p.value,2)
		l1=paste0("r=",rr,",p=",pp)
		mtext(side=3,text=l1,at=par('usr')[1],adj=0,col=col[j],cex=1.5,line=0.5)
	}
	
	plot(x=x,y=y,col=NA,ylab="",xlab="",bty='n',xaxt='n',yaxt='n')
	legend("left",bty='n',pch=c(3,8),legend=c('Gout','HUA'),cex=2,title=colnames(d1)[i])
	dev.off()

}







