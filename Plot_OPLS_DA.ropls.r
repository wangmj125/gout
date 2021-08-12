arg <- commandArgs(T)

library('ropls')
library("ggplot2")
library("car")
library("rgl");
library("scatterplot3d");
library("plot3D")

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);



grpInfo <- read.table(groupFile,header=T,sep='\t',comment.char='',quote='');
rownames(grpInfo)=as.character(grpInfo[,1])
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
Data <- read.table(dataFile,header=T,sep="\t",row.names=1,comment.char='',quote='');

ii=intersect(rownames(grpInfo),colnames(Data))
if(length(ii)<nrow(grpInfo)){
	i1=setdiff(rownames(grpInfo),ii)
	print(paste0('Data Lost:',i1))
}

grpInfo=grpInfo[ii,]
Data=Data[,ii]


Data.tm <- t(as.matrix(Data)); 
groupname <-as.character(grpInfo[colnames(Data),2])


Groups <- unique(groupname)

if(length(Groups) == 2){
	Data.oplsda <- try(opls(Data.tm ,groupname, predI =1, orthoI = NA,plotL = FALSE, printL =FALSE),silent=TRUE);
	if ('try-error' %in% class(Data.oplsda))
	{
		Data.oplsda <- opls(Data.tm ,groupname, predI =1, orthoI = 2,plotL = FALSE, printL =FALSE);
	} 
	pdf(paste(outDir,"/",SamID,".OPLSDA.raw.pdf",sep=""),height=8,width=8)
	plot(Data.oplsda,typeVc="summary")
	plot(Data.oplsda,typeVc="x-score",parLabVc=groupname)
	VIP <- getVipVn(Data.oplsda);
	write.table(as.matrix(VIP),file=paste(outDir,"/",SamID,".All.OPLSDA.VIP.txt",sep=""),quote = F,sep="\t",col.names=NA,row.names=T);
	VIP <- VIP[VIP>1];
	write.table(as.matrix(VIP),file=paste(outDir,"/",SamID,".OPLSDA.VIP.txt",sep=""),quote = F,sep="\t",col.names=NA,row.names=T);
	
	dd <- attributes(Data.oplsda)
	dd.s <-dd$suppLs

	lab <-dd$summaryDF
	flag=dd.s$permMN

	xx=c(flag[,'sim'],flag[,'sim']) 
	yy=flag[,c(2,3)] 
	cc=c(rep('#008B8B',nrow(flag)),rep("#483D8B",nrow(flag))) 
	ll=paste0("pR2Y=",dd$summaryDF[7],',pQ2=',dd$summaryDF[8])

	pdf(paste(outDir,"/",SamID,".permutation.pdf",sep=""),height=8,width=8)
	par(mar=c(4,4,4,8))
	plot(x=xx,y=yy,col=cc,pch=16,xlab='',ylab='',main=ll,cex=1.5)
	abline(h=0)
	abline(v=0)
	legend("topright",bty='n',pch=16,col=c('#008B8B',"#483D8B"),legend=c('R2Y(Cum)','Q2(Cum)'),xpd=T,inset=-0.1)
	dev.off()
	
	
	
}


if(length(Groups) == 2){
	dd <- attributes(Data.oplsda)
	nn <- attributes(dd)
	dd.s <-dd$suppLs
	nn.s <- attributes(dd.s)
	lab <-dd$summaryDF
	p.info <- read.table(arg[5],header=T,row.names=1,comment.char='',quote='',sep='\t')
	p.nn <- rownames(dd$scoreMN)
	ii <- intersect(p.nn,rownames(p.info))
	p.info <- p.info[ii,]
	x <- dd$scoreMN
	y <- dd$orthoScoreMN[,1]	
	col <- as.character(p.info[,2])
	

	ff=paste(outDir,"/",SamID,".OPLSDA.point",sep="")
	xx1=as.numeric(arg[6])
	xx2=as.numeric(arg[7])
	yy1=as.numeric(arg[8])
	yy2=as.numeric(arg[9])


	ff=paste(outDir,"/",SamID,".OPLSDA.area",sep="")
	pdf(paste0(ff,".pdf"),w=8,h=8)
	p_pdf<-dev.cur()   
	png(paste0(ff,'.png'),w=8/7*480,h=8/7*480)
	dev.control("enable")
	
	par(mgp=c(2.5,1,0),mar=c(10,5,2,2))
	pp <-signif(dd$modelDF[1,1],2)*100
	plot(x=x,y=y,pch=rep(20,length(col)),col=col,main=paste(SamID," Scores(OPLS-DA)",sep=""),xlab=paste("t1(",pp,"%)",sep=""),ylab="to1",cex.lab=2,cex.main=2,cex.axis=2,xlim=c(xx1,xx2),ylim=c(yy1,yy2))

	
	dataEllipse(x=as.numeric(x),y=as.numeric(y),groups=factor(groupname,levels=unique(groupname),labels=unique(groupname)),levels = c(0.95), add=TRUE, col =unique(col) , lwd = 2,plot.points=FALSE,fill=TRUE,center.cex=0.2,group.labels="",fill.alpha = 0.3)

	legend("topright",pch=rep(20,length(unique(col))),col=unique(col),bty="n",legend=as.character(unlist(unique(p.info[,1]))),cex=1.8,pt.cex=1.8)

	legend( "bottom",horiz=T,xpd=T,inset=-0.3,legend=l3,bty="n",cex=2)

	dev.copy(which=p_pdf)  
	dev.off()
	dev.off()



}









