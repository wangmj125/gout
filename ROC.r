arg <- commandArgs(T)
library(ROCR)
library(pROC)
library(plotrix)


dataFile <- arg[1];
groupFile <- arg[2];
valFile<-arg[3]
valGroupFile<-arg[4]
outDir <- arg[5];

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];

data=read.table(dataFile ,sep='\t',header=T,row.names=1)
group=read.table(groupFile,sep='\t',header=T)
rownames(group)=as.character(group[,1])

g1=rep(0,nrow(group))
t=unique(as.character(group[,2]))
i1=grep(paste0('^',t[1],'$'),as.character(group[,2]))
g1[i1]=1

d.roc<-cbind(data,group=g1)

model1 <-glm(group~.,data=d.roc,family="binomial") 
pre <- predict(model1,type="response")
modelroc <- roc(d.roc$group,pre)
modelauc <-modelroc$auc
sen <- modelroc$sensitivities
spe <-modelroc$specificities

pdf(paste(outDir,"/",SamID,".Dis_roc.pdf",sep=""))
par(cex=1,mgp=c(3,0.5,0),lwd=2,cex.main=1,mar=c(4,4,4,4))
plot(modelroc, auc.polygon=F, grid=c(0.1, 0.2),grid.col=c("grey", "grey"), max.auc.polygon=F,auc.polygon.col="white",print.thres=F,col="red",print.auc.adj=c(0,-8),main=paste("2"," Metabolites"),xaxt="n",yaxt="n",xlab="",ylab="",cex.main=1.5)
axis(side=1,at=seq(0,1,length.out=3),line=0)
axis(side=2,at=seq(0,1,length.out=3),line=0,las=1)
mtext(side=1,adj=0.5,text="Specificity",cex=1.5,line=2.5,font=2)
mtext(side=2,adj=0.5,text="Sensitivity",cex=1.5,line=1.8,font=2)
p.auc <-round(modelauc,4)
p.ci <- round(as.numeric(ci.auc(modelauc))[c(1,3)],4)
t1 <- which.max(sen+spe)
p.sen <-round(sen[t1],4)*100
p.spe <-round(spe[t1],4)*100
l1 <- paste0("AUC = ",p.auc)
l2.1 <- paste0(p.ci,collapse="-")
l2 <- paste0("CI: ",l2.1)
l3 <- paste0("Sensitivity: ",p.sen,"%")
l4 <- paste0("Specificity: ",p.spe,"%")
legend("bottomright",legend=c(l1,l2,l3,l4),xpd=T,bty="n",cex=1.5)
dev.off()

d_val=read.table(valFile,sep='\t',header=T,row.names=1)
g_val=read.table(valGroupFile,sep='\t',header=T)
rownames(g_val)=as.character(g_val[,1])

g1=rep(0,nrow(g_val))
i1=grep(paste0('^',t[1],'$'),as.character(g_val[,2]))
g1[i1]=1

v.roc<-cbind(d_val,group=g1)


pre <- predict(model1,newdata=v.roc,type="response")
modelroc <- roc(v.roc$group,pre)
pre.dis <- modelroc$predictor

sen <- modelroc$sensitivities
spe <-modelroc$specificities
cutoff <- max(sen+spe)-1
if (cutoff ==1) cutoff=0.5


d1 <- pre[rownames(v.roc)[which(v.roc[,"group"]==1)]]
d2 <- pre[rownames(v.roc)[which(v.roc[,"group"]==0)]]

c.pred1 <- length(which(d1 >= cutoff))/length(d1)
c.pred2 <- length(which(d2<=cutoff))/length(d2)
c.pred<-c(c.pred1*100,c.pred2*100)
d.pdf <- vector(mode = "list", length = 2)
d.pdf[[1]] <-d1
d.pdf[[2]] <- d2
d.names <-t
names(d.pdf)<-d.names

pdf(paste(outDir,"/",SamID,".Val_cutoff.pdf",sep=""))
par(cex=1)
xx<-stripchart(d.pdf, vertical = TRUE,method = "jitter", pch = 16,col =c("#78BF45","#56AED2"))
abline(h=cutoff,col="red",lty=2,lwd=2,xpd=F)
text(x=1.5,y=cutoff+0.05,labels=paste("Cut off value = ",round(cutoff,3),sep=""),cex=1.5)	
mtext(side=1,line = 2,adj=0.5,at=seq(1,length(d.names),1),text=paste0("n = ",lengths(d.pdf[c(1:2)])),cex=1) 	
mtext(side=3,line = 0,adj=0.5,at=seq(1,length(d.names),1),text=paste0(round(as.numeric(c.pred[1:2]),2),"%"),cex=1) 
mtext(side=3,line = 0,adj=0,at=0,text="Pred. :",cex=1,xpd=T) 
dev.off()

