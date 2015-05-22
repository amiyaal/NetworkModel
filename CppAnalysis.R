N=50

mean.deg=array(0,c(9,11))
sd.deg=array(0,c(9,11))
kurt=array(0,c(9,11))
skew=array(0,c(9,11))

den=array(0,c(9,11))
cc=array(0,c(9,11))

#model="BasicModel"
model="DeadByDegree"
for (i in 1:99){
  file=read.csv(paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/",model,"/summary",i,".csv"),sep=";")
  cc1=mean(file[,1])
  den1=mean(file[,2]/((N*(N-1))/2))
  pr=(i-1)%%9
  pn=(i-(pr+1))/9+1
  den[pr+1,pn]=den1
  cc[pr+1,pn]=cc1
  if (i>80) cat("pn=",pn," ")
}

prseq=c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5)
pnseq=seq(from=0,to=1,by=0.1)
library("gplots")
colfunc <- colorRampPalette(c("black", "red"))
heatmap.2(den,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Density",
          labCol=pnseq,labRow=prseq,xlab="Pn")
heatmap.2(cc,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="CC",
          labCol=pnseq,labRow=prseq,xlab="Pn")

s1=read.csv(paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/summary",1,".csv"),sep=";")
s2=read.csv(paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/summary",2,".csv"),sep=";")
s3=rbind(s1,s2)
s3$x=rep(1:50,2)
g=ggplot(s3, aes(x,cum_deg))
g+geom_point(aes(colour=factor(kill_by_degree)))+scale_y_log10()+scale_x_log10()+ggtitle("Pn=0.8 Pr=0.01")

#simulation analysis
for (i in 0:9){
  sim=read.csv(paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/simulation",i,".csv"),sep=";",header=F)
  plot(sim$V1,type="l")
}
all=read.csv(paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/all.csv"),sep=";",header=F)
all.mean=colSums(all)/100
all.mean=all.mean[1:100]
all.mean=t(all.mean)
all.mean=t(all.mean)
plot(all.mean,type="l")
