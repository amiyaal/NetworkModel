library("plyr")
library("igraph")
library("network")
library("pls")

parameters <- read.csv("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/1990/pnpr.csv", header=FALSE, sep=";", stringsAsFactors=FALSE)
colnames(parameters)=c("pn","pr","pi","img")
summary(parameters)

hyena1995network <- read.csv("~/Dropbox/Hyenas/Hyenas2/hyena1995network.csv", stringsAsFactors=FALSE)
hyena=as.matrix(hyena1995network)
N=dim(hyena)[1]
hyena=hyena[1:N,2:(N+1)]
hyena.net <- network(hyena, matrix.type="adjacency", directed=F)
plot(hyena.net)
hyena.inet=graph_from_adjacency_matrix(hyena,mode="undirected")
plot.igraph(hyena.inet, vertex.label=NA, vertex.size=5)
degs=igraph::degree(hyena.inet)
clus=igraph::transitivity(hyena.inet, type="local",isolates = "zero")
#verconn=vertex.connectivity(hyena.inet)
#diam=igraph::diameter(hyena.inet, directed=F)
laplace=graph.laplacian(hyena.inet)
eiglap=eigen(laplace, symmetric=T, only.values = T)$values
hyena.vec=c(sort(degs),sort(clus),eiglap)

hyena.vec89=hyena.vec
hyena.vec90=hyena.vec
hyena.vec91=hyena.vec
hyena.vec92=hyena.vec
hyena.vec93=hyena.vec
hyena.vec94=hyena.vec
hyena.vec95=hyena.vec
hyena.vec96=hyena.vec
hyena.vec97=hyena.vec
hyena.vec98=hyena.vec


predictor=matrix(data = 0,nrow=10000,ncol=N*3+4)
for (i in 1:10000){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/1998/Matrix",i,".csv")
  file=read.delim(filename,sep="\t",header = F,skip = 3,stringsAsFactors = F,col.names = c("V1","V2"))
  u1=unique(file$V1)
  u2=unique(file$V2)
  union=union(u1,u2)
  isolates=N-length(union)
  if (isolates > 0) union = c(union,1:isolates)
  v1=mapvalues(file$V1, from = union, to = 1:N,warn_missing = F)
  v2=mapvalues(file$V2, from = union, to = 1:N,warn_missing = F)
  file$V1=v1
  file$V2=v2
  #list <- strsplit(file$V1, "\t")
  #df <- ldply(file)
  if (isolates == N){
    net=network.initialize(N, directed=F)
  } else {
    net <- network(as.matrix(file), matrix.type="edgelist", directed=F)
    if (isolates > 0) network::add.vertices(net, isolates)
  }
  netm=as.matrix(net)
  inet=graph.adjacency(netm,mode="undirected")
  degs=igraph::degree(inet)
  clus=igraph::transitivity(inet, type="local",isolates = "zero")
  #verconn=vertex.connectivity(inet)
  #diam=igraph::diameter(inet, directed=F)
  laplace=graph.laplacian(inet)
  eiglap=eigen(laplace, symmetric=T, only.values = T)$values
  vec=c(sort(degs),sort(clus),parameters$pn[i],parameters$pr[i],parameters$pi[i],parameters$img[i],
    eiglap)
  predictor[i,] = vec
}
system("say 'done'")

dependent=predictor[,(N*2+1):(N*2+4)]
dependent2=dependent[,1:2]
predictor1=predictor[,-((N*2+1):(N*2+4))]
qpredictor=predictor1^2
predictor1=cbind(predictor1,qpredictor)
dat=as.data.frame(cbind(predictor1,dependent2))

dat90=dat
dat89=dat
dat91=dat
dat92=dat
dat93=dat
dat94=dat
dat95=dat
dat96=dat
dat97=dat
dat98=dat


xnam <- paste("V", 1:(N*6), sep="") #prepare formula
fmla <- as.formula(paste0("cbind(V",N*6+1,",V",N*6+2,") ~ ", paste(xnam, collapse= "+"))) #prepare formula
pls90=plsr(fmla, data=dat,ncomp=100,validation="CV")

summary(pls92)
plot(RMSEP(pls90), legendpos = "topright")
plot(pls90, ncomp = 58, asp = 1, line = TRUE,cex=0.25)

plot(pls91, plottype = "scores", comps = 1:3,cex=0.25)
explvar(pls92)

hy=as.data.frame(t(hyena.vec))
hy2=hy^2
hy3=cbind(hy,hy2,deparse.level =2)
colnames(hy3) = xnam

predict(pls90, ncomp = 56:60, newdata = hy3)

year=1989:1998
Ns=c(66,51,58,54,53,52,60,66,67,71)
Pn=c(0.77,.95,.51,.6,.72,.55,.42,.52,.57,.59)
Pr=c(.11,.035,.052,.041,.063,.032,.036,.058,.031,.047)
plot(Pn ~ year,type = "l",col="blue",ylim=c(0,1),lwd=3,ylab="")
lines(Pr ~ year,col="red",lwd=3)
legend ("topright",c("Pn","Pr"),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))


##############################################################################
# Confidence intervals for networks fits

nets=500
conf.deg=matrix(data = 0,nrow=nets,ncol=N)
conf.clus=matrix(data = 0,nrow=nets,ncol=N)
for (i in 1:nets){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/conf1995/Matrix",i,".csv")
#   file=read.delim(filename,sep="\t",header = F,skip = 3,stringsAsFactors = F,col.names = c("V1","V2"))
#   u1=unique(file$V1)
#   u2=unique(file$V2)
#   union=union(u1,u2)
#   isolates=N-length(union)
#   if (isolates > 0) cat(isolates," ")
#   if (isolates > 0) union = c(union,1:isolates)
#   v1=mapvalues(file$V1, from = union, to = 1:N,warn_missing = F)
#   v2=mapvalues(file$V2, from = union, to = 1:N,warn_missing = F)
#   file$V1=v1
#   file$V2=v2
#   #list <- strsplit(file$V1, "\t")
#   #df <- ldply(file)
#   if (isolates == N){
#     net=network.initialize(N, directed=F)
#   } else {
#     net <- network(as.matrix(file), matrix.type="edgelist", directed=F)
#     if (isolates > 0) network::add.vertices(net, isolates)
#   }
#   netm=as.matrix(net)
#   inet=graph.adjacency(netm,mode="undirected")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  degs=igraph::degree(inet)
  #plot(igraph::degree.distribution(inet,cumulative = T))
  #cum.degs=cumsum(sort(degs, decreasing = T))
  cum.degs=my.deg.dist(degs)
  clus=igraph::transitivity(inet, type="local",isolates = "zero")
  cum.clus=my.cc.dist(clus,cumulative = T)
  conf.deg[i,] = cum.degs
  conf.clus[i,] = cum.clus
  #conf.clus[i,] = sort(clus)
}
system("say 'done'")

#Degree dist
deg.res=apply(conf.deg, 2, function(x){quantile(x, c(0.025, 0.975))})
deg.ave=apply(conf.deg, 2, function(x){mean(x)})

cum.degs.98=my.deg.dist(degs)
cum.degs.95=my.deg.dist(degs)

dist.df=data.frame(x=1:N,upper=deg.res[2,],lower=deg.res[1,],ave=deg.ave,real=cum.degs.95)
dist.df[dist.df==0]=0.01
ggplot(subset(dist.df,x<=25),aes(x=x)) + geom_line(aes(y=ave)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=2)) +
  scale_y_log10() + theme_bw(base_size=15) + guides(size=FALSE) + 
  labs(y="Cumulative Probability",x="Degree")

#CC dist
cc.res=apply(conf.clus, 2, function(x){quantile(x, c(0.025, 0.975))})
cc.ave=apply(conf.clus, 2, function(x){mean(x)})

clus95=clus
cum.cc.98=my.cc.dist(clus98)/N
cum.cc.95=my.cc.dist(clus95)/N

cc.dist.df=data.frame(x=seq((1/N),1,1/N),upper=cc.res[2,],lower=cc.res[1,],ave=cc.ave,real=sort(clus98))
ggplot(subset(cc.dist.df,x<=1),aes(x=x)) + geom_line(aes(y=ave)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=2)) +
  theme_bw(base_size=15) + guides(size=FALSE) + 
  labs(y="Probability",x="Clustering Coefficient")

# cumulative:
cc.dist.df=data.frame(x=seq((1/N),1,1/N),upper=cc.res[2,]/N,lower=cc.res[1,]/N,ave=cc.ave/N,real=cum.cc.95)
#cc.dist.df[cc.dist.df==0]=0.001
ggplot(subset(cc.dist.df,x<=1),aes(x=x)) + geom_line(aes(y=ave)) +
geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=2)) +
scale_y_log10() + theme_bw(base_size=15) + guides(size=FALSE) +
labs(y="Cumulative Probability",x="Clustering Coefficient")


# deg in this function should be a series of degrees from a network
my.deg.dist=function (deg, cumulative = T, ...) 
{
  hi <- hist(deg, 0:length(deg), plot = FALSE)$density
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}

my.cc.dist=function (cc, cumulative = T, ...) 
{
  hi <- hist(cc, seq(0,1,1/length(cc)), plot = FALSE)$density
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}
#=================================================================
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
all.mean=t(all)
all.mean=t(all.mean)
plot(all.mean,type="l")

all=read.csv(paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/all.csv"),sep=";",header=F)
all.mean=t(all)
plot(all.mean,type="l")
