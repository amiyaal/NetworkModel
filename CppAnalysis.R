library("plyr")
library("igraph")
library("network")
library("pls")
library("ggplot2")

lizard=read.csv(paste0("~/Dropbox/Erol/Published networks/lizard.csv"),sep=",",header=F)
lizard=as.matrix(lizard)
quantile(lizard[lizard>0], na.rm=T)
thre=0.0089 #threshold for binarization
lizard[lizard<thre]=0
lizard[lizard>=thre]=1
hyena=lizard
n4=graph_from_adjacency_matrix(lizard,mode="undirected")
n4=simplify(n4)
plot.igraph(n4, vertex.label=NA, vertex.size=5, edge.color="black")
c.lizard=cluster_walktrap(n4)
modularity(c.lizard)
plot(c.lizard,n4,vertex.label=NA, vertex.size=5, edge.color="black")
graph.density(n4)
transitivity(n4, type="global")
clus= igraph::transitivity(n4, type="local",isolates = "zero")
  
dolphin=read.csv(paste0("~/Dropbox/Erol/Published networks/dolphins/out.dolphins"),sep="\t",header=F)
g<-network.initialize(62, directed=F)
network.edgelist(dolphin, g)
n3=as.matrix(g)
n3=graph_from_adjacency_matrix(n3,mode="undirected")
n3=simplify(n3)
plot.igraph(n3, vertex.label=NA, vertex.size=5, edge.color="black")
hyena=as.matrix(g)
c.dolphin=cluster_walktrap(n3)
modularity(c.dolphin)
plot(c.dolphin,n3,vertex.label=NA, vertex.size=5, edge.color="black")
graph.density(n3)
transitivity(n3, type="global")
clus=igraph::transitivity(n3, type="local",isolates = "zero")

filename="~/Dropbox/Erol/Published networks/David 2009.csv"
David2009=read.csv(filename, stringsAsFactors = F)
groups=strsplit(David2009$ID," ")
dates=as.list(David2009$Date)
##getgroupbyindividualmatrix 
gbiDavid2009 = get_group_by_individual(groups, data_format = "groups")
networkDavid2009 = get_network(gbiDavid2009, data_format= "GBI",times=dates)
quantile(networkDavid2009[networkDavid2009>0],0.25)
hyrax=networkDavid2009
thre=0.158 #threshold for binarization :25% precentile
hyrax[hyrax<thre]=0
hyrax[hyrax>=thre]=1
hyrax.inet=graph_from_adjacency_matrix(hyrax,mode="undirected")

hyraxD2009network <- read.csv("~/Downloads/hyrax networks/David 2009.csv", header=FALSE, stringsAsFactors=FALSE)
hyrax=as.matrix(hyraxD2009network)
N=dim(hyrax)[1]
hyrax.inet=graph_from_adjacency_matrix(hyrax,mode="undirected")
hyrax.inet=simplify(hyrax.inet)
plot.igraph(hyrax.inet, vertex.label=NA, vertex.size=5, edge.color="black")
degs=igraph::degree(hyrax.inet)
clus=igraph::transitivity(hyrax.inet, type="local",isolates = "zero")

thre=0.064
hyena1997network <- read.csv("~/Copy/Hyenas/Hyenas2/hyena1997network.csv", stringsAsFactors=FALSE)
hyena=as.matrix(hyena1997network)
hyraxD2009network <- read.csv("~/Downloads/hyrax networks/David 2009.csv", header=FALSE, stringsAsFactors=FALSE)
hyena=as.matrix(hyraxD2009network)
N=dim(hyena)[1]
hyena=hyena[1:N,2:(N+1)]
class(hyena) <- "numeric"
hyena[hyena<thre]=0
hyena[hyena>0]=1
#hyena.net <- network(hyena, matrix.type="adjacency", directed=F)
#plot(hyena.net)
hyena.inet=graph_from_adjacency_matrix(hyena,mode="undirected")
hyena.inet=simplify(hyena.inet)
hyrax.inet=hyena.inet
plot.igraph(hyena.inet, vertex.label=NA, vertex.size=5, edge.color="black")
degs=igraph::degree(hyena.inet)
clus=igraph::transitivity(hyena.inet, type="local",isolates = "zero")
#verconn=vertex.connectivity(hyena.inet)
#diam=igraph::diameter(hyena.inet, directed=F)
#laplace=graph.laplacian(hyena.inet)
#eiglap=eigen(laplace, symmetric=T, only.values = T)$values
#hyena.vec=c(sort(degs),sort(clus),eiglap)
hyena.vec=c(sort(degs),sort(clus))
c.hyena=cluster_walktrap(hyena.inet)
modularity(c.hyena)
plot(c.hyena,hyena.inet,vertex.label=NA, vertex.size=5, edge.color="black")
c.hyrax=cluster_walktrap(hyrax.inet)
modularity(c.hyrax)
plot(c.hyrax,hyrax.inet,vertex.label=NA, vertex.size=5, edge.color="black")
graph.density(hyrax.inet)
graph.density(hyena.inet)
transitivity(hyrax.inet, type="global")
transitivity(hyena.inet, type="global")

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
hyraxD09=hyena.vec
dolphinLusseau=hyena.vec
lizard.vec=hyena.vec

parameters <- read.csv("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/hyena_pls3/pnpr.csv", header=FALSE, sep=";", stringsAsFactors=FALSE)
#parameters <- read.csv("/Volumes/HDD/Documents/Social Inheritance/Data/1997/pnpr.csv", header=FALSE, sep=";", stringsAsFactors=FALSE)
colnames(parameters)=c("pn","pr")
summary(parameters)

predictor=matrix(data = 0,nrow=10000,ncol=N*2+2)
for (i in 1:10000){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/hyena_pls3/Matrix",i,".csv")
  #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/1997/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  degs=igraph::degree(inet)
  clus=igraph::transitivity(inet, type="local",isolates = "zero")
  #verconn=vertex.connectivity(inet)
  #diam=igraph::diameter(inet, directed=F)
  #laplace=graph.laplacian(inet)
  #eiglap=eigen(laplace, symmetric=T, only.values = T)$values
  #vec=c(sort(degs),sort(clus),parameters$pn[i],parameters$pr[i],parameters$pi[i],parameters$img[i],
  #  eiglap)
  vec=c(sort(degs),sort(clus),parameters$pn[i],parameters$pr[i])
  predictor[i,] = vec
}
system("say 'done'")

dependent=predictor[,(N*2+1):(N*2+2)]
predictor1=predictor[,-((N*2+1):(N*2+2))]
#qpredictor=predictor1^2
#predictor1=cbind(predictor1,qpredictor)
dat=as.data.frame(cbind(predictor1,dependent))

dat.lizard=dat
dat.dolphinLusseau=dat
dat.hyrax.D09=dat
dat89=dat
dat90=dat
dat91=dat
dat92=dat
dat93=dat
dat94=dat
dat95=dat
dat96=dat
dat97=dat
dat98=dat

xnam <- paste("V", 1:(N*2), sep="") #prepare formula
fmla <- as.formula(paste0("cbind(V",N*2+1,",V",N*2+2,") ~ ", paste(xnam, collapse= "+"))) #prepare formula
pls.test3=plsr(fmla, data=dat,ncomp=68,validation="CV")

#summary(pls.test)
plot(RMSEP(pls.test3), legendpos = "topright")
fitted=plot(pls.test3, ncomp = 32, asp = 1, line = TRUE,cex=0.25)

plot(parameters$pn,pls.test2$fitted.values[,1,25])
x=parameters$pn
y=pls.test2$fitted.values[,1,25]
lm=lm(y ~ poly(x, 2, raw=TRUE))
lm$coefficients
plot(lm)

#plot(pls91, plottype = "scores", comps = 1:3,cex=0.25)
#explvar(pls92)

hy=as.data.frame(t(hyena.vec))
#hy2=hy^2
#hy3=cbind(hy,hy2,deparse.level =2)
colnames(hy) = xnam
predict(pls.test, ncomp = 10:14, newdata = hy)

########################### 
#Verification of PLS
nets=500
col.res=matrix(data = 0,nrow=nets,ncol=2)
for (i in 1:nets){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/pls_simulate/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  degs=igraph::degree(inet)
  clus=igraph::transitivity(inet, type="local",isolates = "zero")
  hyena.vec=c(sort(degs),sort(clus))
  hy=as.data.frame(t(hyena.vec))
  colnames(hy) = xnam
  col=predict(pls.test, ncomp = 24, newdata = hy)
  col.res[i,1]=col[1]
  col.res[i,2]=col[2]
}
summary(col.res[,1])
col.res.df=as.data.frame(col.res)
colnames(col.res.df)=c("Pn","Pr")

gs1=ggplot(col.res.df, aes(x=Pn)) + geom_histogram(binwidth=.008, colour="black", fill="white") +
  geom_vline(aes(xintercept=0.82), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Pn")
ggsave(file="pn_dist.pdf")

gs2=ggplot(col.res.df, aes(x=Pr)) +geom_histogram(binwidth=.002, colour="black", fill="white") +
  geom_vline(aes(xintercept=0.014), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Pr")

pdf(file="combined_dist.pdf",width = 8, height = 5)
grid.arrange(gs1,gs2,ncol=2)
dev.off()

########################### 
# Verification of PLS for 10 values of Pn
nets=500
all.dist=matrix(data = 0,nrow=nets*10,ncol=3)
for (pn in 1:10){
  for (i in 1:nets){
    filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/simulation",pn,"/Matrix",i,".csv")
    inet=read_graph(filename,format = "pajek")
    inet=simplify(inet)
    degs=igraph::degree(inet)
    clus=igraph::transitivity(inet, type="local",isolates = "zero")
    hyena.vec=c(sort(degs),sort(clus))
    hy=as.data.frame(t(hyena.vec))
    colnames(hy) = xnam
    col=predict(pls.test3, ncomp = 24, newdata = hy)
    all.dist[i+(pn-1)*nets,1]=col[1]
    all.dist[i+(pn-1)*nets,2]=col[2]
    all.dist[i+(pn-1)*nets,3]=pn
  }
}
system("say 'ok'")
summary(all.dist[1:500,3])
hist(all.dist[4501:5000,1])
hist(all.dist[1:500,1])
all.dist.df=as.data.frame(all.dist)
colnames(all.dist.df)=c("Pn","Pr","Simulated")

pdf("pn_all_sim.pdf")
val=seq(0.6,0.9,length.out = 10)
i=1:10
hg=ggplot(all.dist.df, aes(factor(Simulated),Pn)) + geom_boxplot() + theme_bw(base_size=20) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  guides(size=FALSE) + labs(y="Predicted Pn",x="Simulated Pn")
hg+scale_x_discrete(breaks=seq(1:10),labels=round(val,digits = 2))+
  geom_segment(aes(x = i-0.4, y = val[i], xend = i+0.4, yend = val[i]),colour="red")
dev.off()

#############
year=1989:1998
Ns=c(66,50,57,51,50,47,53,58,58,57)
Pn=c(0.94,.93,.81,.8,.87,.79,.64,.82,.82,.78)
Pr=c(.031,.024,.025,.028,.023,.016,.008,.015,.014,.011)
plot(Pn ~ year,type = "l",col="blue",ylim=c(0,1),lwd=3,ylab="")
lines(Pr ~ year,col="red",lwd=3)
legend ("topright",c("Pn","Pr"),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))


##############################################################################
# Confidence intervals for networks fits

N=60 # hyena=58, hyrax=34, dolphin=62, lizard=60
nets=500
conf.deg=matrix(data = 0,nrow=nets,ncol=N)
conf.clus=matrix(data = 0,nrow=nets,ncol=N)
conf.mod=matrix(data = 0,nrow=nets,ncol=1) # stores modularity values
for (i in 1:nets){
  #filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/hyena_pls/Matrix",i,".csv")
  filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  degs=igraph::degree(inet)
  #plot(igraph::degree.distribution(inet,cumulative = T))
  #cum.degs=cumsum(sort(degs, decreasing = T))
  cum.degs=my.deg.dist(degs)
  clus=igraph::transitivity(inet, type="local",isolates = "zero")
  cum.clus=my.cc.dist(clus,cumulative = T)
  mod=cluster_walktrap(inet) %>% modularity()
  conf.deg[i,] = cum.degs
  conf.clus[i,] = cum.clus
  conf.mod[i,1] = mod
  #conf.clus[i,] = sort(clus)
}
system("say 'done'")

#Modularity
conf.mod.dolphin=conf.mod
conf.mod.hyena=conf.mod
conf.mod.hyrax=conf.mod
conf.mod.lizard=conf.mod
conf.mod1=data.frame(conf.mod)
md1=ggplot(conf.mod1, aes(x=conf.mod)) + geom_histogram(binwidth=.05, colour="black", fill="white") +
  geom_vline(aes(xintercept=0.23,y=0), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Modularity") + ggtitle('Spotted hyena')
md2=ggplot(conf.mod1, aes(x=conf.mod)) + geom_histogram(binwidth=.05, colour="black", fill="white") +
  geom_vline(aes(xintercept=0.40,y=0), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Modularity") + ggtitle('Rock hyrax')
md3=ggplot(conf.mod1, aes(x=conf.mod)) + geom_histogram(binwidth=.02, colour="black", fill="white") +
  geom_vline(aes(xintercept=0.49,y=0), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Modularity") + ggtitle('Bottlenose dolphin')
md4=ggplot(conf.mod1, aes(x=conf.mod)) + geom_histogram(binwidth=.03, colour="black", fill="white") +
  geom_vline(aes(xintercept=0.80,y=0), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Modularity") + ggtitle('Sleepy lizard')

require(gridExtra)
pdf(file="combined_modularity.pdf",width = 8, height = 5)
grid.arrange(md1,md2,md3,md4,ncol=2)
dev.off()


#Degree dist
deg.res=apply(conf.deg, 2, function(x){quantile(x, c(0.025, 0.975))})
deg.ave=apply(conf.deg, 2, function(x){mean(x)})

cum.degs.97=my.deg.dist(degs)
cum.degs.98=my.deg.dist(degs)
cum.degs.95=my.deg.dist(degs)
cum.degs.hyrax2=my.deg.dist(degs)
cum.degs.dolphin=my.deg.dist(degs)
cum.degs.lizard=my.deg.dist(degs)

dolphinmf=rev(cumsum(rev(dolphinMeanfield$V2)))
lizardmf=rev(cumsum(rev(lizardMeanfield$V2)))
hyenamf=rev(cumsum(rev(hyena97Meanfield$V2)))
hyraxmf=rev(cumsum(rev(hyraxd09Meanfield$V2)))

hyena.dist.df=data.frame(x=0:(N-1),upper=deg.res[2,],lower=deg.res[1,],ave=deg.ave,real=cum.degs.97)
hyrax.dist.df=data.frame(x=0:(N-1),upper=deg.res[2,],lower=deg.res[1,],ave=deg.ave,real=cum.degs.hyrax2)
dolphin.dist.df=data.frame(x=0:(N-1),upper=deg.res[2,],lower=deg.res[1,],ave=deg.ave,real=cum.degs.dolphin,mf=dolphinmf)
lizard.dist.df=data.frame(x=0:(N-1),upper=deg.res[2,],lower=deg.res[1,],ave=deg.ave,real=cum.degs.lizard,mf=lizardmf)

update_geom_defaults("point",   list(shape = 21))

# hyena degree dist.
dist.df=hyena.dist.df
hyena.dist.df$mf=hyenamf
hyena.dist.df$x=0:57
hyena.dist.df[hyena.dist.df==0]=0.0001
pdf(file="p1.pdf", width = 25/4,height = 10/2)
g1=ggplot(subset(hyena.dist.df,x<=30),aes(x=x)) + geom_line(aes(y=ave,size=3)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  geom_point(shape=22,aes(y=mf,size=4),colour="blue") + scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) + 
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

# hyrax degree dist.
dist.df=hyrax.dist.df
hyrax.dist.df$mf=hyraxmf
hyrax.dist.df$x=0:33
hyrax.dist.df[hyrax.dist.df==0]=0.001
pdf(file="p2.pdf", width = 25/4,height = 10/2)
gh1=ggplot(subset(hyrax.dist.df,x<=18),aes(x=x)) + geom_line(aes(y=ave,size=3)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  geom_point(shape=22,aes(y=mf,size=4),colour="blue") + scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) + 
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

# dolphin degree dist.
dolphin.dist.df[dolphin.dist.df==0]=0.0001
pdf(file="p3.pdf", width = 25/4,height = 10/2)
gd1=ggplot(subset(dolphin.dist.df,x<=17),aes(x=x)) + geom_line(aes(y=ave,size=3)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  geom_point(shape=22,aes(y=mf,size=4),colour="blue") + scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) + 
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

# lizard degree dist.
lizard.dist.df[lizard.dist.df==0]=0.0001
pdf(file="p4.pdf", width = 25/4,height = 10/2)
gl1=ggplot(subset(lizard.dist.df,x<=13),aes(x=x)) + geom_line(aes(y=ave,size=3)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  geom_point(shape=22,aes(y=mf,size=4),colour="blue") + scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) + 
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

hyena.dist.df=dist.df
hyrax.dist.df=dist.df
dolphin.dist.df=dist.df
lizard.dist.df=dist.df

#CC dist
cc.res=apply(conf.clus, 2, function(x){quantile(x, c(0.025, 0.975))})
cc.ave=apply(conf.clus, 2, function(x){mean(x)})

clus97=clus
clus.hyrax=clus
cum.cc.98=my.cc.dist(clus98)/N
cum.cc.97=my.cc.dist(clus97)/N
cum.cc.95=my.cc.dist(clus95)/N
cum.cc.hyrax=my.cc.dist(clus.hyrax)/N
cum.cc.dolphin=my.cc.dist(clus)/N
cum.cc.lizard=my.cc.dist(clus)/N

# to go back to filled points: replace shape=21 with alpha=0.7
#hyena
hyena.cc.dist.df=data.frame(x=seq((1/N),1,1/N),upper=cc.res[2,]/N,lower=cc.res[1,]/N,ave=cc.ave/N,real=cum.cc.97)
hyena.cc.dist.df[hyena.cc.dist.df==0]=0.001
odd_indexes<-seq(1,length(hyena.cc.dist.df$x),2)

pdf(file="p5.pdf", width = 25/4,height = 10/2)
g2=ggplot(subset(hyena.cc.dist.df[odd_indexes,],x<=1),aes(x=x)) + geom_line(aes(y=ave,size=3)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + 
  geom_point(aes(y=real,size=4)) +
  scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) + 
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

# hyrax
cc.dist.df=data.frame(x=seq((1/N),1,1/N),upper=cc.res[2,]/N,lower=cc.res[1,]/N,ave=cc.ave/N,real=cum.cc.hyrax)
cc.dist.df[cc.dist.df==0]=0.001
odd_indexes<-seq(1,length(cc.dist.df$x),2)
pdf(file="p6.pdf", width = 25/4,height = 10/2)
gh2=ggplot(subset(cc.dist.df[odd_indexes,],x<=1),aes(x=x)) + geom_line(aes(y=ave,size=3)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) +
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

# dolphin
dolphin.cc.dist.df=data.frame(x=seq((1/N),1,1/N),upper=cc.res[2,]/N,lower=cc.res[1,]/N,ave=cc.ave/N,real=cum.cc.dolphin)
dolphin.cc.dist.df[dolphin.cc.dist.df==0]=0.01
odd_indexes<-seq(1,length(dolphin.cc.dist.df$x),2)
pdf(file="p7.pdf", width = 25/4,height = 10/2)
gd2=ggplot(subset(dolphin.cc.dist.df[odd_indexes,],x<=1),aes(x=x)) + geom_line(aes(y=ave,size=3)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) +
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

#lizard
lizard.cc.dist.df=data.frame(x=seq((1/N),1,1/N),upper=cc.res[2,]/N,lower=cc.res[1,]/N,ave=cc.ave/N,real=cum.cc.lizard)
lizard.cc.dist.df[lizard.cc.dist.df==0]=0.001
odd_indexes<-seq(1,length(lizard.cc.dist.df$x),2)
pdf(file="p8.pdf", width = 25/4,height = 10/2)
gl2=ggplot(subset(lizard.cc.dist.df[odd_indexes,],x<=1),aes(x=x)) + geom_line(aes(y=ave,size=3)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",alpha=0.4) + geom_point(aes(y=real,size=4)) +
  scale_y_log10() + theme_bw(base_size=15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=20)) + 
  guides(size=FALSE) + labs(y="",x="")
grid.edit("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

cc.dist.df=hyrax.cc.dist.df
hyena.cc.dist.df=cc.dist.df
hyrax.cc.dist.df=cc.dist.df
dolphin.cc.dist.df=cc.dist.df
lizard.cc.dist.df=cc.dist.df

# Model fit figure
require(gridExtra)
pdf(file="FitLower2.pdf", width = 25,height = 10)
eight=grid.arrange(g1, gh1, gd1, gl1, g2, gh2, gd2, gl2,ncol=4,name="all")
all.grob=arrangeGrob(g1, gh1, gd1, gl1, g2, gh2, gd2, gl2,ncol=4,name="all")
all.grob = ggplotGrob(eight)
grid.ls(fullNames = T)
grid.draw(all.grob)
grid.edit(gPath("all.grob", "geom_point.points"), grep = T, gp = gpar(lwd = 4),global = T)
editGrob("geom_point.points", grep = T, gp = gpar(lwd = 4),global = T)
dev.off()

par(mfrow=c(1,4))
par(mai=c(.5,.5,.5,.5))
plot.igraph(hyena.inet, vertex.label=NA, vertex.size=7, edge.color="black")
plot.igraph(hyrax.inet, vertex.label=NA, vertex.size=7, edge.color="black")
plot.igraph(n3, vertex.label=NA, vertex.size=7, edge.color="black") #dolphin
plot.igraph(n4, vertex.label=NA, vertex.size=7, edge.color="black", margin=0) #lizard

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# deg in this function should be a series of degrees from a network
my.deg.dist=function (deg, cumulative = T, ...) 
{
  hi <- hist(deg, 0:length(deg), right=F,plot = FALSE)$density
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
  hi <- hist(cc, seq(0,1,1/length(cc)), right=F,plot = FALSE)$density
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}

###############
# Heritability figure
regressdata <- read.csv("~/Dropbox/Erol/NetworkModel/regressdata.csv", stringsAsFactors=FALSE)
hg=ggplot(regressdata, aes(factor(pn),regression.coefficient)) + geom_boxplot() + theme_bw(base_size=20) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15),axis.title.x = element_text(face="italic")) + 
  guides(size=FALSE) + labs(y="Parent-offspring regression coefficient",x="Pn")
hg+scale_x_discrete(breaks=seq(.1,.9,.2),labels=seq(.1,.9,.2))

##################################
# Assortativity
nets=500
assort=matrix(data = 0,nrow=nets,ncol=9)
assort.random=matrix(data = 0,nrow=nets,ncol=9)
for (j in 1:9){
  for (i in 1:nets){
    filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/assort",j,"/Matrix",i,".csv")
    inet=read_graph(filename,format = "pajek")
    inet=simplify(inet)
    V(inet)$name=V(inet)$id
    traitfile=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/assort",j,"/Traits",i,".csv")
    traits <- read.csv(traitfile, header=FALSE, sep=";", stringsAsFactors=FALSE)
    inet=set.vertex.attribute(inet, "trait", index=as.character(traits$V1), value=traits$V2)
    assort[i,j]=igraph::assortativity(inet, types1=V(inet)$trait, directed = F)
    # Now shuffling the traits assigned
    inet=set.vertex.attribute(inet, "trait", index=as.character(traits$V1), value=sample(traits$V2,100))
    assort.random[i,j]=igraph::assortativity(inet, types1=V(inet)$trait, directed = F)
  }
}
system("say 'done'")
summary(assort)
summary(assort.random)
which.max(assort[,4])

as.df <- data.frame(matrix(ncol = 4, nrow = 9))
colnames(as.df)=c("model.means","model.error","null.means","null.error")
for (i in 1:9){
  as.df$model.means[i]=mean(assort[,i])
  as.df$model.error[i]=std.error(assort[,i])
  as.df$null.means[i]=mean(assort.random[,i])
  as.df$null.error[i]=std.error(assort.random[,i])
}
#df.as=data.frame(model.means,model.sd,null.means,null.sd)

# Assortativity figure:
gg=ggplot(as.df, aes(x=seq(0.1,0.9,0.1))) +
  geom_line(aes(y=null.means)) + geom_point(aes(y=null.means, size=3)) +
  geom_errorbar(aes(ymin=model.means-model.error, ymax=model.means+model.error), colour="red", width=.01) +
  geom_line(aes(y=model.means),colour="red") + geom_point(aes(y=model.means), colour="red",size=3) +
  geom_errorbar(aes(ymin=null.means-null.error, ymax=null.means+null.error), colour="black", width=.01) +
  theme_bw(base_size=15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15),axis.title.x = element_text(face="italic")) + guides(size=FALSE) +
  labs(y="Assortativity Coefficient",x="Pn")
gg+ scale_x_continuous(breaks=seq(0.1,0.9,by=0.1))

trait.col=as.numeric(factor(V(inet)$trait),levels=1:100)
cols=viridis(100)
plot(inet, vertex.color=cols[trait.col], vertex.label=NA, vertex.size=10, edge.color="black")


##################################
# Explicit Assortativity Model
nets=500
explicit.assort=matrix(data = 0,nrow=nets,ncol=1)
explicit.assort.random=matrix(data = 0,nrow=nets,ncol=1)
explicit.clus=matrix(data = 0,nrow=nets,ncol=1)
explicit.mod=matrix(data = 0,nrow=nets,ncol=1)
for (i in 1:nets){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/ExplicitAssort/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  V(inet)$name=V(inet)$id
  traitfile=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/ExplicitAssort/Traits",i,".csv")
  traits <- read.csv(traitfile, header=FALSE, sep=";", stringsAsFactors=FALSE)
  inet=set.vertex.attribute(inet, "trait", index=as.character(traits$V1), value=traits$V2)
  explicit.assort[i,1]=igraph::assortativity(inet, types1=V(inet)$trait, directed = F)
  # Now shuffling the traits assigned
  inet=set.vertex.attribute(inet, "trait", index=as.character(traits$V1), value=sample(traits$V2,100))
  explicit.assort.random[i,1]=igraph::assortativity(inet, types1=V(inet)$trait, directed = F)
  explicit.clus[i,1]=igraph::transitivity(inet, type="global",isolates = "zero")
  explicit.mod[i,1]=cluster_walktrap(inet) %>% modularity()
}
system("say 'done'")
summary(explicit.assort)
summary(explicit.assort.random)
summary(explicit.clus)
summary(explicit.mod)

inet=set.vertex.attribute(inet, "trait", index=as.character(traits$V1), value=traits$V2)
inet=delete_edges(inet,sample(E(inet),1173))
graph.density(inet)
igraph::assortativity(inet, types1=V(inet)$trait, directed = F)
trait.col=as.numeric(factor(V(inet)$trait),levels=1:100)
cols=viridis(100)
pdf(file="ExplicitAssortExample.pdf", width = 5,height = 5)
plot(inet, vertex.color=cols[trait.col], vertex.label=NA, vertex.size=10, edge.color="black")
dev.off()
curve(exp(-x)-exp(-1))

explicit.clus1=data.frame(explict.clus)
exp1=ggplot(explicit.clus1, aes(x=explicit.clus)) + geom_histogram(binwidth=.02, colour="black", fill="white") +
  geom_vline(aes(xintercept=c(0.31,0.49,0.62,0.94)), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Global Clustering Coefficent") + 
  annotate("text", x=0.64, y=200, label="Hyena", color="red",angle = 270) +
  annotate("text", x=0.96, y=200, label="Hyrax", color="red", angle = 270) +
  annotate("text", x=0.33, y=200, label="Dolphin", color="red",angle = 270) +
  annotate("text", x=0.51, y=200, label="Lizard", color="red",angle = 270)

explicit.mod1=data.frame(explicit.mod)
exp2=ggplot(explicit.mod1, aes(x=explicit.mod)) + geom_histogram(binwidth=.02, colour="black", fill="white") +
  geom_vline(aes(xintercept=c(0.23,0.4,0.49,0.8)), color="red", linetype="dashed", size=1) +
  theme_bw(base_size = 15) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.text = element_text(size=15)) + 
  labs(y="Count",x="Modularity") + 
  annotate("text", x=0.25, y=70, label="Hyena", color="red",angle = 270) +
  annotate("text", x=0.42, y=70, label="Hyrax", color="red", angle = 270) +
  annotate("text", x=0.51, y=70, label="Dolphin", color="red",angle = 270) +
  annotate("text", x=0.82, y=70, label="Lizard", color="red",angle = 270)

pdf("ExplicitAssortComparison.pdf",width = 12,height = 5)
grid.arrange(exp1,exp2,ncol=2)
dev.off()
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

#########################
# Unstable population size

# Shrinking network size
nets=500
shrink.dens=matrix(data = 0,nrow=nets,ncol=1)
shrink.clus=matrix(data = 0,nrow=nets,ncol=1)
shrink.mod=matrix(data = 0,nrow=nets,ncol=1) # stores modularity values
Ns=matrix(data = 0,nrow=nets,ncol=1) # Store final population size
for (i in 1:nets){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/LizardShrink/Matrix",i,".csv")
  #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  shrink.dens[i,1]=graph.density(inet)
  shrink.clus[i,1]=igraph::transitivity(inet, type="global",isolates = "zero")
  shrink.mod[i,1] = cluster_walktrap(inet) %>% modularity()
  Ns[i,1]=vcount(inet)
}
system("say 'done'")

summary(Ns)
hist(shrink.mod)
hist(shrink.dens)
hist(shrink.clus)

# Stable network size
stable.dens=matrix(data = 0,nrow=nets,ncol=1)
stable.clus=matrix(data = 0,nrow=nets,ncol=1)
stable.mod=matrix(data = 0,nrow=nets,ncol=1) # stores modularity values
Ns=matrix(data = 0,nrow=nets,ncol=1) # Store final population size
for (i in 1:nets){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/LizardStable2/Matrix",i,".csv")
  #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  stable.dens[i,1]=graph.density(inet)
  stable.clus[i,1]=igraph::transitivity(inet, type="global",isolates = "zero")
  stable.mod[i,1] = cluster_walktrap(inet) %>% modularity()
  Ns[i,1]=vcount(inet)
}
system("say 'done'")

# Growing network size
nets=500
grow.dens=matrix(data = 0,nrow=nets,ncol=1)
grow.clus=matrix(data = 0,nrow=nets,ncol=1)
grow.mod=matrix(data = 0,nrow=nets,ncol=1) # stores modularity values
Ns=matrix(data = 0,nrow=nets,ncol=1) # Store final population size
for (i in 1:nets){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/LizardGrow/Matrix",i,".csv")
  #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  grow.dens[i,1]=graph.density(inet)
  grow.clus[i,1]=igraph::transitivity(inet, type="global",isolates = "zero")
  grow.mod[i,1] = cluster_walktrap(inet) %>% modularity()
  Ns[i,1]=vcount(inet)
}
system("say 'done'")

summary(Ns)
hist(stable.mod)
hist(stable.dens)
hist(stable.clus)
  
mean(stable.dens)
sd(stable.dens)
mean(shrink.dens)
sd(shrink.dens)

t.test(stable.mod,shrink.mod)
t.test(stable.dens,shrink.dens)
t.test(stable.clus,shrink.clus)

mean(stable.clus)
sd(stable.clus)
mean(shrink.clus)
sd(shrink.clus)

mean(stable.mod)
sd(stable.mod)
mean(shrink.mod)
sd(shrink.mod)

#===============
mean(stable.dens)
sd(stable.dens)
mean(grow.dens)
sd(grow.dens)

t.test(stable.mod,grow.mod)
t.test(stable.dens,grow.dens)
t.test(stable.clus,grow.clus)

mean(stable.clus)
sd(stable.clus)
mean(grow.clus)
sd(grow.clus)

mean(stable.mod)
sd(stable.mod)
mean(grow.mod)
sd(grow.mod)
#==========================================
# New parameter set - shrink & growth - Oct 5 2015
nets=100
#shrink
shrink.densities=data.frame(matrix(ncol = 5, nrow = 5*3*nets))
shrink.ccs=data.frame(matrix(ncol = 5, nrow = 5*3*nets))
shrink.modularities=data.frame(matrix(ncol = 5, nrow = 5*3*nets))
colnames(shrink.densities)=c("Pn","Pr","i","value","type")
colnames(shrink.ccs)=c("Pn","Pr","i","value","type")
colnames(shrink.modularities)=c("Pn","Pr","i","value","type")
loc=0
for (a in 1:5){
  for (b in 1:3){
    #shrink.dens=matrix(data = 0,nrow=nets,ncol=1)
    #shrink.clus=matrix(data = 0,nrow=nets,ncol=1)
    #shrink.mod=matrix(data = 0,nrow=nets,ncol=1) # stores modularity values
    for (i in 1:nets){
      filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/ShrinkSim/Matrix",(a-1)*3+b,"A",
                      i,".csv")
      #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
      inet=read_graph(filename,format = "pajek")
      inet=simplify(inet)
      loc=loc+1
      shrink.densities[loc,]=c(a,b,i,graph.density(inet),1)
      shrink.ccs[loc,]=c(a,b,i,igraph::transitivity(inet, type="global",isolates = "zero"),1)
      shrink.modularities[loc,]=c(a,b,i,cluster_walktrap(inet) %>% modularity(),1)
    }
  }
}

for (a in 1:5){
  for (b in 1:3){
    for (i in 1:nets){
      filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/ShrinkSimStable/Matrix",(a-1)*3+b,"A",
                      i,".csv")
      #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
      inet=read_graph(filename,format = "pajek")
      inet=simplify(inet)
      loc=loc+1
      shrink.densities[loc,]=c(a,b,i,graph.density(inet),2)
      shrink.ccs[loc,]=c(a,b,i,igraph::transitivity(inet, type="global",isolates = "zero"),2)
      shrink.modularities[loc,]=c(a,b,i,cluster_walktrap(inet) %>% modularity(),2)
    }
  }
}

shrink.den.dif=(shrink.densities[1:1500,]-shrink.densities[1501:3000,])
shrink.den.dif[,c(1,2,3,5)]=shrink.densities[1:1500,c(1,2,3,5)]
shrink.cc.dif=(shrink.ccs[1:1500,]-shrink.ccs[1501:3000,])
shrink.cc.dif[,c(1,2,3,5)]=shrink.ccs[1:1500,c(1,2,3,5)]
shrink.mod.dif=(shrink.modularities[1:1500,]-shrink.modularities[1501:3000,])
shrink.mod.dif[,c(1,2,3,5)]=shrink.modularities[1:1500,c(1,2,3,5)]

# growth
growth.densities=data.frame(matrix(ncol = 5, nrow = 5*3*nets))
growth.ccs=data.frame(matrix(ncol = 5, nrow = 5*3*nets))
growth.modularities=data.frame(matrix(ncol = 5, nrow = 5*3*nets))
colnames(growth.densities)=c("Pn","Pr","i","value","type")
colnames(growth.ccs)=c("Pn","Pr","i","value","type")
colnames(growth.modularities)=c("Pn","Pr","i","value","type")
loc=0
for (a in 1:5){
  for (b in 1:3){
    for (i in 1:nets){
      filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/GrowthSim/Matrix",(a-1)*3+b,"A",
                      i,".csv")
      #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
      inet=read_graph(filename,format = "pajek")
      inet=simplify(inet)
      loc=loc+1
      growth.densities[loc,]=c(a,b,i,graph.density(inet),1)
      growth.ccs[loc,]=c(a,b,i,igraph::transitivity(inet, type="global",isolates = "zero"),1)
      growth.modularities[loc,]=c(a,b,i,cluster_walktrap(inet) %>% modularity(),1)
    }
  }
}

for (a in 1:5){
  for (b in 1:3){
    for (i in 1:nets){
      filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/GrowthSimStable/Matrix",(a-1)*3+b,"A",
                      i,".csv")
      #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
      inet=read_graph(filename,format = "pajek")
      inet=simplify(inet)
      loc=loc+1
      growth.densities[loc,]=c(a,b,i,graph.density(inet),2)
      growth.ccs[loc,]=c(a,b,i,igraph::transitivity(inet, type="global",isolates = "zero"),2)
      growth.modularities[loc,]=c(a,b,i,cluster_walktrap(inet) %>% modularity(),2)
    }
  }
}
system("say 'done'")
growth.den.dif=(growth.densities[1:1500,]-growth.densities[1501:3000,])
growth.den.dif[,c(1,2,3,5)]=growth.densities[1:1500,c(1,2,3,5)]
growth.cc.dif=(growth.ccs[1:1500,]-growth.ccs[1501:3000,])
growth.cc.dif[,c(1,2,3,5)]=growth.ccs[1:1500,c(1,2,3,5)]
growth.mod.dif=(growth.modularities[1:1500,]-growth.modularities[1501:3000,])
growth.mod.dif[,c(1,2,3,5)]=growth.modularities[1:1500,c(1,2,3,5)]

funmin <- function(x) {
  r <- mean(x)-std.error(x)
  names(r) <- "ymin"
  r
}
funmax <- function(x) {
  r <- mean(x)+std.error(x)
  names(r) <- "ymax"
  r
}

sccg=ggplot(shrink.mod.dif,aes(x=Pn+Pr/5-0.4,y=value,colour=as.factor(Pr))) +
  stat_summary(size=1.5,fun.y = mean, fun.ymin = funmin, fun.ymax = funmax) +
  theme_bw(base_size=10) + 
  theme(axis.text = element_text(size=13),axis.title = element_text(face="bold",size=15),
        legend.text= element_text(size=13), legend.title= element_text(size=13),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + 
  labs(y="Stable - Shrinking",x="Pn") +
  scale_x_continuous(breaks=1:5,labels=seq(0.5,0.9,0.1)) +
  scale_colour_discrete(name="Pr",labels=c(0.01, 0.05, 0.1))
  g+coord_cartesian(ylim=c(-0.9,0.05))


gmod=ggplot(growth.mod.dif,aes(x=Pn+Pr/5-0.4,y=value,colour=as.factor(Pr))) +
  stat_summary(size=1.5,fun.y = mean, fun.ymin = funmin, fun.ymax = funmax) +
  theme_bw(base_size=10) + 
  theme(axis.text = element_text(size=13),axis.title = element_text(face="bold",size=15),
        legend.text= element_text(size=13), legend.title= element_text(size=13),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + 
  labs(y="Stable - Growing",x="Pn") +
  scale_x_continuous(breaks=1:5,labels=seq(0.5,0.9,0.1)) +
  scale_colour_discrete(name  ="Pr",labels=c(0.01, 0.05, 0.1))

pdf("density_diff.pdf",width = 16,height = 7)
grid.arrange(sden,gden,ncol=2)
dev.off()

pdf("cc_diff.pdf",width = 16,height = 7)
grid.arrange(scc,gcc,ncol=2)
dev.off()

pdf("modularity_diff.pdf",width = 16,height = 7)
grid.arrange(smod,gmod,ncol=2)
dev.off()

# stable for comparison to shrinking
stable.densities=matrix(data = 0,nrow=5,ncol=6) #making room for mean+SD for 15 param sets
stable.ccs=matrix(data = 0,nrow=5,ncol=6)
stable.modularities=matrix(data = 0,nrow=5,ncol=6)
for (a in 1:5){
  for (b in 1:3){
    stable.dens=matrix(data = 0,nrow=nets,ncol=1)
    stable.clus=matrix(data = 0,nrow=nets,ncol=1)
    stable.mod=matrix(data = 0,nrow=nets,ncol=1) # stores modularity values
    for (i in 1:nets){
      filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/ShrinkSimStable/Matrix",(a-1)*3+b,"A",
                      i,".csv")
      #filename=paste0("/Volumes/HDD/Documents/Social Inheritance/Data/conflizard/Matrix",i,".csv")
      inet=read_graph(filename,format = "pajek")
      inet=simplify(inet)
      stable.dens[i,1]=graph.density(inet)
      stable.clus[i,1]=igraph::transitivity(inet, type="global",isolates = "zero")
      stable.mod[i,1] = cluster_walktrap(inet) %>% modularity()
    }
    stable.densities[a,b*2-1]=mean(stable.dens)
    stable.densities[a,b*2]=sd(stable.dens)
    stable.ccs[a,b*2-1]=mean(stable.clus)
    stable.ccs[a,b*2]=sd(stable.clus)
    stable.modularities[a,b*2-1]=mean(stable.mod)
    stable.modularities[a,b*2]=sd(stable.mod)
  }
}

res.den=round((stable.densities-shrink.densities)/stable.densities,3)[,c(1,3,5)]
res.cc=round((stable.ccs-shrink.ccs)/stable.ccs,3)[,c(1,3,5)]
res.mod=round((stable.modularities-shrink.modularities)/stable.modularities,3)[,c(1,3,5)]

library("reshape")
df.den<-melt(res.den)
gd=
  ggplot(df.den, aes(x=Var1, y=value)) + geom_point(aes(colour=as.factor(Var2)), size=6,position=position_dodge(width=0.4)) 
  + theme_bw(base_size=10) + 
  theme(panel.border=element_blank(),axis.text = element_text(size=8),axis.title = element_text(face="italic",size=15),axis.ticks = element_blank()) + 
  labs(y="Pr",x="Pn") +
  scale_x_continuous(expand = c(0, 0),breaks=1:10,labels=seq(0,.9,.1)) + 
  scale_y_continuous(expand = c(0, 0),breaks=1:10,labels=prseq) + coord_equal()
