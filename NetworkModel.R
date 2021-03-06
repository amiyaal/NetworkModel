library("network")
library("sna")
library("igraph")
library("beepr")
library("plotrix")
library("Matrix")


N=62 # number of individuals
E=1 # environmental factor. Will make ties more likely if higher
b=0.333 # probability to not connect to mother when born
pn=0.7 # probability to connect to friend of mother
pr=0.01 # probability to connect to random individual
r=0.5
n=1000 #simulation iterations
born.to.central=0 # Do more central mothers have more newborns?
mr=25 # mutation rate is 1/mr
tthre=0.05 # trait similarity threshold. More likely to connect if difference between trait values are below this threshold
tbased = 0 # should connection be based on trait similarity
init.ties=0.1 # density of initial network

simulate = function(n, N, mr, b, pn, pr, tbased, init.ties, kill.method=1){
  # Initializing
  init=sample(0:1,N*N,replace=T,prob=c(1-init.ties,init.ties)) # start from random graph
  m=matrix(init, ncol=N,nrow=N)
  m=matrix(forceSymmetric(m),ncol=N,nrow=N)
  trait=runif(N) # vector of random initial trait values
  cc=rep(0,n) # will record clustering coefficient every generation
  density=rep(0,n) # will record network density
  degree.of.dead=rep(0,n) # will record the degree of the dead individual
  assort=rep(0,n) # will record trait assortativity
  rc.prob=E*pr
  fc.prob=E*pn
  tbased=E*pr
  for (i in 1:n){
    if (born.to.central == 0) mother=sample(1:N,1)
    else mother=sample(1:N,1,prob=rowSums(m))
    vecpn=m[,mother]*sample(0:1,N,prob=c(1-fc.prob,fc.prob),replace=T)
    vecpr=(1-m[,mother])*sample(0:1,N,prob=c(1-rc.prob,rc.prob),replace=T)
    vec=vecpn+vecpr
    vec[vec==2]=1
    vec=c(vec,1)
    #vec=rep(0,N+1)
    m=cbind(m,vec)
    m=rbind(m,vec)
    if (runif(1) > b){
      m[mother,N+1]=1
      m[N+1,mother]=1
    }
    trait=c(trait,trait[mother])
    if (sample(1:mr,1)==1){ # mutation happens
      #trait[N+1]=i+N
      trait[N+1]=rnorm(1,mean=trait[mother],sd=0.1) # newborn gets a trait similar to his mother
      if (trait[N+1]>1) trait[N+1]=trait[N+1]-1 # making a torus
      if (trait[N+1]<0) trait[N+1]=1+trait[N+1]
    } 
#    for (j in 1:N){ # connecting the newborn
 #     if (j==mother) next
      #       if(trait[mother]-trait[j]==tthre) { 
      #         rc.prob=sqrt(E*pr) #random connection increases when sharing traits
      #         fc.prob=sqrt(E*pn) #connections to mother's friends increase with sharing traits
      #         tbased=E*pn
      #       }
      #       else {
      
      #       }
#      if (tbased == pr){
 #       if (m[mother,j]==1) connect=sample(0:1,1,prob=c(1-fc.prob,fc.prob))
  #      else connect=sample(0:1,1,prob=c(1-rc.prob,rc.prob))
  #    }
  #    else connect=sample(0:1,1,prob=c(1-tbased,tbased))
  #    m[j,N+1]=connect
  #    m[N+1,j]=connect
  #    m[N+1,N+1]=1
  #  }
    #  j=sample(1:N,1) #modifying connections among existing individuals
    #  k=sample(1:N,1)
    # connect according to correlation in friends
    #   if (k!=j && m[j,k]==0) {
    #     cor=1-sum(abs(m[j,]-m[k,]))/N # the correlation in their friends
    #     #if (m[j,k]==0) connect=sample(0:1,1,prob=c(1-cor,cor))
    #     connect=sample(0:1,1,prob=c(1-cor,cor))
    #     #m[j,k]=connect
    #     #m[k,j]=connect
    #   }
    if (kill.method==1) dead=sample(1:(N+1),1) # choosing who to kill randomly
    else dead=sample(1:(N+1),1,prob=1/(sqrt(rowSums(m))+1)) # Higher degree -> less chance to die
    #degree.of.dead[i]=rowSums(m)[dead]
    m=m[-dead,-dead]
    trait=trait[-dead]
    #E=E-0.00005
    #newnet=network(m,directed=F)
    #if (n %% 1000 == 0) plot(newnet,displayisolates=T,edge.col="black",vertex.border="black")
    #cc[i]=gtrans(newnet,mode="graph")
    #density[i]=network.density(newnet)
    #traits[i]=trait[mother]
    igraph=graph.adjacency(m,mode="undirected")
    density[i]=graph.density(igraph)
    cc[i]=transitivity(igraph,type="global")
    assort[i]=assortativity(igraph,trait,directed=F)
    if (is.nan(assort[i])==T) assort[i]=0
  }
  res = list(m=m,trait=trait,degree.of.dead=degree.of.dead,cc=cc,density=density
             ,assort=assort)
  res
}

sim.res=simulate(n, N, mr, pn, pr, tbased, init.ties, 2)
beep()

sim.res=simulate(n, 47, mr, 0.7, 0, tbased, 0.2, 2)
sim.net=graph.adjacency(sim.res$m,mode="undirected")


b=0.01
# Values for b (add 0 and remove 1):
prseq=(c(10^-(seq(0,2,by=0.25)),0))
prseq=round(prseq[order(prseq)],digits=3)
res.arr.k1=array(list(),c(10,10,100))
TRY: array(list(NULL), c(10,100))
for (a in 1:10){
#  for (b in seq(from=0,to=0.9,by=0.1)){
    for (iter in 1:100){
      cat(sprintf("%d %d\n", a, iter))
      res = simulate(n, N, mr, (a/10)-0.1, b, tbased, init.ties,kill.method=2)
      res.arr.k1[[a]][[1]][[iter]]=res
    }
#  }
}

save(res.arr.k3,file="resarrk3.RData")

#====================================================================
# summarizing results
den=array(0,c(10,10))
cc=array(0,c(10,10))
assort=array(0,c(10,10))
for (i in 0:9){
  name=load(paste0("~/Dropbox/Erol/NetworkModel/resarr",i,".RData"))
  res=eval(parse(text=name))
  for (j in 1:10){
    den1=0
    cc1=0
    assort1=0
    for (iter in 1:100){
      den1=den1+mean(res[[j]][[1]][[iter]]["density"][[1]])
      cc1=cc1+mean(res[[j]][[1]][[iter]]["cc"][[1]])
      assort1=assort1+mean(res[[j]][[1]][[iter]]["assort"][[1]])
    }
      den[i+1,j]=den1/100
      cc[i+1,j]=cc1/100
      assort[i+1,j]=assort1/100
    }
}
rm(res)

# Degree Distribution
mean.deg=array(0,c(10,10))
sd.deg=array(0,c(10,10))
kurt=array(0,c(10,10))
skew=array(0,c(10,10))
for (i in 0:9){
  name=load(paste0("~/Dropbox/Erol/NetworkModel/resarr",i,".RData"))
  res=eval(parse(text=name))
  for (j in 1:10){
    skew1=0
    kurt1=0
    mean.deg1=0
    sd1=0
    for (iter in 1:100){
      net=res[[j]][[1]][[iter]]["newnet"][[1]]
      deg=sna::degree(net,gmode="graph")
      mean.deg1=mean.deg1+mean(deg)
      sd1=sd1+sd(deg)
      skew1=skew1+skewness(deg)
      kurt1=kurt1+kurtosis(deg)
    }
    mean.deg[i+1,j]=mean.deg1/100
    sd.deg[i+1,j]=sd1/100
    kurt[i+1,j]=kurt1/100
    skew[i+1,j]=skew1/100
  }
  rm(list=ls(pattern="res.arr"))
}

cov=sd.deg/mean.deg # Coefficient of variation
library("gplots")
colfunc <- colorRampPalette(c("black", "red"))
heatmap.2(den,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Density",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")
heatmap.2(cc,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="CC",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")
heatmap.2(assort,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Assortativity",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")

heatmap.2(mean.deg,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Degree Mean",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")
heatmap.2(sd.deg,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Degree SD",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")
heatmap.2(cov,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Coefficient of Variation",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")
heatmap.2(kurt,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Kurtosis",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")
heatmap.2(skew,Rowv=F,Colv=F,trace="none",col=colfunc(15),ylab="Pr",main="Skewness",
          labCol=seq(from=0,to=0.9,by=0.1),labRow=prseq,xlab="Pn")


filled.contour(sd.deg,color.palette=topo.colors,x=seq(from=0,to=0.9,by=0.1),y=prseq,
               ylab="Pn",xlab="Pr")
#col.l <- colorRampPalette(c('white', 'black'))
#levelplot(cc,col.regions=col.l)
# 3D
library(rgl)
persp3d(x=seq(from=0,to=0.9,by=0.1), y=prseq, den, col="skyblue")
persp3d(x=seq(from=0,to=0.9,by=0.1), y=prseq, cc, col="blue")
persp3d(x=seq(from=0,to=0.9,by=0.1), y=prseq, assort, col="blue")

par(xpd=F)
res1den=array(0,c(10,10,100))
for (a in 1:10){
  for (iter in 1:100){
    res1den[a,1,iter]=mean(res.arr.k3[[a]][[1]][[iter]]["density"][[1]])
  }
}
boxplot(t(res1den[1:10,1,1:100]),main="Density",xaxt="n",xlab="pn")
axis(1,at=1:10,labels=seq(from=0,to=0.9,by=0.1))

res0cc=array(0,c(10,10,100))
for (a in 1:10){
  for (iter in 1:100){
    res0cc[a,1,iter]=mean(res.arr0[[a]][[1]][[iter]]["cc"][[1]])
  }
}
boxplot(t(res0cc[1:10,1,1:100]),main="CC",xaxt="n",xlab="pn")
axis(1,at=1:10,labels=seq(from=0,to=0.9,by=0.1))

res0assort=array(0,c(10,10,100))
for (a in 1:10){
  for (iter in 1:100){
    res0assort[a,1,iter]=mean(res.arr0[[a]][[1]][[iter]]["assort"][[1]])
  }
}
boxplot(t(res0assort[1:10,1,1:100]),main="Assortativity",xaxt="n",xlab="pn")
axis(1,at=1:10,labels=seq(from=0,to=0.9,by=0.1))


# Autocorrelation
ac=acf(res.arr[[10]][[1]][[1]]["assort"][[1]],lag.max=100)

i=15
# Continuous trait coloring
newnet=network(res.arr.k1[[8]][[1]][[i]]["m"][[1]],directed=F)
#cols=color.scale(res.arr[[1]][[1]][[1]]["trait"][[1]],cs1=1)
t=res.arr.k1[[8]][[1]][[i]]["trait"][[1]]
cols = rainbow(length(t))[rank(t)]
plot(newnet,displayisolates=T,vertex.col=cols,edge.col="black",vertex.border="black",vertex.cex=3)
im=res.arr.k1[[8]][[1]][[i]]["m"][[1]]
newneti=graph.adjacency(im,mode="undirected")
assortativity(newneti,t,directed=F)

par(mfrow=c(3,4),mar=c(1,1,1,1))
for (i in 1:10){
  t=res.arr0[[i]][[1]][[1]]["trait"][[1]]
  cols = rainbow(length(t))[rank(t)]
  plot(res.arr0[[i]][[1]][[1]]["newnet"][[1]],displayisolates=T,vertex.col=cols,edge.col="black",vertex.border="black",vertex.cex=2)  
}

#====================================================================
par(xpd=T)
plot(res["cc"][[1]],type="l",col="blue", ylim=c(min(res["density"][[1]],res["cc"][[1]],res["assort"][[1]]),1),xlab="Time",ylab="")
lines(res["density"][[1]],col="red")
lines(res["assort"][[1]],col="green")
legend ("top",inset=c(0,-0.275),c("Density","Clustering coeficient","Assortativty"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("red","blue","green"))

mean(res["degree.of.dead"][[1]])
mean(res["density"][[1]])
mean(res["cc"][[1]])
mean(res["assort"][[1]])


# Initializing
init.empty=rep(0,N*N) # start from empty matrix
init=sample(0:1,N*N,replace=T,prob=c(0.95,0.05)) # start from random graph
init=init.empty
initm=matrix(init, ncol=N,nrow=N)
diag(initm)=1

initn=network(initm,directed=F)
plot(initn,displayisolates=T)
network.density(initn)
gtrans(initn,mode="graph")
cug.test(initn,gtrans,cmode="edges",mode="graph")

# random relatedness matrix
relat=runif(N*N,0,0.5)
# generating trait values
#trait=1:N # nominal trait
trait=runif(N,0.1) # vector of random initial trait values
m=initm
cc=rep(0,n) # will record clustering coefficient every generation
density=rep(0,n) # will record network density
degree.of.dead=rep(0,n) # will record the degree of the dead individual
traits=rep(0,n) 
assort=rep(0,n) # will record trait assortativity


# mean(degree.of.dead)
# mean(density)
# mean(cc)
# mean(assort)

trait=as.numeric(factor(trait),levels=1:length(unique(trait)))
cols=brewer.pal(length(unique(trait)),"Paired")
plot(newnet,displayisolates=T,vertex.col=cols[trait],edge.col="black",vertex.border="black")
par(xpd=T)
plot(cc,type="l",col="blue", ylim=c(min(density,cc,assort),1),xlab="Time",ylab="")
lines(density,col="red")
lines(assort,col="green")
legend ("top",inset=c(0,-0.275),c("Density","Clustering coeficient","Assortativty"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("red","blue","green"))
legend (2000,1.5,c("Density","Clustering coeficient","Assortativty"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("red","blue","green"))

plot(density,cc,main=("Born to central mothers"))
cor.test(density,cc)
density.random=density
cc.random=cc
density09=density
cc09=cc
density05=density
cc05=cc
density01=density
cc01=cc


newnet=network(m,directed=F)

network.density(newnet)
gtrans(newnet,mode="graph")
cug.test(newnet,gtrans,cmode="edges",mode="graph")

#curve(1/N + E*(pn*x*(1+r)/N + pr), from=1,to=10)