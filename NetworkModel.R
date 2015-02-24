library("network")
library("sna")
library("igraph")
library("beepr")
library("plotrix")

N=50 # number of individuals
E=1 # environmental factor. Will make ties more likely if higher
pn=0.7 # probability to connect to friend of mother
pr=0.01 # probability to connect to random individual
r=0.5
n=5000 #simulation iterations
born.to.central=0 # Do more central mothers have more newborns?
mr=25 # mutation rate is 1/mr
tthre=0.05 # trait similarity threshold. More likely to connect if difference between trait values are below this threshold
tbased = 0 # should connection be based on trait similarity

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
#### MAKE CONTINUOUS TRAIT
#trait=1:N # nominal trait
trait=runif(N,0.1) # vector of random initial trait values
m=initm
cc=rep(0,n) # will record clustering coefficient every generation
density=rep(0,n) # will record network density
degree.of.dead=rep(0,n) # will record the degree of the dead individual
traits=rep(0,n) 
assort=rep(0,n) # will record trait assortativity

simulate = function(n, N, trait, mr, pn, pr, tbased){
  # Initializing
  init.empty=rep(0,N*N) # start from empty matrix
  init=sample(0:1,N*N,replace=T,prob=c(0.95,0.05)) # start from random graph
  init=init.empty
  m=matrix(init, ncol=N,nrow=N)
  for (i in 1:n){
    if (born.to.central == 0) mother=sample(1:N,1)
    else mother=sample(1:N,1,prob=rowSums(m))
    vec=rep(0,N+1)
    m=cbind(m,vec)
    m=rbind(m,vec)
    m[mother,N+1]=1
    m[N+1,mother]=1
    trait=c(trait,trait[mother])
    if (sample(1:mr,1)==1){ # mutation happens
      #trait[N+1]=i+N
      trait[N+1]=rnorm(1,mean=trait[mother],sd=0.1) # newborn gets a trait similar to his mother
      if (trait[N+1]>1) trait[N+1]=trait[N+1]-1 # making a torus
      if (trait[N+1]<0) trait[N+1]=1+trait[N+1]
    } 
    for (j in 1:N){ # connecting the newborn
      if (j==mother) next
      if(trait[mother]-trait[j]==tthre) { 
        rc.prob=sqrt(E*pr) #random connection increases when sharing traits
        fc.prob=sqrt(E*pn) #connections to mother's friends increase with sharing traits
        tbased=E*pn
      }
      else {
        rc.prob=E*pr
        fc.prob=E*pn
        tbased=E*pr
      }
      if (tbased == 0){
        if (m[mother,j]==1) connect=sample(0:1,1,prob=c(1-fc.prob,fc.prob))
        else connect=sample(0:1,1,prob=c(1-rc.prob,rc.prob))
      }
      else connect=sample(0:1,1,prob=c(1-tbased,tbased))
      m[j,N+1]=connect
      m[N+1,j]=connect
      m[N+1,N+1]=1
    }
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
    dead=sample(1:(N+1),1) # choosing who to kill
    degree.of.dead[i]=rowSums(m)[dead]
    m=m[-dead,-dead]
    trait=trait[-dead]
    #E=E-0.00005
    newnet=network(m,directed=F)
    #if (n %% 1000 == 0) plot(newnet,displayisolates=T,edge.col="black",vertex.border="black")
    cc[i]=gtrans(newnet,mode="graph")
    density[i]=network.density(newnet)
    #traits[i]=trait[mother]
    igraph=graph.adjacency(m,mode="undirected")
    assort[i]=assortativity(igraph,trait,directed=F)
  }
  res = list(m=m,newnet=newnet,trait=trait,degree.of.dead=degree.of.dead,cc=cc,density=density
             ,assort=assort)
  res
}

res=simulate(n, N, trait, mr, pn, pr, tbased)
beep()

mean(res["degree.of.dead"][[1]])
mean(res["density"][[1]])
mean(res["cc"][[1]])
mean(res["assort"][[1]])

# Continuous trait coloring
cols=color.scale(res["trait"][[1]],cs1=1)
plot(res["newnet"][[1]],displayisolates=T,vertex.col=cols,edge.col="black",vertex.border="black")

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