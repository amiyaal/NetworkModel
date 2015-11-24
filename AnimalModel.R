library(MCMCglmm)
library(plyr)

pedigree = data.frame(animal=NA,sire=NA,dam=NA)
data.all=data.frame(degree=NA,cc=NA,betweenness=NA,animal=NA,network=NA,mother=NA)
for (i in 1:10) {
  pedigree.temp = data.frame(animal=0:99,sire=NA,dam=NA)
  for (j in 1:100) {
    pedigree.temp$animal[j]=paste0(letters[i],pedigree.temp$animal[j])
  }
  pedigree=rbind(pedigree,pedigree.temp)
}

for (i in 1:10){
  filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/AnimalModel1/Matrix",i,".csv")
  inet=read_graph(filename,format = "pajek")
  inet=simplify(inet)
  V(inet)$name=V(inet)$id
  parentsfile=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/AnimalModel1/Parents",i,".csv")
  parents <- read.csv(parentsfile, header=FALSE, sep=";", stringsAsFactors=FALSE)
  parents=data.frame(apply(parents,c(1,2),add.letter))
  colnames(parents)=c("animal","dam","sire")
  pedigree=rbind(pedigree,parents)
  data.deg=as.data.frame(degree(inet))
  data.deg$cc=transitivity(inet,type="local",isolates = "zero")
  data.deg$bet=betweenness(inet,directed = F)
  data.deg$animal=add.letter(V(inet)$name)
  data.deg$network=i
  data.deg$mother=NA
  for (j in 1:100){
    data.deg$mother[j]=parents$dam[parents$animal==data.deg$animal[j]]    
  }
  colnames(data.deg)=c("degree","cc","betweenness","animal","network","mother")
  data.all=rbind(data.all,data.deg)
}
pedigree=pedigree[-1,]
data.all=data.all[-1,]

data.all.08=data.all
pedigree.08=pedigree
data.all.rnd=data.all
pedigree.rnd=pedigree



add.letter=function(x){paste0(letters[i],x)}

prior <- list(R = list(V = 1, nu = 1), G = list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
                                                G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
                                                G3=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))
model1.deg <- MCMCglmm(degree ~ 1, random = ~animal + network + mother, family = "gaussian",
                       prior = prior, pedigree = pedigree, data = data.all, nitt = 60000,
                       burnin = 10000, thin = 10)

model.rnd.deg2=model1.deg
model.08.deg=model1.deg

herit.deg.08 <- model1.deg$VCV[, "animal"]/(model1.deg$VCV[, "animal"] + model1.deg$VCV[, "units"] + 
                                         model1.deg$VCV[, "network"]+ model1.deg$VCV[, "mother"])
mean(herit.deg.08)
HPDinterval(heritEC2)
plot(heritEC2)

i=1
filename=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/AnimalModelRandom/Matrix",i,".csv")
inet=read_graph(filename,format = "pajek")
inet=simplify(inet)
V(inet)$name=V(inet)$id
parentsfile=paste0("~/Dropbox/Erol/Cpp/Snap/examples/testgraph/AnimalModelRandom/Parents",i,".csv")
parents <- read.csv(parentsfile, header=FALSE, sep=";", stringsAsFactors=FALSE)
colnames(parents)=c("animal","dam","sire")
pedigree = data.frame(animal=0:99,sire=NA,dam=NA)
pedigree=rbind(pedigree,parents)
data.deg=as.data.frame(degree(inet))
data.deg$cc=transitivity(inet,type="local",isolates = "zero")
transitivity(inet)
data.deg$animal=V(inet)$name
colnames(data.deg)=c("degree","cc","animal")
plot(inet,vertex.label=NA, vertex.size=5, edge.color="black")

prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))
prior <- list(R = list(V = 1, nu = 1), G = list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))
model1.deg <- MCMCglmm(degree ~ 1, random = ~animal, family = "gaussian",
                  prior = prior, pedigree = pedigree, data = data.deg, nitt = 100000,
                  burnin = 10000, thin = 10)
model1.cc <- MCMCglmm(cc ~ 1, random = ~animal, family = "gaussian",
                       prior = prior, pedigree = pedigree, data = data.deg, nitt = 100000,
                       burnin = 10000, thin = 10)

system("say 'done'")

plot(model1.deg$Sol)
plot(model1.deg$VCV)
effectiveSize(model1.deg$Sol)
effectiveSize(model1.deg$VCV)
autocorr.diag(model1.deg$VCV)
heidel.diag(model1.deg$VCV)
summary(model1.deg)

herit <- model1.deg$VCV[, "animal"]/(model1.deg$VCV[, "animal"] + model1.deg$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) 
plot(herit)

#==== CC
plot(model1.cc$Sol)
plot(model1.cc$VCV)
effectiveSize(model1.cc$Sol)
effectiveSize(model1.cc$VCV)
autocorr.diag(model1.cc$VCV)
heidel.diag(model1.cc$VCV)
summary(model1.cc)

herit <- model1.cc$VCV[, "animal"]/(model1.cc$VCV[, "animal"] + model1.cc$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) 
plot(herit)
