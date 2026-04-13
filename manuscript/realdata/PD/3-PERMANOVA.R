library(vegan)
library(compositions)
pre_cut=0.1
age_cut=70
age_cut2=100
##################################################################################################
######################################## Wallen 1 ################################################
##################################################################################################
##### Wallen1
Wallen1 = readRDS("~/manuscript/realdata/PD/Wallen1.rds")
otutable1=Wallen1$otu
meta1 = Wallen1$meta
meta1$Age=as.numeric(meta1$Age)
keepsample = which(meta1$Age >= age_cut & meta1$Age < age_cut2)
otutable1=otutable1[keepsample,]
meta1=meta1[keepsample,]



datasize <- nrow(otutable1)
prevalence = apply(as(otutable1, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  pre_cut)
otutable1 = otutable1[,keepOTUs]


###### Wallen2
Wallen2 = readRDS("~/manuscript/realdata/PD/Wallen2.rds")
otutable=Wallen2$otu
meta = Wallen2$meta
meta$Age=as.numeric(meta$Age)
keepsample = which(meta$Age >= age_cut & meta$Age < age_cut2)
otutable=otutable[keepsample,]
meta=meta[keepsample,]


datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  pre_cut)
otutable = otutable[,keepOTUs]





taxs=union(colnames(otutable1),colnames(otutable))
# #
otutable2=otutable1[which(meta1$PD=="PD"),]
meta2=meta1[which(meta1$PD=="PD"),]
otutable3=otutable1[which(meta1$PD=="HC"),]
meta3=meta1[which(meta1$PD=="HC"),]

otu=rbind(otutable2,otutable3)
meta0=rbind(meta2,meta3)
grp=c(rep(1,nrow(otutable2)),rep(0,nrow(otutable3)))

Y <- t(otu)
u=grp
z1=as.numeric(factor(meta0$sex))
z2=as.numeric(factor(meta0$national))
z3=as.numeric((meta0$Age))

Z=cbind(as.matrix(u),as.matrix(z1),as.matrix(z2),as.matrix(z3))
ZZ=data.frame(factor(u),factor(z1),factor(z2),z3)
confoun=as.matrix(cbind(as.matrix(z1),as.matrix(z2),as.matrix(z3)))
colnames(confoun)=c("z1","z2","z3")
colnames(Z)=c("u","z1","z2","z3")
colnames(ZZ)=c("u","z1","z2","z3")
formula <- 'u+z1+z2+z3'
formulaa <- '~u+z1+z2+z3'
n=ncol(Y)
m=nrow(Y)


set.seed(123)



clr_data <- clr(as.matrix(t(Y) + 0.5))
dist_aitchison <- dist(clr_data)



set.seed(123)

res1 <- adonis2(
  dist_aitchison ~ u + z1 + z2 + z3,
  data = data.frame(Z),
  permutations = 99999,
  by = "margin"  
)








##################################################################################################
######################################## Wallen 2 ################################################
##################################################################################################

otutable2=otutable[which(meta$PD=="PD"),]
meta2=meta[which(meta$PD=="PD"),]
otutable3=otutable[which(meta$PD=="HC"),]
meta3=meta[which(meta$PD=="HC"),]

otu=rbind(otutable2,otutable3)
meta0=rbind(meta2,meta3)
grp=c(rep(1,nrow(otutable2)),rep(0,nrow(otutable3)))

Y <- t(otu)
u=grp
z1=as.numeric(factor(meta0$sex))
z3=as.numeric((meta0$Age))


Z=cbind(as.matrix(u),as.matrix(z1),as.matrix(z3))
ZZ=data.frame(factor(u),factor(z1),z3)
confoun=as.matrix(cbind(as.matrix(z1),as.matrix(z3)))
colnames(confoun)=c("z1","z3")
colnames(Z)=c("u","z1","z3")
colnames(ZZ)=c("u","z1","z3")
formula <- 'u+z1+z3'
formulaa <- '~u+z1+z3'
n=ncol(Y)
m=nrow(Y)



set.seed(123)



clr_data <- clr(as.matrix(t(Y) + 0.5))
dist_aitchison <- dist(clr_data)


set.seed(123)

res2 <- adonis2(
  dist_aitchison ~ u + z1 + z3,
  data = data.frame(Z),
  permutations = 99999,
  by = "margin"  
)


print(res1)
print(res2)
