library(LinDA)
library(mvtnorm)
library(foreach)
library(doParallel)
library(LOCOM)
library(microbiome)
library(corncob)
library(ANCOMBC)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidytree)
source("~/manuscript/helper_function.R")
source("~/manuscript/PODA.R")
source("~/manuscript/imputation.R")

##################################################################################################
##################################################################################################
##################################################################################################
pre_cut=0.1
age_cut=70
age_cut2=100


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


taxs=union(colnames(otutable),colnames(otutable1))


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



Z=cbind(as.matrix(u),as.matrix(z1))
ZZ=data.frame(factor(u),factor(z1))
confoun=as.matrix(z1)
colnames(confoun)=c("z1")
colnames(Z)=c("u","z1")
colnames(ZZ)=c("u","z1")
formula <- 'u+z1'
formulaa <- '~u+z1'
n=ncol(Y)
m=nrow(Y)

######## perform different methods ########
# linda
lib_cut=1
res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut,imputation = TRUE)
rej.ld <- which(res$output[[1]]$reject)

# poda
rej=PODA( Y = Y, u = u, confoun = confoun,rep_time=500,tau = 1,fdrnomial=0.05 ,lib_cut=0,
          pre_cut=0,nCore=9)
rej2=rej$rej

# locom
tt = LOCOM::locom(otu.table = t(Y), Y = data.frame(factor(u)), C =confoun, fdr.nominal = 0.05,seed=100, n.cores = 10,filter.thresh = 0)
rej.locom = match(colnames(tt$q.otu)[which(tt$q.otu<=0.05)],rownames(Y))

#ancombc
rej.ancombc <- ancombc1.fun(m,Y, ZZ, formula, alpha=0.05,p_adj_method = "BH")

#ancombc2
rej.ancombc2 <- ancombc2.fun(m,Y, ZZ, formula, alpha=0.05,p_adj_method="BH")

rej5=match(rownames(Y)[rej2],taxs)
rej5.ld=match(rownames(Y)[rej.ld],taxs)
rej5.locom=match(rownames(Y)[rej.locom],taxs)
rej5.ancombc2=match(rownames(Y)[rej.ancombc2],taxs)
rej5.ancombc=match(rownames(Y)[rej.ancombc],taxs)

resres2$r100=rej5
lapply(resres2, length)
