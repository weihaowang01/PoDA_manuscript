library(LinDA)
library(mvtnorm)
library(foreach)
library(doParallel)
library(dacomp)
library(LOCOM)
library(microbiome)
library(corncob)
library(ANCOMBC)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidytree)
source("~/simulation/helpler_function.R")
source("~/poda-function/R/poda.R")
##################################################################################################
######################################## Wallen 1 ################################################
##################################################################################################

pre_cut=0.1
age_cut=70 
age_cut2=100


###### Wallen2
CDIlist2 = readRDS("~/Downloads/Wallen2.rds")
otutable=CDIlist2$otu
meta = CDIlist2$meta
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
CDIlist = readRDS("~/Downloads/Wallen1.rds")
otutable1=CDIlist$otu
meta1 = CDIlist$meta
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

# ##Wallen2
# otutable2=otutable[which(meta$PD=="PD"),]
# meta2=meta[which(meta$PD=="PD"),]
# otutable3=otutable[which(meta$PD=="HC"),]   
# meta3=meta[which(meta$PD=="HC"),]

#Wallen1
otutable2=otutable1[which(meta1$PD=="PD"),]
meta2=meta1[which(meta1$PD=="PD"),]
otutable3=otutable1[which(meta1$PD=="HC"),]
meta3=meta1[which(meta1$PD=="HC"),]
print("Wallen1")

otu=rbind(otutable2,otutable3)
meta0=rbind(meta2,meta3)
grp=c(rep(1,nrow(otutable2)),rep(0,nrow(otutable3)))


num=100
Error1=matrix(0,ncol=5,nrow=ncol(otu))
for(i in 1:num){
  print(i)
  Y <- t(otu)
  u=grp
  set.seed(i)
  u=sample(grp)
  Z=as.matrix(u)
  ZZ=data.frame(factor(u))
  confun=NULL
  
  colnames(Z)=c("u")
  colnames(ZZ)=c("u")
  formula <- 'u'
  formulaa <- '~u'

  res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = 0)
  rej.ld <- which(res$output[[1]]$reject)
  if(length(rej.ld)!=0){
    Error1[rej.ld,1]=1+Error1[rej.ld,1]
  }
  
  rej=PODA( Y = Y, u = u, confun = confun,lmdpi=0.5,rep_time=500,tau = 1,fdrnomial=0.05 ,lib_cut=0,
            pre_cut=0,nCore=9)
  rej=rej$rej
  if(length(rej)!=0){
    Error1[rej,2]=1+Error1[rej,2]
  }
  
  tt = LOCOM::locom(otu.table = t(Y), Y = data.frame(factor(u)), C =NULL, fdr.nominal = 0.05,seed=100, n.cores = 10,filter.thresh = 0)
  rej.locom = match(colnames(tt$q.otu)[which(tt$q.otu<=0.05)],rownames(Y))
  print(rej.locom)
  if(length(rej.locom)!=0){
    Error1[rej.locom,3]=1+Error1[rej.locom,3]
  }
  
  
  rej.ancombc <- ancombc1.fun(m,Y, ZZ, formula, alpha=0.05,p_adj_method = "BH")
  if(length(rej.ancombc)!=0){
    Error1[rej.ancombc,4]=1+Error1[rej.ancombc,4]
  }
  
  rej.ancombc2 <- ancombc2.fun(m,Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
  if(length(rej.ancombc2)!=0){
    Error1[rej.ancombc2,5]=1+Error1[rej.ancombc2,5]
  }
  
}


# ##################################################################################################
# ######################################## Wallen 2 ################################################
# ##################################################################################################
# pre_cut=0.1
# age_cut=70 
# age_cut2=100
# 
# 

# ###### Wallen2
# CDIlist2 = readRDS("~/Downloads/Wallen2.rds")
# otutable=CDIlist2$otu
# meta = CDIlist2$meta
# meta$Age=as.numeric(meta$Age)
# keepsample = which(meta$Age >= age_cut & meta$Age < age_cut2)
# otutable=otutable[keepsample,]
# meta=meta[keepsample,]
# 
# 
# datasize <- nrow(otutable)
# prevalence = apply(as(otutable, "matrix"), 2, function(x) {
#   return(sum(x > 0))
# })/(datasize)
# 
# keepOTUs = which(prevalence>  pre_cut)
# otutable = otutable[,keepOTUs]
# 
# 
# ##### Wallen1
# CDIlist = readRDS("~/Downloads/Wallen1.rds")
# otutable1=CDIlist$otu
# meta1 = CDIlist$meta
# meta1$Age=as.numeric(meta1$Age)
# keepsample = which(meta1$Age >= age_cut & meta1$Age < age_cut2)
# otutable1=otutable1[keepsample,]
# meta1=meta1[keepsample,]
# 
# 
# 
# datasize <- nrow(otutable1)
# prevalence = apply(as(otutable1, "matrix"), 2, function(x) {
#   return(sum(x > 0))
# })/(datasize)
# 
# keepOTUs = which(prevalence>  pre_cut)
# otutable1 = otutable1[,keepOTUs]
# 
# taxs=union(colnames(otutable),colnames(otutable1))
# 
# ##Wallen2
# otutable2=otutable[which(meta$PD=="PD"),]
# meta2=meta[which(meta$PD=="PD"),]
# otutable3=otutable[which(meta$PD=="HC"),]   
# meta3=meta[which(meta$PD=="HC"),]
# 
# ##Wallen1
# # otutable2=otutable1[which(meta1$PD=="PD"),]
# # meta2=meta1[which(meta1$PD=="PD"),]
# # otutable3=otutable1[which(meta1$PD=="HC"),]   
# # meta3=meta1[which(meta1$PD=="HC"),]
# # print("Wallen1")
# 
# otu=rbind(otutable2,otutable3)
# meta0=rbind(meta2,meta3)
# grp=c(rep(1,nrow(otutable2)),rep(0,nrow(otutable3)))
# 
# 
# num=100
# Error2=matrix(0,ncol=5,nrow=ncol(otu))
# for(i in 1:num){
#   print(i)
#   Y <- t(otu)
#   u=grp
#   set.seed(i)
#   u=sample(grp)
#   Z=as.matrix(u)
#   ZZ=data.frame(factor(u))
#   confun=NULL
#   
#   colnames(Z)=c("u")
#   colnames(ZZ)=c("u")
#   formula <- 'u'
#   formulaa <- '~u'
#   
#   res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = 0)
#   rej.ld <- which(res$output[[1]]$reject)
#   if(length(rej.ld)!=0){
#     Error2[rej.ld,1]=1+Error2[rej.ld,1]
#   }
#   
#   rej=PODA( Y = Y, u = u, confun = confun,lmdpi=0.5,rep_time=500,tau = 1,fdrnomial=0.05 ,lib_cut=0,
#             pre_cut=0,nCore=9)
#   rej=rej$rej
#   if(length(rej)!=0){
#     Error2[rej,2]=1+Error2[rej,2]
#   }
#   
#   tt = LOCOM::locom(otu.table = t(Y), Y = data.frame(factor(u)), C =NULL, fdr.nominal = 0.05,seed=100, n.cores = 10,filter.thresh = 0)
#   rej.locom = match(colnames(tt$q.otu)[which(tt$q.otu<=0.05)],rownames(Y))
#   print(rej.locom)
#   if(length(rej.locom)!=0){
#     Error2[rej.locom,3]=1+Error2[rej.locom,3]
#   }
#   
#   
#   rej.ancombc <- ancombc1.fun(m,Y, ZZ, formula, alpha=0.05,p_adj_method = "BH")
#   if(length(rej.ancombc)!=0){
#     Error2[rej.ancombc,4]=1+Error2[rej.ancombc,4]
#   }
#   
#   rej.ancombc2 <- ancombc2.fun(m,Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
#   if(length(rej.ancombc2)!=0){
#     Error2[rej.ancombc2,5]=1+Error2[rej.ancombc2,5]
#   }
#   
# }
# 
