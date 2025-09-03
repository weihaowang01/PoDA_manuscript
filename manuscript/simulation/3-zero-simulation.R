library(LinDA)
library(mvtnorm)
library(foreach)
library(doParallel)
library(LOCOM)
library(microbiome)
library(ANCOMBC)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidytree)
library(zCompositions)
result_all=NULL
source("~/manuscript/helpler_function.R")
source("~/manuscript/PODA.R")
source("~/manuscript/imputation.R")
source("~/manuscript/linda_modify.R")
source("~/manuscript/poda_modify.R")
for(mumu in c(2)){
  for(gammain in c(0.05,0.1,0.2)){
    for(nnum in c(50,100,200)){
      q=1
      num=100
      beta.tpr.test<-matrix(0,nrow=num,ncol=5)
      beta.fdr.test<-matrix(0,nrow=num,ncol=5)

      for(q in 1:num){
        formula <- 'u' # for A1 & A2
        formula <- 'u+z1' # for A10

        source("~/manuscript/simulation/data_simulator_scene1.R")
        # source("~/manuscript/simulation/data_simulator_scene2.R")
        ####################### sample filtration ########################
        sample_fil = which(colSums(Y)>=1000)
        Y=Y[,sample_fil]
        confoun=confoun[sample_fil,]
        # confoun=NULL # if there are no confounders
        u=u[sample_fil]

        name=colnames(Z)
        Z=as.matrix(Z[sample_fil,])
        colnames(Z)=name
        n=ncol(Y)
        m=nrow(Y)

        ####################### taxa filtration ########################
        datasize <- ncol(Y)
        prevalence = apply(as(Y, "matrix"), 1, function(x) {
          return(sum(x > 0))
        })/(datasize)

        keepOTUs = which(prevalence> 0.2)
        Y = Y[keepOTUs,]
        print(dim(Y))
        n=ncol(Y)
        m=nrow(Y)

        ################## zero-replacement strategies ##################
        # # pseudo=count
        Y=Y+0.5 # 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1

        # # multRepl
        # tY=t(Y)
        # Y=t(multRepl(tY,label=0,z.delete = FALSE))
        # colnames(Y)=rownames(tY)

        # # multLN
        # tY=t(Y)
        # Y=t(multLN(tY,label=0,z.delete = FALSE))

        # # adaImpute
        # Y<-imputation(Y,Z,formulaa)

        ##################################
        ############## PoDA ##############
        ##################################
        rej=PODA1(Y = Y, u = u, confoun = confoun ,rep_time=50,fdrnomial=0.05 ,lib_cut=0, pre_cut=0)
        rej = as.numeric(gsub('[taxon]',' ',rownames(Y)[rej$rej]))
        if(length(rej) == 0){
          beta.fdr.test[q,1]<- 0
          beta.tpr.test[q,1] <- 0
        }else{
          if(length(rej) == 0) {
            beta.fdr.test[q,1]<- 0
            beta.tpr.test[q,1] <- 0
          } else {
            beta.fdr.test[q,1]<- sum(H[rej] == 0) / length(rej)
            beta.tpr.test[q,1] <- sum(H[rej] == 1) / sum(H)
          }
        }


        ###################################
        ############## LinDA ##############
        ###################################
        res <- lindaa(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = 0,prev.cut = 0)
        rej= as.numeric(gsub('[taxon]',' ',rownames(res$otu.tab.use)[which(res$output[[1]]$reject)]))
        if(length(rej) == 0) {
          beta.fdr.test[q,2]<- 0
          beta.tpr.test[q,2] <- 0
        } else {
          beta.fdr.test[q,2]<- sum(H[rej] == 0) / length(rej)
          beta.tpr.test[q,2] <- sum(H[rej] == 1) / sum(H)
        }



        ###################################
        ############ ANCOMBC #############
        ###################################
        # Since ANCOM-BC automatically handles zeros by adding 1 to the count matrix,
        # the microbial count matrix needs to be subtracted by 1 before input here.
        rej.ancombc <- ancombc1.fun(m,Y-1, Z, formula, alpha=0.05,p_adj_method = "BH")
        rej.ancombc=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc]))
        sum(H[rej.ancombc] == 1)
        length(rej.ancombc)
        sum(H[rej.ancombc] == 1) / sum(H)
        sum(H[rej.ancombc] == 0) / length(rej.ancombc)
        if(length(rej) == 0) {
          beta.fdr.test[q,4]<- 0
          beta.tpr.test[q,4] <- 0
        } else {
          beta.fdr.test[q,4]<- sum(H[rej.ancombc] == 0) / length(rej.ancombc)
          beta.tpr.test[q,4] <- sum(H[rej.ancombc] == 1) / sum(H)
        }


      }

      result_single<-c(colMeans(beta.tpr.test,na.rm = TRUE),colMeans(beta.fdr.test,na.rm = TRUE))
      result_all=rbind(result_all,result_single)
      print(result_all)
    }
  }
}
print(result_all)
# saveRDS(result_all,"~/manuscript/simulation/Results/zero/A1_pseudo/A1pesudo5.rds")
