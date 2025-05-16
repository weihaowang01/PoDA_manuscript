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
result_all=NULL
source("~/helpler_function.R")
source("~/PoDA.R")

for(mumu in c(2)){
  for(gammain in c(0.05,0.1,0.2)){
    for(nnum in c(50,100,200)){
      q=1
      num=100
      beta.tpr.test<-matrix(0,nrow=num,ncol=5)
      beta.fdr.test<-matrix(0,nrow=num,ncol=5)
      time=matrix(0,nrow=num,ncol=5)
      
      for(q in 1:num){
        formula <- 'u' # for A1-A3 & A5-A7 
        # formula <- 'u+z1+z2' # for A4
        # formula <- 'u+z1' # for A8-A10

        source("~/simulation/data_simulator_scene1.R")
        # source("~/simulation/data_simulator_scene2.R")
        ####################### sample filtration ########################
        sample_fil = which(colSums(Y)>=1000)
        Y=Y[,sample_fil]
        confoun=confoun[sample_fil,]
        confoun=NULL # if there are no confounders
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
        
        ###################################
        ############### PoDA ##############
        ###################################
        rej=PODA(Y = Y, u = u, confoun = confoun ,lmdpi=0.5,rep_time=50,fdrnomial=0.05 ,lib_cut=0, pre_cut=0)
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
            print(sum(H[rej] == 0) / length(rej))
            print(sum(H[rej] == 1) / sum(H))
            beta.tpr.test[q,1] <- sum(H[rej] == 1) / sum(H)
          }
        }
        
        ###################################
        ############## LinDA ##############
        ###################################
        res <- LinDA::linda(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = sample_cut,prev.cut = 0.2)
        rej= as.numeric(gsub('[taxon]',' ',rownames(res$otu.tab.use)[which(res$output[[1]]$reject)]))
        print(sum(H[rej] == 0) / length(rej))
        print(sum(H[rej] == 1) / sum(H))
        if(length(rej) == 0) {
          beta.fdr.test[q,2]<- 0
          beta.tpr.test[q,2] <- 0
        } else {
          beta.fdr.test[q,2]<- sum(H[rej] == 0) / length(rej)
          beta.tpr.test[q,2] <- sum(H[rej] == 1) / sum(H)
        }
        
        
        ###################################
        ############# ANCOMBC2 ############
        ###################################
        start_time <- as.numeric(Sys.time())
        rej.ancombc2 <- ancombc2.fun.cont(m,Y, Z, formula, alpha=0.05,p_adj_method="BH")
        end_time <- as.numeric(Sys.time())
        time[q,3] <- end_time - start_time
        sum(H[rej.ancombc2] == 1) / sum(H)
        sum(H[rej.ancombc2] == 0) / length(rej.ancombc2)
        #
        rej.ancombc2=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc2]))
        if(length(rej) == 0) {
          beta.fdr.test[q,3]<- 0
          beta.tpr.test[q,3] <- 0
        } else {
          beta.fdr.test[q,3]<- sum(H[rej.ancombc2] == 0) / length(rej.ancombc2)
          beta.tpr.test[q,3] <- sum(H[rej.ancombc2] == 1) / sum(H)
        }

        
        ###################################
        ############ ANCOMBC #############
        ###################################
        start_time <- as.numeric(Sys.time())
        rej.ancombc <- ancombc1.fun(m,Y, Z, formula, alpha=0.05,p_adj_method = "BH")
        end_time <- as.numeric(Sys.time())
        time[q,4] <- end_time - start_time
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


        ###################################
        ############## locom ##############
        ###################################
        start_time <- as.numeric(Sys.time())
        tt <- try(locom(otu.table = t(Y), Y = u, C=NULL, fdr.nominal = 0.05, seed = 1, n.cores = 10,filter.thresh=0)
                  , silent = TRUE)
        end_time <- as.numeric(Sys.time())
        time[q,1] <- end_time - start_time
        if(class(tt)=="try-error"){
          beta.fdr.test[q,5]<- 0
          beta.tpr.test[q,5] <- 0
        }else{
          rej.locom = as.numeric(gsub('[taxon]',' ',colnames(tt$q.otu)[which(tt$q.otu<0.05)]))
          
          if(length(rej.locom) == 0) {
            beta.fdr.test[q,5]<- 0
            beta.tpr.test[q,5] <- 0
          } else {
            beta.fdr.test[q,5]<- sum(H[rej.locom] == 0) / length(rej.locom)
            beta.tpr.test[q,5] <- sum(H[rej.locom] == 1) / sum(H)
          }
        }
        
      }

      result_single<-c(colMeans(beta.tpr.test,na.rm = TRUE),colMeans(beta.fdr.test,na.rm = TRUE),colMeans(time))
      result_all=rbind(result_all,result_single)
    }
  }
}
print(result_all)
# saveRDS(result_all,"~/simulation/dependency.rds")