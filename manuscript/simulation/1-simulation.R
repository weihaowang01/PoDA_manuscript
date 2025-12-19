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
  result_all=NULL
  source("~/manuscript/helpler_function.R")
  source("~/manuscript/PODA.R")
  source("~/manuscript/imputation.R")
  for(mumu in c(2)){ # effect size
    for(gammain in c(0.05,0.1,0.2)){ # signal proportion
      for(nnum in c(50,100,200)){ # sample size
        q=1
        num=100
        beta.tpr.test<-matrix(0,nrow=num,ncol=5)
        beta.fdr.test<-matrix(0,nrow=num,ncol=5)
        time=matrix(0,nrow=num,ncol=5)

        for(q in 1:num){
          
          source("~/manuscript/simulation/data_simulator_scene1.R")
          # source("~/manuscript/simulation/data_simulator_scene2.R")
          
          #################################################################
          ########################## preparation ##########################
          #################################################################
          # Please comment or retain the following code modules based on the corresponding simulation setting,
          # with the default set to setting A1.
          
          # ############### settings A1-A3 & A5-A7 & A11-A12################
          # # To generate data under settings  A1-A3 & A5-A7 & A11-A12, please comment out this section of the code.
          formula <- 'u' # for A1-A3 & A5-A7 & A11-A12
          # sample filtration 
          sample_fil = which(colSums(Y)>=1000)
          Y=Y[,sample_fil]
          confoun=NULL # if there are no confounders
          u=u[sample_fil]
          name=colnames(Z)
          Z=as.data.frame(Z[sample_fil,])
          Z[,1]=factor(Z[,1])
          colnames(Z)=name
          n=ncol(Y)
          m=nrow(Y)
          
          
          # ######################## setting A4 ##########################
          # # To generate data under setting A4, please comment out this section of the code.
          # formula <- 'u+z1+z2' =
          # # sample filtration 
          # sample_fil = which(colSums(Y)>=1000)
          # Y=Y[,sample_fil]
          # confoun=confoun[sample_fil,]
          # u=u[sample_fil]
          # name=colnames(Z)
          # Z=as.data.frame(Z[sample_fil,])
          # Z[,1]=factor(Z[,1])
          # colnames(Z)=name
          # n=ncol(Y)
          # m=nrow(Y)
          
          # ######################## settings A8-A10 ##########################
          # # To generate data under settings A8-A10, please comment out this section of the code.
          # formula <- 'u+z1' 
          # # sample filtration 
          # sample_fil = which(colSums(Y)>=1000)
          # Y=Y[,sample_fil]
          # confoun=confoun[sample_fil,]
          # u=u[sample_fil]
          # name=colnames(Z)
          # Z=as.data.frame(Z[sample_fil,])
          # Z[,1]=factor(Z[,1])
          # colnames(Z)=name
          # n=ncol(Y)
          # m=nrow(Y)
          
          
          # taxa filtration
          datasize <- ncol(Y)
          prevalence = apply(as(Y, "matrix"), 1, function(x) {
            return(sum(x > 0))
          })/(datasize)
          
          keepOTUs = which(prevalence> 0.01)
          Y = Y[keepOTUs,]
          print(dim(Y))
          n=ncol(Y)
          m=nrow(Y)
          
          
          #################################################################
          ###################### competing methods ########################
          #################################################################
          

          ###################################
          ############### PoDA ##############
          ###################################
          rej=PODA(Y = Y, u = u, confoun = confoun ,rep_time=50,fdrnomial=0.05 ,lib_cut=0, pre_cut=0)
          rej = as.numeric(gsub('[taxon]',' ',rownames(Y)[rej$rej]))
          rej
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
          res <- LinDA::linda(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = 0,prev.cut = 0)
          rej= as.numeric(gsub('[taxon]',' ',rownames(res$otu.tab.use)[which(res$output[[1]]$reject)]))
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
          ####### without filteration under (except for A7) #######
          rej.ancombc2 <- ancombc2.fun(m,Y, Z, formula, alpha=0.05,p_adj_method="BH")
          rej.ancombc2=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc2]))

          if(length(rej.ancombc2) == 0) {
            beta.fdr.test[q,3]<- 0
            beta.tpr.test[q,3] <- 0
          } else {
            beta.fdr.test[q,3]<- sum(H[rej.ancombc2] == 0) / length(rej.ancombc2)
            beta.tpr.test[q,3] <- sum(H[rej.ancombc2] == 1) / sum(H)
          }
          
          
          # ####### without filteration under A7 #######
          # rej.ancombc2 <- ancombc2.fun.con(m,Y, Z, formula, alpha=0.05,p_adj_method="BH")
          # rej.ancombc2=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc2]))
          # 
          # if(length(rej) == 0) {
          #   beta.fdr.test[q,3]<- 0
          #   beta.tpr.test[q,3] <- 0
          # } else {
          #   beta.fdr.test[q,3]<- sum(H[rej.ancombc2] == 0) / length(rej.ancombc2)
          #   beta.tpr.test[q,3] <- sum(H[rej.ancombc2] == 1) / sum(H)
          # }
          
          
          # ####### with filteration under (except for A7) #######
          # rej.ancombc2 <- ancombc2.fun.ss(m,Y, Z, formula, alpha=0.05,p_adj_method="BH")
          # rej.ancombc2=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc2]))
          # 
          # if(length(rej.ancombc2) == 0) {
          #   beta.fdr.test[q,3]<- 0
          #   beta.tpr.test[q,3] <- 0
          # } else {
          #   beta.fdr.test[q,3]<- sum(H[rej.ancombc2] == 0) / length(rej.ancombc2)
          #   beta.tpr.test[q,3] <- sum(H[rej.ancombc2] == 1) / sum(H)
          # }
          
          
          # ####### with filteration under A7 #######
          # rej.ancombc2 <- ancombc2.fun.conss(m,Y, Z, formula, alpha=0.05,p_adj_method="BH")
          # rej.ancombc2=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc2]))
          # 
          # if(length(rej) == 0) {
          #   beta.fdr.test[q,3]<- 0
          #   beta.tpr.test[q,3] <- 0
          # } else {
          #   beta.fdr.test[q,3]<- sum(H[rej.ancombc2] == 0) / length(rej.ancombc2)
          #   beta.tpr.test[q,3] <- sum(H[rej.ancombc2] == 1) / sum(H)
          # }


          # ###################################
          # ############ ANCOMBC #############
          # ###################################
          rej.ancombc <- ancombc1.fun(m,Y, Z, formula, alpha=0.05,p_adj_method = "BH")
          rej.ancombc=as.numeric(gsub('[taxon]',' ',rownames(Y)[rej.ancombc]))
          if(length(rej.ancombc) == 0) {
            beta.fdr.test[q,4]<- 0
            beta.tpr.test[q,4] <- 0
          } else {
            beta.fdr.test[q,4]<- sum(H[rej.ancombc] == 0) / length(rej.ancombc)
            beta.tpr.test[q,4] <- sum(H[rej.ancombc] == 1) / sum(H)
          }
          print(beta.fdr.test[q,4])
          print(beta.tpr.test[q,4])

          ###################################
          ############## LOCOM ##############
          ###################################
        tt <- try(locom(otu.table = t(Y), Y = u, C=confoun, fdr.nominal = 0.05, seed = 1, n.cores = 10,filter.thresh=0,
                        n.perm.max = 100000), silent = TRUE)
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
          print(beta.fdr.test[q,5])
          print(beta.tpr.test[q,5])
        print(q)
        }

        result_single1<-c(colMeans(beta.tpr.test,na.rm = TRUE),colMeans(beta.fdr.test,na.rm = TRUE))
        result_single2<-c(apply(beta.tpr.test, 2, sd, na.rm = TRUE)/sqrt(num),apply(beta.fdr.test, 2, sd, na.rm = TRUE)/sqrt(num))
        result_single=c(result_single1,result_single2)
        result_all=rbind(result_all,result_single)
      }
    }
  }
  print(result_all)
  # saveRDS(result_all,"~/manuscript/simulation/Results/main/A1binomialSE.rds")
