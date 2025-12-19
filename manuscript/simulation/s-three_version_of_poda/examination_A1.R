
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


#################################################################
############################ PODA ###############################
#################################################################
for(mumu in c(2)){ # effect size
  for(gammain in c(0.05,0.1,0.2)){ # signal proportion
    for(nnum in c(50,100,200)){ # sample size
      q=1
      num=100
      beta.tpr.test<-matrix(0,nrow=num,ncol=3)
      beta.fdr.test<-matrix(0,nrow=num,ncol=3)
      time=matrix(0,nrow=num,ncol=5)
      
      for(q in 1:num){
        
        source("~/manuscript/simulation/data_simulator_scene1.R")

        formula <- 'u' 
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

      }
      
      result_single1<-c(colMeans(beta.tpr.test,na.rm = TRUE),colMeans(beta.fdr.test,na.rm = TRUE))
      result_single2<-c(apply(beta.tpr.test, 2, sd, na.rm = TRUE)/sqrt(num),apply(beta.fdr.test, 2, sd, na.rm = TRUE)/sqrt(num))
      result_single=c(result_single1,result_single2)
      result_all_median=rbind(result_all,result_single)
    }
  }
}


#################################################################
########################## PODA-wmedian #########################
#################################################################
source("~/Desktop/newPODA/manuscript/simulation/s-three_version_of_poda/poda-wmedian.R")
for(mumu in c(2)){ # effect size
  for(gammain in c(0.05,0.1,0.2)){ # signal proportion
    for(nnum in c(50,100,200)){ # sample size
      q=1
      num=100
      beta.tpr.test<-matrix(0,nrow=num,ncol=3)
      beta.fdr.test<-matrix(0,nrow=num,ncol=3)
      time=matrix(0,nrow=num,ncol=5)
      
      for(q in 1:num){
        
        source("~/manuscript/simulation/data_simulator_scene1.R")
        
        formula <- 'u' 
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
        
      }
      
      result_single1<-c(colMeans(beta.tpr.test,na.rm = TRUE),colMeans(beta.fdr.test,na.rm = TRUE))
      result_single2<-c(apply(beta.tpr.test, 2, sd, na.rm = TRUE)/sqrt(num),apply(beta.fdr.test, 2, sd, na.rm = TRUE)/sqrt(num))
      result_single=c(result_single1,result_single2)
      result_all_wmedian=rbind(result_all,result_single)
    }
  }
}



#################################################################
########################### PODA-mean ###########################
#################################################################
source("~/Desktop/newPODA/manuscript/simulation/s-three_version_of_poda/poda-wmean.R")
for(mumu in c(2)){ # effect size
  for(gammain in c(0.05,0.1,0.2)){ # signal proportion
    for(nnum in c(50,100,200)){ # sample size
      q=1
      num=100
      beta.tpr.test<-matrix(0,nrow=num,ncol=3)
      beta.fdr.test<-matrix(0,nrow=num,ncol=3)
      time=matrix(0,nrow=num,ncol=5)
      
      for(q in 1:num){
        
        source("~/manuscript/simulation/data_simulator_scene1.R")
        
        formula <- 'u' 
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
        
      }
      
      result_single1<-c(colMeans(beta.tpr.test,na.rm = TRUE),colMeans(beta.fdr.test,na.rm = TRUE))
      result_single2<-c(apply(beta.tpr.test, 2, sd, na.rm = TRUE)/sqrt(num),apply(beta.fdr.test, 2, sd, na.rm = TRUE)/sqrt(num))
      result_single=c(result_single1,result_single2)
      result_all_wmean=rbind(result_all,result_single)
    }
  }
}

abs0=matrix(0,nrow=9,ncol=6)
abs0[,c(1,4)]=result_all_wmedian[,c(1,4)]
abs0[,c(2,5)]=result_all_wmean[,c(1,4)]
abs0[,c(3,6)]=result_all_median[,c(1,4)]
saveRDS(abs0,"~/manuscript/simulation/s-three_version_of_poda/three_version.rds")

# #################################################################
# ############################ plot ###############################
# #################################################################
# 
# abs0=readRDS("~/manuscript/simulation/s-three_version_of_poda/three_version.rds")
# colind = c(2,1,3)
# colind = c(colind, colind + 3)
# 
# simlist = t(abs0[1:3, colind])
# simlist = simlist * 100
# rownames(simlist) = c(
#   "power.IPOD-mean","power.IPOD-wm","power.IPOD-median",
#   "fdr.IPOD-mean","fdr.IPOD-wm","fdr.IPOD-median"
# )
# 
# f = function(i){
#   data.frame(
#     class  = rep(c("Power(%)","Empirical FDR(%)"), each = 3),
#     method = rep(c('PoDA-wmean','PoDA-wmedian','PoDA'), 2),
#     value  = simlist[, i]
#   )
# }
# 
# power.fdr = foreach(i = 1:3, .combine = rbind) %do% f(i)
# 
# out = data.frame(
#   sig = rep(c(50,100,200), each = 6),
#   q   = rep(c(0.5,5), each = 3),
#   power.fdr
# )
# 
# colnames(out) = c("sig.prob","q","class","method","value")
# 
# out$sig.prob = factor(out$sig.prob, levels = c("50","100","200"))
# out$class    = factor(out$class, levels = c("Power(%)","Empirical FDR(%)"))
# out$method   = factor(out$method, levels = c('PoDA','PoDA-wmedian','PoDA-wmean'))
# out$q        = factor(out$q, levels = c("0.5","5"))
# 
# line = data.frame(
#   class = factor(c("Power(%)","Empirical FDR(%)"),
#                  levels = c("Power(%)","Empirical FDR(%)")),
#   y = c(NA, 5)
# )
# 
# p1 = ggplot(out[out$q=="5",], aes(x = sig.prob, y = value, fill = method)) +
#   theme_bw()
# 
# ptest1 = p1 +
#   geom_bar(stat="identity", position=position_dodge(0.6), width=0.5) +
#   scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284','#E64B3599')) +
#   facet_grid(vars(class), labeller = labeller(.cols = label_both)) +
#   ylim(c(0,100)) +
#   geom_hline(data=line, aes(yintercept=y), linetype="dashed") +
#   theme(
#     legend.position="top",
#     legend.title=element_blank(),
#     legend.text=element_text(size=6),
#     legend.key.size=unit(0.3,"cm"),
#     axis.title=element_text(size=8),
#     axis.text=element_text(size=6),
#     panel.spacing=unit(0.5,"lines")
#   ) +
#   guides(fill=guide_legend(ncol=6)) +
#   ggnewscale::new_scale_fill() +
#   geom_bar(
#     data=out[out$q=="0.5",],
#     aes(x=sig.prob, y=value, fill=method),
#     stat="identity",
#     position=position_dodge(0.6),
#     width=0.5,
#     show.legend=FALSE
#   ) +
#   scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284','#E64B3599')) +
#   facet_grid(vars(class), labeller = labeller(.cols = label_both)) +
#   ylim(c(0,100)) +
#   labs(y=" ", x=expression(n))
# 
# ptest1 <- ptest1 + theme(
#   axis.title.x = element_blank(),
#   legend.text = element_text(size = 10),
#   legend.key.size = unit(0.5, "cm"),
#   axis.title.y = element_text(size = 10),
#   plot.margin = margin(1,1,1,1),
#   strip.background = element_blank(),
#   strip.text = element_blank()
# )
# 
# p.comb <- ggarrange(
#   ptest1, ptest2, ptest3,
#   ncol = 3,
#   common.legend = TRUE,
#   legend = "top"
# )
# 
# p.comb <- annotate_figure(
#   p.comb,
#   bottom = text_grob(expression(n), size = 12)
# )
# 
# p.comb
