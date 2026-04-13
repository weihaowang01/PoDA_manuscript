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
source("~/manuscript/helper_function.R")
source("~/manuscript/PODA.R")
source("~/manuscript/imputation.R")


#################################################################
######################### examination ###########################
#################################################################
mumu=2
gammain=0.1
for(rt in c(1,2,6,10,30,50,80,100,150,200,250)){ # rep_time
  for(nnum in c(50,100,200)){
      q=1
      num=100
      beta.tpr.test<-matrix(0,nrow=num,ncol=5)
      beta.fdr.test<-matrix(0,nrow=num,ncol=5)
      time=matrix(0,nrow=num,ncol=5)
      
      for(q in 1:num){
        
        source("~/manuscript/simulation/data_simulator_scene1.R")
        
        #################################################################
        ########################## preparation ##########################
        #################################################################
        # Please comment or retain the following code modules based on the corresponding simulation setting,
        # with the default set to setting A1.
        
        # ############### setting A1################
        # # To generate data under setting  A1, please comment out this section of the code.
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
        
        
        # ######################## setting A4 ##########################
        # # To generate data under setting A4, please comment out this section of the code.
        # formula <- 'u+z1+z2' 
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
        
        


        ###################################
        ############### PoDA ##############
        ###################################
        rej=PODA(Y = Y, u = u, confoun = confoun ,rep_time=rt,fdrnomial=0.05 ,lib_cut=0, pre_cut=0)
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
      result_all=rbind(result_all,result_single)
  }
}

print(result_all)

# if the current execution is based on setting A1.
tmp_bin=t(data.frame(result_all[,1]))
colnames(tmp_bin)=rep(c("1","2","4","6","10","30",
                   "50","80","100","150","200","250"),each=3)
saveRDS(tmp_bin,"~/manuscript/simulation/s-rep_time_examination/tmp_bin.rds")

# if the current execution is based on setting A4.
# tmp_con=t(data.frame(result_all[,1]))
# colnames(tmp_con)=rep(c("1","2","4","6","10","30",
#                         "50","80","100","150","200","250"),each=3)
# saveRDS(tmp_con,"~/manuscript/simulation/tmp_con.rds")

# #################################################################
# ############################ plot ###############################
# #################################################################
# 
# library(ggplot2)
# library(gridExtra)
# library(grid)
# 
# tmp_bin = readRDS("~/manuscript/simulation/tmp_bin.rds")
# tmp_con = readRDS("~/manuscript/simulation/tmp_con.rds")
# # 数据准备
# data <- expand.grid(
#   SampleSize = c("n = 50", "n = 100", "n = 200"),
#   SignalStrength = seq(1, 11, 1),
#   SignalType = c("A1", "A4")
# )
# 
# 
# data$EmpiricalFalseDiscoveryRate <- c( tmp_bin,tmp_con)*100
# 
# 
# data_n50 <- data %>% filter(SampleSize == "n = 50")
# data_n50_s1 <- data_n50 %>% filter(SignalType == "A1")
# data_n50_s7 <- data_n50 %>% filter(SignalType == "A4")
# 
# data_n100 <- data %>% filter(SampleSize == "n = 100")
# data_n100_s1 <- data_n100 %>% filter(SignalType == "A1")
# data_n100_s7 <- data_n100 %>% filter(SignalType == "A4")
# 
# 
# data_n200 <- data %>% filter(SampleSize == "n = 200")
# data_n200_s1 <- data_n200 %>% filter(SignalType == "A1")
# data_n200_s7 <- data_n200 %>% filter(SignalType == "A4")
# 
# 
# p_n50_s1 <- ggplot(data_n50_s1, aes(x = SignalStrength, 
#                                     y = EmpiricalFalseDiscoveryRate,
#                                     group = SignalType)) +
#   geom_line(size = 0.8, color = "black") +   
#   geom_point(aes(shape = SignalStrength == 6), size = 4, color = "black") + 
#   scale_shape_manual(values = c(1, 16)) +  
#   facet_grid(~ SampleSize) + 
#   scale_x_continuous(
#     breaks = seq(1, 11, 1),
#     labels = c("1", "2", "6", "10", "30", "50", "80", "100", "150", "200", "250")
#   ) + 
#   labs(
#     x = "Signal Strength",
#     y = "Empirical False Discovery Rate"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black"), 
#     strip.text = element_text(size = 12), 
#     panel.border = element_rect(color = "black", fill = NA), 
#     axis.title.y = element_blank(),       # 去掉纵坐标标题
#     axis.text.y = element_text(size = 10, angle = 90), 
#     axis.text.x = element_text(size = 10), 
#     legend.position = "none"                
#   )
# 
# p_n100_s1 <-  ggplot(data_n100_s1, aes(x = SignalStrength, 
#                                        y = EmpiricalFalseDiscoveryRate,
#                                        group = SignalType)) +
#   geom_line(size = 0.8, color = "black") +    # 黑色线条
#   geom_point(aes(shape = SignalStrength == 6), size = 4, color = "black") +  
#   scale_shape_manual(values = c(1, 16)) +  
#   facet_grid(~ SampleSize) + 
#   scale_x_continuous(
#     breaks = seq(1, 11, 1),
#     labels = c("1", "2", "6", "10", "30", "50", "80", "100", "150", "200", "250")
#   ) +  
#   labs(
#     x = "Signal Strength",
#     y = "Empirical False Discovery Rate"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black"),
#     strip.text = element_text(size = 12),
#     panel.border = element_rect(color = "black", fill = NA), 
#     axis.title.y = element_blank(),      
#     axis.text.y = element_text(size = 10, angle = 90),  
#     axis.text.x = element_text(size = 10), 
#     legend.position = "none"              
#   )
# 
# 
# 
# p_n200_s1 <-  ggplot(data_n200_s1, aes(x = SignalStrength, 
#                                        y = EmpiricalFalseDiscoveryRate,
#                                        group = SignalType)) +
#   geom_line(size = 0.8, color = "black") +    
#   geom_point(aes(shape = SignalStrength == 6), size = 4, color = "black") +  
#   scale_shape_manual(values = c(1, 16)) + 
#   facet_grid(SignalType ~ SampleSize) +  
#   scale_x_continuous(
#     breaks = seq(1, 11, 1),
#     labels = c("1", "2", "6", "10", "30", "50", "80", "100", "150", "200", "250")
#   ) +  
#   labs(
#     x = "Signal Strength",
#     y = "Empirical False Discovery Rate"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black"), 
#     strip.text = element_text(size = 12),
#     panel.border = element_rect(color = "black", fill = NA), 
#     axis.title.y = element_blank(),       
#     axis.text.y = element_text(size = 10, angle = 90), 
#     axis.text.x = element_text(size = 10), 
#     legend.position = "none"               
#   )
# 
# 
# p_n50_s7 <-  ggplot(data_n50_s7, aes(x = SignalStrength, 
#                                      y = EmpiricalFalseDiscoveryRate,
#                                      group = SignalType)) +
#   geom_line(size = 0.8, color = "black") +    # 黑色线条
#   geom_point(aes(shape = SignalStrength == 6), size = 4, color = "black") + 
#   scale_shape_manual(values = c(1, 16)) +  
#   scale_x_continuous(
#     breaks = seq(1, 11, 1),
#     labels = c("1", "2", "6", "10", "30", "50", "80", "100", "150", "200", "250")
#   ) + 
#   labs(
#     x = "Signal Strength",
#     y = "Empirical False Discovery Rate"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black"),
#     strip.text = element_text(size = 12),
#     panel.border = element_rect(color = "black", fill = NA),
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(size = 10, angle = 90),
#     axis.text.x = element_text(size = 10),
#     legend.position = "none"
#   )
# 
# p_n100_s7 <- ggplot(data_n100_s7, aes(x = SignalStrength, 
#                                       y = EmpiricalFalseDiscoveryRate,
#                                       group = SignalType)) +
#   geom_line(size = 0.8, color = "black") +
#   geom_point(aes(shape = SignalStrength == 6), size = 4, color = "black") +
#   scale_shape_manual(values = c(1, 16)) +
#   scale_x_continuous(
#     breaks = seq(1, 11, 1),
#     labels = c("1", "2", "6", "10", "30", "50", "80", "100", "150", "200", "250")
#   ) +
#   scale_y_continuous(
#     limits = c(80, 89)
#   ) +
#   labs(
#     x = "Signal Strength",
#     y = "Empirical False Discovery Rate"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black"),
#     strip.text = element_text(size = 12),
#     panel.border = element_rect(color = "black", fill = NA),
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(size = 10, angle = 90),
#     axis.text.x = element_text(size = 10),
#     legend.position = "none"
#   )
# 
# 
# p_n200_s7 <- ggplot(data_n200_s7, aes(x = SignalStrength, 
#                                       y = EmpiricalFalseDiscoveryRate,
#                                       group = SignalType)) +
#   geom_line(size = 0.8, color = "black") +
#   geom_point(aes(shape = SignalStrength == 6), size = 4, color = "black") +
#   scale_shape_manual(values = c(1, 16)) +
#   facet_grid(SignalType ~ .) +
#   scale_x_continuous(
#     breaks = seq(1, 11, 1),
#     labels = c("1", "2", "6", "10", "30", "50", "80", "100", "150", "200", "250")
#   ) +
#   labs(
#     x = "Signal Strength",
#     y = "Empirical False Discovery Rate"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black"),
#     strip.text = element_text(size = 12),
#     panel.border = element_rect(color = "black", fill = NA),
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(size = 10, angle = 90),
#     axis.text.x = element_text(size = 10),
#     legend.position = "none"
#   )
# 
# common_theme1 <- theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   axis.text.x = element_blank(),
#   panel.grid.minor = element_blank(),
#   panel.grid.major = element_line(size = 0.25, color = "gray90")
# )
# 
# common_theme <- theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   panel.grid.minor = element_blank(),
#   panel.grid.major = element_line(size = 0.25, color = "gray90")
# )
# p_n50_s1 <- p_n50_s1 + common_theme1
# p_n100_s1 <- p_n100_s1 + common_theme1
# p_n200_s1 <- p_n200_s1 + common_theme1
# p_n50_s7 <- p_n50_s7 + common_theme
# p_n100_s7 <- p_n100_s7 + common_theme
# p_n200_s7 <- p_n200_s7 + common_theme
# 
# 
# 
# 
# final_plot <- grid.arrange(
#   arrangeGrob(
#     p_n50_s1, p_n100_s1, p_n200_s1, 
#     p_n50_s7, p_n100_s7, p_n200_s7,
#     ncol = 3
#   )
# )
# 
# 
# 
# final_plot <- annotate_figure(
#   final_plot,
#   bottom = text_grob(expression(T), size = 12),  
#   left = text_grob("Power (%)", size = 12, rot = 90)
# )
# 
# final_plot




