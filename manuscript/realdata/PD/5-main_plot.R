# source("~/manuscript/realdata/PD/3-wallen1_detection.R")
# source("~/manuscript/realdata/PD/3-wallen2_detection.R")
library(patchwork)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)
rm(count)
sdsd=setdiff(rej4,Reduce(union,list(rej4.ld,rej4.locom,rej4.ancombc,rej4.ancombc2)))
taxs[sdsd]
xixi1=taxs[sdsd]
sdsd2=setdiff(rej5,Reduce(union,list(rej5.ld,rej5.locom,rej5.ancombc,rej5.ancombc2)))
taxs[sdsd2]
xixi2=taxs[sdsd2]
xixi3=taxs[intersect(rej4,rej5)]

findings=matrix(0,nrow=200,ncol=5)
rej.list=list()
rej.list[[1]]=rej4.ld
rej.list[[2]]=rej4
rej.list[[3]]=rej4.locom
rej.list[[4]]=rej4.ancombc
rej.list[[5]]=rej4.ancombc2
for(i in 1:5){
  findings[rej.list[[i]],i]=1
}


findings = as.data.frame(findings)
colnames(findings)=c("LinDA","PoDA","LOCOM","ANCOM-BC","ANCOM-BC2")

findings=findings[which(rowSums(findings)!=0),]


bar_colors <- c("#3D5c6f", rep("blue", 10))

upset_plot1 <-ComplexUpset::upset(
  findings,
  intersect = c("LinDA","PoDA","LOCOM","ANCOM-BC","ANCOM-BC2"),
  name = NULL,
  queries=list(
    ComplexUpset::upset_query(
      intersect=c('PoDA'),
      color="#E47159",
      fill="#E47159",
      only_components=c('intersections_matrix', 'Intersection size')
    ),
    ComplexUpset::upset_query("LinDA", fill="#F6BCA9"),
    ComplexUpset::upset_query("LOCOM", fill="#F6BCA9"),
    ComplexUpset::upset_query("ANCOM-BC", fill="#F6BCA9"),
    ComplexUpset::upset_query("ANCOM-BC2", fill="#F6BCA9"),
    ComplexUpset::upset_query("PoDA", fill="#E47159")
  ),
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      counts = TRUE,
      aes(fill = as.factor(seq_along(after_stat(count))))
    ) + scale_fill_manual(values = bar_colors, guide = 'none')
  )
)


upset_plot1
upset_plot1_wrp <- wrap_elements(full = upset_plot1)


##################################################################################################
##################################################################################################
##################################################################################################




findings=matrix(0,nrow=300,ncol=5)
rej.list=list()
rej.list[[1]]=rej5.ld
rej.list[[2]]=rej5
rej.list[[3]]=rej5.locom
rej.list[[4]]=rej5.ancombc
rej.list[[5]]=rej5.ancombc2
for(i in 1:5){
  findings[rej.list[[i]],i]=1
}
findings = as.data.frame(findings)
colnames(findings)=c("LinDA","PoDA","LOCOM","ANCOM-BC","ANCOM-BC2")

findings=findings[which(rowSums(findings)!=0),]


bar_colors <- c("#3D5c6f", rep("blue", 10))

upset_plot2 <-ComplexUpset::upset(
  findings,
  intersect = c("LinDA","PoDA","LOCOM","ANCOM-BC","ANCOM-BC2"),
  name = NULL,
  queries=list(
    ComplexUpset::upset_query(
      intersect=c('PoDA'),
      color="#E47159",
      fill="#E47159",
      only_components=c('intersections_matrix', 'Intersection size')
    ),
    ComplexUpset::upset_query("LinDA", fill="#F6BCA9"),
    ComplexUpset::upset_query("LOCOM", fill="#F6BCA9"),
    ComplexUpset::upset_query("ANCOM-BC", fill="#F6BCA9"),
    ComplexUpset::upset_query("ANCOM-BC2", fill="#F6BCA9"),
    ComplexUpset::upset_query("PoDA", fill="#E47159")
  ),
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      counts = TRUE,
      aes(fill = as.factor(seq_along(after_stat(count))))
    ) + scale_fill_manual(values = bar_colors, guide = 'none')
  )
)


upset_plot2
upset_plot2_wrp <- wrap_elements(full = upset_plot2)


##################################################################################################
##################################################################################################
##################################################################################################

.ratio_calcu=function(vec1,vec2){
  inter=intersect(vec1,vec2)
  union_set=union(vec1,vec2)
  setdiff1=setdiff(vec1,vec2)
  setdiff2=setdiff(vec2,vec1)

  inter_ratio <- length(inter) / length(union_set)
  setdiff_rej4_ratio <- length(setdiff1) / length(union_set)
  setdiff_rej4_ld_ratio <- length(setdiff2) / length(union_set)
  return(c(setdiff_rej4_ratio,inter_ratio,setdiff_rej4_ld_ratio))
}


ra=.ratio_calcu(rej4,rej5)
ra.ld=.ratio_calcu(rej4.ld,rej5.ld)
ra.locom=.ratio_calcu(rej4.locom,rej5.locom)
ra.ancombc=.ratio_calcu(rej4.ancombc,rej5.ancombc)
ra.ancombc2=.ratio_calcu(rej4.ancombc2,rej5.ancombc2)

data <- data.frame(
  Method = c("PoDA", "LinDA", "LOCOM", "ANCOM-BC", "ANCOM-BC2"),
  Ratio1 = c(ra[1], ra.ld[1], ra.locom[1], ra.ancombc[1], ra.ancombc2[1]),
  Ratio2 = c(ra[2], ra.ld[2], ra.locom[2], ra.ancombc[2], ra.ancombc2[2]),
  Ratio3 = c(ra[3], ra.ld[3], ra.locom[3], ra.ancombc[3], ra.ancombc2[3])
)

data_long <- pivot_longer(
  data,
  cols = starts_with("Ratio"),
  names_to = "Ratio",
  values_to = "Value"
)


library(scales)
data_long$custom_label <- percent(data_long$Value, accuracy = 0.01)

data_long$Method <- factor(data_long$Method, levels = c("PoDA", "LinDA", "LOCOM", "ANCOM-BC", "ANCOM-BC2"))
p.venn=ggplot(data_long, aes(x = Method, y = Value, fill = Ratio)) +
  geom_bar(stat = "identity",width = 0.7) +
  theme_minimal() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("#F9AE78", "#3D5C6F", "#E47159"),
                    labels = c(expression(omega[1] * " \\ " * omega[2]),
                               expression(omega[1] * intersect(omega[2])),
                               expression(omega[2] * " \\ " * omega[1])), name = NULL ) +
  scale_y_continuous(breaks = c(0, 1))+
  theme(
    legend.position = "none",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(t = 0, b = 0, l = 0, r = 1),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 5,
                               size=8),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,
    plot.margin = margin(0, 0, 0, 0)
  ) +
  geom_text(aes(label = custom_label), position = position_stack(vjust = 0.5), size = 3, color = "white")

p.venn







##################################################################################################
#####################################  replicability #############################################
##################################################################################################
###### Wallen2
# CDIlist2 = readRDS("~/PD/Wallen2.rds")
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

#
##### Wallen1
# CDIlist = readRDS("~/PD/Wallen1.rds")
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

repro = match(intersect(colnames(otutable),colnames(otutable1)),taxs)
repro1=length(intersect(intersect(rej4,rej5),repro))/length(intersect(union(rej4,rej5),repro))
repro2=length(intersect(intersect(rej4.ld,rej5.ld),repro))/length(intersect(union(rej4.ld,rej5.ld),repro))
repro3=length(intersect(intersect(rej4.ancombc,rej5.ancombc),repro))/length(intersect(union(rej4.ancombc,rej5.ancombc),repro))
repro4=length(intersect(intersect(rej4.ancombc2,rej5.ancombc2),repro))/length(intersect(union(rej4.ancombc2,rej5.ancombc2),repro))
repro5=length(intersect(intersect(rej4.locom,rej5.locom),repro))/length(intersect(union(rej4.locom,rej5.locom),repro))




methods <- c("PoDA", "LinDA", "LOCOM", "ANCOM-BC", "ANCOM-BC2")
jaccard_index <- c(repro1, repro2, repro5, repro3, repro4)

data <- data.frame(Method = methods, Jaccard_Index = jaccard_index)


data$Method <- factor(data$Method, levels = c("PoDA", "LinDA", "LOCOM", "ANCOM-BC", "ANCOM-BC2"))


p.sig2=ggplot(data, aes(x = Method, y = Jaccard_Index)) +
  geom_bar(stat = "identity", aes(fill = Method), position=position_dodge(0.1),width = 0.5) +
  scale_fill_manual(values = c("PoDA" = "#E47159", "LinDA" = "#3D5c6f", "LOCOM" = "#3D5c6f", "ANCOM-BC" = "#3D5c6f", "ANCOM-BC2" = "#3D5c6f")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1,size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.title.y = element_text(size = 10)
  )+
  ylab("Jaccard index")+xlab("")

main1=upset_plot1_wrp+upset_plot2_wrp+plot_layout(ncol = 2,widths = c(1, 1))
main2=p.venn+p.sig2+plot_layout(ncol = 2,widths = c(1.1, 1))
p.pd.main=main1/main2+plot_annotation(tag_levels = 'a') +  plot_layout( heights = c(1.3, 1))&
  theme(plot.tag = element_text(face = "bold", size = 12))
p.pd.main

# ggsave(filename = "~/manuscript/realdata/Plots/pd_main.png",
#        width = 10,
#        height = 10,
#        limitsize = FALSE,
#        plot = p.pd.main)
#

