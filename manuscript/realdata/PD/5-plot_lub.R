# source("~/manuscript/realdata/PD/3-wallen1_detection.R")
# source("~/manuscript/realdata/PD/3-wallen2_detection.R")
# rej4 is the result from ~/PD/wallen1.R, and rej5 is the result from ~/manuscript/realdata/PD/3-wallen2_detection.R

library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)

pdcommon=taxs[intersect(rej4,rej5)]
pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
Lub6 = readRDS("~/manuscript/realdata/PD/Lub6.rds")
otutable=Lub6$otu
meta = Lub6$meta
meta$Age=as.numeric(meta$Age)
keepsample = which(meta$Age >= age_cut & meta$Age < age_cut2)
otutable=otutable[keepsample,]
meta=meta[keepsample,]


datasize = nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  pre_cut)
otutable = otutable[,keepOTUs]

pdcommon[which(pdcommon%in%colnames(otutable))]
otu_lub6=otutable[,which(colnames(otutable)%in%pdcommon)]




pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
Lub12 = readRDS("~/manuscript/realdata/PD/Lub12.rds")
otutable=Lub12$otu
meta = Lub12$meta
meta$Age=as.numeric(meta$Age)
keepsample = which(meta$Age >= age_cut & meta$Age < age_cut2)
otutable=otutable[keepsample,]
meta=meta[keepsample,]


datasize = nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  pre_cut)
otutable = otutable[,keepOTUs]


pdcommon[which(pdcommon%in%colnames(otutable))]


otu_lub12=otutable[,which(colnames(otutable)%in%pdcommon)]


filter_taxa_by_percentage = function(count_table, threshold = 0.5) {
  presence_count = colSums(count_table > 0) 
  filtered_taxa = colnames(count_table)[presence_count / nrow(count_table) > threshold]
  
  return(filtered_taxa)
}

filtered_taxa6 = filter_taxa_by_percentage(otu_lub6)
filtered_taxa12 = filter_taxa_by_percentage(otu_lub12)
# filtered_taxa6==filtered_taxa12
filtered_taxa = filtered_taxa12



##################################################################################################
########################################## Lub 6 #################################################
##################################################################################################

pdcommon=taxs[intersect(rej4,rej5)]
pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
Lub6 = readRDS("~/manuscript/realdata/PD/Lub6.rds")
otutable=Lub6$otu
meta = Lub6$meta
meta$Age=as.numeric(meta$Age)
keepsample = which(meta$Age >= age_cut & meta$Age < age_cut2)
otutable=otutable[keepsample,]
meta=meta[keepsample,]


datasize = nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  pre_cut)
otutable = otutable[,keepOTUs]






pdcommon[which(pdcommon%in%colnames(otutable))]




u=as.numeric(factor(meta$PD))-1
Y=t(otutable)
z1=as.numeric(factor(meta$sex))
u[which(z1==2)]=u[which(z1==2)]+2









hh = match(filtered_taxa,pdcommon)
plot_list = list()

for (i in 1:19) {
  kkk = which(colnames(otutable) == pdcommon[hh[i]])
  dataset = data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                       group = factor(u),facet2=pdcommon[hh[i]])
  
  colnames(dataset) = c("value", "group","facet2")
  
  p = ggplot(dataset, aes(x = group, y = value, color = group)) +
    theme_bw() +
    geom_boxplot() +
    scale_x_discrete(labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")) +
    scale_color_manual(values = c("#3D5C6F", "#E47159", "#3D5C6F", "#E47159")) +
    facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    ) +
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
    labs(x = "Group", y = "Relative abundance")
  
  plot_list[[i]] = p
}


plot_list

for (i in c(7,14,19)) {
  kkk = which(colnames(otutable) == pdcommon[hh[i]])
  dataset = data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                       group = factor(u),facet2=pdcommon[hh[i]],facet1="Lub6")
  
  colnames(dataset) = c("value", "group","facet2","facet1")
  p = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
    geom_boxplot() +
    facet_grid(facet1~facet2  , labeller = labeller(facet2 = label_value))+
    scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank()  , axis.text.y = element_text(
            angle = 90,
            hjust = 0.5,
            vjust = 0.5        ))+
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
    labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)
  
  plot_list[[i]] = p
}




##################################################################################################
######################################### Lub 12 #################################################
##################################################################################################

pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
Lub12 = readRDS("~/manuscript/realdata/PD/Lub12.rds")
otutable=Lub12$otu
meta = Lub12$meta
meta$Age=as.numeric(meta$Age)
keepsample = which(meta$Age >= age_cut & meta$Age < age_cut2)
otutable=otutable[keepsample,]
meta=meta[keepsample,]


datasize = nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  pre_cut)
otutable = otutable[,keepOTUs]






pdcommon[which(pdcommon%in%colnames(otutable))]




u=as.numeric(factor(meta$PD))-1
Y=t(otutable)
z1=as.numeric(factor(meta$sex))
u[which(z1==2)]=u[which(z1==2)]+2



hh = match(filtered_taxa,pdcommon)
plot_list2 = list()

for (i in 1:19) {
  kkk = which(colnames(otutable) == pdcommon[hh[i]])
  
  dataset = data.frame(
    value = c(t(as.matrix(Y[kkk, ])) / colSums(Y)),
    group = factor(u)
  )
  colnames(dataset) = c("value", "group")
  
  p = ggplot(dataset, aes(x = group, y = value, color = group)) +
    theme_bw() +
    geom_boxplot() +
    scale_x_discrete(labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")) +
    scale_color_manual(values = c("#3D5C6F", "#E47159", "#3D5C6F", "#E47159")) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    ) +
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
    labs(x = "Group", y = "Relative abundance")
  
  plot_list2[[i]] = p
}


for (i in c(7,14,19)) {
  kkk = which(colnames(otutable) == pdcommon[hh[i]])
  dataset = data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                       group = factor(u),facet1="Lub12")
  
  colnames(dataset) = c("value", "group","facet1")
  p = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
    geom_boxplot() +
    scale_x_discrete(  
      labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
    ) +
    scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank()  , axis.text.y = element_text(
            angle = 90,
            hjust = 0.5,
            vjust = 0.5        ))+
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
    facet_grid(facet1~.  , labeller = labeller(facet2 = label_value))+
    labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)
  
  
  plot_list2[[i]] = p
}





library(gridExtra)

combined_plots = list()

for (i in 1:19) {
  combined_plots[[i]] = grid.arrange(plot_list[[i]], plot_list2[[i]], ncol = 1)
}

plub=grid.arrange(
  grobs = combined_plots, 
  ncol = 7, 
  nrow = 3, 
  layout_matrix = rbind(
    c(1, 2, 3, 4, 5, 6, 7),
    c(8, 9, 10, 11, 12, 13, 14),
    c(15, 16, 17, 18, 19,NA, NA)  
  )
)

plub = annotate_figure(
  plub,
  left = text_grob("Relative abundance",
                   rot = 90,
                   # vjust = 1,
                   # hjust = 1,
                   size = 12)
)

# ggsave(filename = "~/manuscript/realdata/Plots/pd_lub.png",
#        width = 16,
#        height = 20,
#        limitsize = FALSE, bg = "white",
#        plot = plub)

