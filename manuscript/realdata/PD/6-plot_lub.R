# source("~/manuscript/realdata/PD/3-wallen1_detection.R")
# source("~/manuscript/realdata/PD/3-wallen2_detection.R")
# rej4 is the result from ~/PD/wallen1.R, and rej5 is the result from ~/manuscript/realdata/PD/3-wallen2_detection.R

library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)

xixi3=taxs[intersect(rej4,rej5)]
pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
CDIlist2 = readRDS("~/manuscript/realdata/PD/Lub6.rds")
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






xixi3[which(xixi3%in%colnames(otutable))]




u=as.numeric(factor(meta$PD))-1
Y=t(otutable)
z1=as.numeric(factor(meta$sex))
u[which(z1==2)]=u[which(z1==2)]+2




hh = 1
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p11 = ggplot(dataset, aes(x = group, y = value, color = group))  + theme_bw()+
  geom_boxplot() +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

p11




hh = 2
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p21 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p21


hh = 3
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p31 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p31

hh = 4
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p41 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p41


hh = 5
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p51 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p51



hh = 6
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p61 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p61


hh = 7
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh],facet1="Lub6")

colnames(dataset) = c("value", "group","facet2","facet1")
p71 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
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

p71


hh = 8
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p81 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p81


hh = 10
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p101 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p101


hh = 11
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p111 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p111


hh = 12
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p121 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p121


hh = 14
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p141 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
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

p141


hh = 15
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh])

colnames(dataset) = c("value", "group","facet2")
p151 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(.~facet2  , labeller = labeller(facet2 = label_value))+
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

p151


hh = 16
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet2=xixi3[hh],facet1="Lub6")

colnames(dataset) = c("value", "group","facet2","facet1")
p161 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  facet_grid(facet1~facet2  , labeller = labeller(facet2 = label_value))+
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

p161






pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
CDIlist2 = readRDS("~/manuscript/realdata/PD/Lub12.rds")
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






xixi3[which(xixi3%in%colnames(otutable))]




u=as.numeric(factor(meta$PD))-1
Y=t(otutable)
z1=as.numeric(factor(meta$sex))
u[which(z1==2)]=u[which(z1==2)]+2




hh = 1
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p110 = ggplot(dataset, aes(x = group, y = value, color = group))  + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

p110



hh = 2
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p210 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
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
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

p210


hh = 3
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p310 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
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
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)



hh = 4
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p410 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
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
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)




hh = 5
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p510 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(  # 修改 x 轴刻度标签
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
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)





hh = 6
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p610 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
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
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)




hh = 7
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet1="Lub12")

colnames(dataset) = c("value", "group","facet1")
p710 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(  # 修改 x 轴刻度标签
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

p710


hh = 8
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p810 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

p810


hh = 10
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p1010 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)



hh = 11
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p1110 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(  # 修改 x 轴刻度标签
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)



# 你的数据准备
hh = 12
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p1210 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(  # 修改 x 轴刻度标签
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")  # 按 group 的原始顺序对应标签
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)




hh = 14
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p1410 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)

# p141


hh = 15
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u))

colnames(dataset) = c("value", "group")
p1510 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  # geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)




hh = 16
kkk = which(colnames(otutable) == xixi3[hh])
dataset <- data.frame(value = c(t(as.matrix(Y[kkk,]))/colSums(Y)),
                      group = factor(u),facet1="Lub12")

colnames(dataset) = c("value", "group","facet1")
p1610 = ggplot(dataset, aes(x = group, y = value, color = group)) + theme_bw()+
  geom_boxplot() +
  scale_x_discrete(  # 修改 x 轴刻度标签
    labels = c("Male-HC", "Male-PD", "Female-HC", "Female-PD")
  ) +
  scale_color_manual(values = c("#3D5C6F", "#E47159","#3D5C6F", "#E47159")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6  ),
        axis.title.y = element_blank()  , axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5        ))+
  stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
  facet_grid(facet1~.  , labeller = labeller(facet2 = label_value))+
  labs( x = "Group", y = "Relative abundance")#+ylim(lima,limb)










plub=ggarrange(p11,p21,p31,p41,p51,p61,p71,
              p110,p210,p310,p410,p510,p610,p710,
              p81,p101,p111,p121,p141,p151,p161,
              p810,p1010,p1110,p1210,p1410,p1510,p1610,
              ncol=7,nrow=4)

plub <- annotate_figure(
  plub,
  left = text_grob("Relative abundance",
                   rot = 90,
                   # vjust = 1,
                   # hjust = 1,
                   size = 12)
)
plub



ggsave(filename = "~/manuscript/realdata/Plots/plub.png",
       width = 16,
       height = 18,
       limitsize = FALSE, bg = "white",
       plot = plub)
