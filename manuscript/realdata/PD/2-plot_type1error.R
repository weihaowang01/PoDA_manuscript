
# Error1 and Error2 are the results from ~/manuscript/realdata/PD/1-type1error.R
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)

#################################################################################
################################### Wallen 1 ####################################
#################################################################################
pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
CDIlist2 = readRDS("~/manuscript/realdata/PD/Wallen2.rds")
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
CDIlist = readRDS("~/manuscript/realdata/PD/Wallen1.rds")
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

##Wallen1
otutable2=otutable1[which(meta1$PD=="PD"),]
meta2=meta1[which(meta1$PD=="PD"),]
otutable3=otutable1[which(meta1$PD=="HC"),]   
meta3=meta1[which(meta1$PD=="HC"),]

otu=rbind(otutable2,otutable3)
meta0=rbind(meta2,meta3)
Y <- t(otu)



modified_string <- rownames(Y)




rejtaxon = Error1[,1]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "LinDA"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig1 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none", axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))



rejtaxon = Error1[,2]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "Wallen1",
  facet2 = "PoDA"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig2 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",  axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text.x = element_text(size = 22),  
    strip.text.y = element_text(size = 16)  
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet2 ~ facet, scales = "free", labeller = labeller(facet = label_value))




rejtaxon = Error1[,3]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "LOCOM"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig3 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",  axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))



rejtaxon = Error1[,4]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "ANCOM-BC"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig4 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",  axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))

# Display the plot
p.sig4




rejtaxon = Error1[,5]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "ANCOM-BC2"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig5 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_text(hjust=1,angle=90,vjust=0.5,size=12,face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.7)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))

# Display the plot
p.sig5




arranged_plots=ggarrange(p.sig2, p.sig1, p.sig3, p.sig4, p.sig5,
                         ncol = 1,
                         heights = c(1.2, 1, 1, 1, 2))



# ptype2=arranged_plots2

ptype <- annotate_figure(
  arranged_plots,
  left = text_grob("Type I error",
                   rot = 90,
                   vjust = 0.5,
                   hjust = 0,
                   size = 22)
)







#################################################################################
################################### Wallen 1 ####################################
#################################################################################
pre_cut=0.1
age_cut=70
age_cut2=100

###### Wallen2
CDIlist2 = readRDS("~/manuscript/realdata/PD/Wallen2.rds")
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
CDIlist = readRDS("~/manuscript/realdata/PD/Wallen1.rds")
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
# #

# otutable2=otutable1[which(meta1$PD=="PD"),]
# meta2=meta1[which(meta1$PD=="PD"),]
# otutable3=otutable1[which(meta1$PD=="HC"),]   ##Wallen1
# meta3=meta1[which(meta1$PD=="HC"),]
# print("Wallen1")
otutable2=otutable[which(meta$PD=="PD"),]
meta2=meta[which(meta$PD=="PD"),]
otutable3=otutable[which(meta$PD=="HC"),]   ##Wallen2
meta3=meta[which(meta$PD=="HC"),]

otu=rbind(otutable2,otutable3)
meta0=rbind(meta2,meta3)
grp=c(rep(1,nrow(otutable2)),rep(0,nrow(otutable3)))
Y <- t(otu)



modified_string <- rownames(Y)




rejtaxon = Error2[,1]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "LinDA"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig1 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",   axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5 ,size=7       ), axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))



rejtaxon = Error2[,2]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet1 = "Wallen2",
  facet = "PoDA"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig2 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",   axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),plot.margin = margin(0, 0, 0, 0),
    strip.text.x = element_text(size = 22),  
    strip.text.y = element_text(size = 16)  
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid( facet ~ facet1, scales = "free", labeller = labeller(facet = label_value))




rejtaxon = Error2[,3]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "LOCOM"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig3 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",   axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))



rejtaxon = Error2[,4]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "ANCOM-BC"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig4 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",   axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5     ,size=7    ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))

# Display the plot
p.sig4




rejtaxon = Error2[,5]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "ANCOM-BC2"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig5 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",   axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_text(hjust=1,angle=90,vjust=0.5,size=12,face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.7)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))

# Display the plot
p.sig5




arranged_plots2=ggarrange(p.sig2, p.sig1, p.sig3, p.sig4, p.sig5,
                          ncol = 1,
                          heights = c(1.2, 1, 1, 1, 2))

# ptype2=arranged_plots2

ptype2 <- annotate_figure(
  arranged_plots2,
  left = text_grob("Type I error",
                   rot = 90,
                   vjust = 0.5,
                   hjust = 0,
                   size = 22)
)











p.t1e.main=ggarrange(ptype,ptype2,ncol=1,
                     align = "h",labels = c("a", "b"),font.label = list(size = 22))




# ggsave(filename = "~/manuscript/realdata/Plots/pd_t1e.png",
#        width = 30,
#        height = 35,
#        limitsize = FALSE,
#        plot = p.t1e.main,bg="white")




