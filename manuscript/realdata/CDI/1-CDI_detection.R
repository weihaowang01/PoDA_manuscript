library(dplyr)
library(tidyr)
library(tibble)
library(parallel)
library(doParallel)
library(ggplot2)
library(ggVennDiagram)
library(ggvenn)
library(LinDA)
library(mvtnorm)
library(foreach)
library(LOCOM)
library(microbiome)
library(ANCOMBC)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(tidytree)

source("~/manuscript/helper_function.R")
source("~/manuscript/PODA.R")
load("~/manuscript/realdata/CDI/CDI.RData")

ind <- data.obj$meta.dat$disease_stat %in% c('Case', 'DiarrhealControl')
CDI.otu <- as.data.frame(data.obj$otu.tab[, ind])
CDI.otu2=CDI.otu
CDI.meta =data.frame(data.obj$meta.dat[ind,'disease_stat'])
rownames(CDI.meta)=rownames(data.obj$meta.dat[ind,])

cdi.name=data.obj$otu.name[,6]




otu_long <- CDI.otu2 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Otu") %>%
  pivot_longer(-Otu, names_to = "DA", values_to = "Abundance") %>%
  mutate(Genus = cdi.name[match(Otu, rownames(CDI.otu2))])



merged_otu <- otu_long %>%
  group_by(Genus, DA) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = "DA", values_from = "Abundance") %>%
  column_to_rownames(var = "Genus")


print(dim(merged_otu))


preprocess.fun <- function(otu.tab, meta, prev.cut = 0, lib.cut = 1000,
                           winsor.quan = 0.97) {
  keep.sam <- colSums(otu.tab) >= lib.cut
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(meta[keep.sam, ])
  names(Z) <- names(meta)
  rownames(Z) <- rownames(meta)[keep.sam]
  keep.tax <- rowSums(Y > 0) / ncol(Y) >= prev.cut
  Y <- Y[keep.tax, ]

  return(list(Y = Y, Z = Z, keep.sam = keep.sam, keep.tax = keep.tax))
}


#################################
########### detection ###########
#################################

res = preprocess.fun(merged_otu,CDI.meta,prev.cut = 0.1, lib.cut = 1000, winsor.quan = 1)
Y.genus <- res$Y
Z <- res$Z

Y=Y.genus
u=as.matrix(Z)



u[which(u=="DiarrhealControl"),]=0
u[which(u=="Case"),]=1
u=as.numeric(u)
Z=as.data.frame(factor(u))
colnames(Z)="u"
formula <- 'u'
formulaa<-'~u'

colnames(Z)='u'

# poda
rej=PoDA::PODA(Y = Y, u = u, confoun = NULL ,rep_time=500,fdrnomial=0.05 ,lib_cut=0, pre_cut=0)
rej=rej$rej
rej

# linda
res <- LinDA::linda(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = 0)
rej.ld <- which(res$output[[1]]$reject)
rej.ld

# locom
tt = locom(otu.table = t(Y), Y = u, C = NULL, fdr.nominal = 0.05,seed = 1234, n.cores = 10,filter.thresh = 0)
rej.locom =  which(colnames(tt$q.otu)%in%tt$detected.otu)
rej.locom

#ancombc
m=nrow(Y)
rej.ancombc <- ancombc1.fun(m,Y, Z, formula, alpha=0.05,p_adj_method = "BH")
rej.ancombc

#ancombc2
rej.ancombc2 <- ancombc2.fun(m,Y, Z, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc2


#################################
########### Venn Plot ###########
#################################


venn_list <- setNames(
  list(
    as.character(rej.ancombc2),
    as.character(rej.ld),
    as.character(rej) #rej

  ),
  c(" "," ", " ")
)

venn_data <- process_data(Venn(venn_list))
venn_data$region$count <- sapply(venn_data$region$element, length)


p=ggVennDiagram(venn_list, label = "count") +
  scale_fill_gradient(low = "white", high = "#2a83a2") +
  theme(legend.position = "none", text = element_text(size = 1))


p_fin <- p +
  annotate("text", x = -4.9, y = 2.8, label = "ANCOM-BC2", size = 3) +
  annotate("text", x = 8.2, y = 2.8,, label = "LinDA", size = 3) +
  annotate("text", x = 2.1, y = -9.4, label = "PoDA&\nANCOM-BC&\nLOCOM", size = 3)


# p_fin1 <- wrap_elements( p_fin)


p_fin11 <- p_fin + coord_cartesian(xlim = c(-6, 9.6), ylim = c(-15, 9))

p_fin11
