# Here, we only consider three methods that need to handle zeros: PoDA, LinDA, and ANCOM-BC.

library(dplyr)
library(tidyr)
library(tibble)

source("~/manuscript/linda_modify.R")
source("~/manuscript/poda_modify.R")
source("~/manuscript/helpler_function.R")
source("~/manuscript/PODA.R")
load("~/manuscript/realdata/CDI/CDI.RData")

ind <- data.obj$meta.dat$disease_stat %in% c('Case', 'DiarrhealControl')
CDI.otu <- as.data.frame(data.obj$otu.tab[, ind])
CDI.otu2=CDI.otu
CDI.meta =data.frame(data.obj$meta.dat[ind,'disease_stat'])
rownames(CDI.meta )=rownames(data.obj$meta.dat[ind,])

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

################## zero-replacement strategies ##################
# # pseudo=count
# Y=Y+0.5 # 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1

# multRepl
tY=t(Y)
Y=t(multRepl(tY,label=0,z.delete = FALSE))
colnames(Y)=rownames(tY)
# # multLN
# tY=t(Y)
# Y=t(multLN(tY,label=0,z.delete = FALSE))

# # adaImpute
# Y<-imputation(Y,Z,formulaa)

# poda
rej=PODA1(Y = Y, u = u, confoun = NULL ,rep_time=500,fdrnomial=0.05 ,lib_cut=0, pre_cut=0)
rej=rej$rej
rej
# 0.1: 15 17
# 0.2-1 & multRepl & multLN & adaImpute: 15

# linda
res <- lindaa(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = 0)
rej.ld <- which(res$output[[1]]$reject)
rej.ld
# 0.1-0.5 & adaImpute & multRepl: 15 17
# 0.6-1 & multLN: 15

#ancombc
# Since ANCOM-BC automatically handles zeros by adding 1 to the count matrix,
# the microbial count matrix needs to be subtracted by 1 before input here.
m=nrow(Y)
rej.ancombc <- ancombc1.fun(m,Y-1, Z, formula, alpha=0.05,p_adj_method = "BH")
rej.ancombc
# 0.1-0.2: 3  5 15 20 26 32 34 40 41 45
# 0.3: 3 15 20 26 32 34 40 41 45
# 0.4: 15 32
# 0.5-1 & adaImput & multLN & multRepl: 15



