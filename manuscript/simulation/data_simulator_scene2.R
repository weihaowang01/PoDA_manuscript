# a data generation mechanism similar to that of Lin et al.
library(LOCOM)
data(throat.otu.table, package = "LOCOM")
prevalence = apply(t(throat.otu.table), 1, function(x)
  sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
tax_keep = which(prevalence >= 0.05)
n = nnum
d = length(tax_keep)
m=d
diff_prop = gammain
iter_num = 1
seed = seq_len(iter_num)
df_sim_params = data.frame(expand.grid(n, diff_prop, seed)) %>%
  dplyr::rename(n = Var1, diff_prop = Var2, seed = Var3) %>%
  arrange(n, diff_prop, seed)
list_sim_params = apply(df_sim_params, 1, paste0, collapse = "_")
lfc_value<-seq(0.5,2,by=0.5)
# lfc_value<-c(-1,-2,1,2)

lfc_bin_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  
  set.seed(q*100+i)
  
  lfc_bin_list[[i]] = sample(c(0, lfc_value), size = d, replace = TRUE,
                             prob = c(1 - diff_prop[i], 
                                      rep(diff_prop[i]/length(lfc_value), length(lfc_value))))
}
names(lfc_bin_list) = diff_prop

# Log-fold-changes for the continuous confounder
lfc_cont_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  set.seed(q*100+i)
  
  lfc_cont_list[[i]] = sample(c(0, 1), size = d, replace = TRUE,
                              prob = c(1 - diff_prop[i], diff_prop[i]))
}
names(lfc_cont_list) = diff_prop


i = list_sim_params
params = strsplit(i, "_")[[1]]

n = as.numeric(params[1])
diff_prop = as.numeric(params[2])

set.seed(q)
abn_data = ANCOMBC::sim_plnm(abn_table = throat.otu.table, taxa_are_rows = FALSE,
                             prv_cut = 0.05, n = n, lib_mean = 1e8, disp = 0.5)
log_abn_data = log(abn_data + 1e-5)
rownames(log_abn_data) = paste0("taxon", seq_len(d))
colnames(log_abn_data) = paste0("sample", seq_len(n))

set.seed(q)
smd = data.frame(sample = paste0("sample", seq_len(n)),
                 ### setting A8
                 # samp_frac = log(c(runif(n/2, min = 1e-3, max = 1e-2), #us=0
                 #                   runif(n/2, min = 1e-3, max = 1e-2))), #us=1
                 ## setting A9
                 # samp_frac = log(c(runif(n/2, min = 1e-4, max = 1e-3), #us=0
                 #                   runif(n/2, min = 1e-4, max = 1e-3))), #us=1
                 # ### setting A10
                 samp_frac = log(c(runif(n/2, min = 5e-4, max = 1e-3),#us=0
                                   runif(n/2, min = 1e-3, max = 5e-3))), #us=1
                 cont_cov = rnorm(n),   
                 bin_cov = as.factor(rep(seq_len(2), each = n/2)))

d = nrow(abn_data) 
lfc_bin = lfc_bin_list[[as.character(diff_prop)]]
lfc_cont = lfc_cont_list[[as.character(diff_prop)]]
set.seed(q)
fmd = data.frame(taxon = paste0("T", seq_len(d)),
                 seq_eff = log(runif(d, min = 0.1, max = 1)),
                 # seq_eff = 0,
                 lfc_cont = lfc_cont,
                 lfc_bin = lfc_bin)

# Add effect sizes of covariates to the true abundances
smd_dmy = model.matrix(~ 0 + cont_cov + bin_cov, data = smd)
log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
log_abn_data = log_abn_data + outer(fmd$lfc_bin, smd_dmy[, "bin_cov2"])
true_abundance = exp(log_abn_data)    #microbioload




# Add sample- and taxon-specific biases
log_otu_data = t(t(log_abn_data) + smd$samp_frac)
log_otu_data = log_otu_data + fmd$seq_eff
otu_data = round(exp(log_otu_data))

gamma=diff_prop 
Y = otu_data
u = as.matrix(as.numeric(smd[,4])-1)
Z = as.numeric(as.matrix((smd[,-(1:2)])))
Z = matrix(Z,ncol=2)
temp = Z[,1]
Z[,1] = u
Z[,2] = temp
ZZ=data.frame(Z)
ZZ[,1]=factor(ZZ[,1])
colnames(Z)<-c("u","z1")
colnames(ZZ)<-c("u","z1")
confoun<-as.matrix(Z[,-1])
z1=confoun
colnames(confoun)<-"z1"
H<-rep(0,m)
H[which(fmd$lfc_bin!=0)]<-1
alpha = lfc_bin_list[[1]]