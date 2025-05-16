# a data generation mechanism similar to that of Zhou et al.
corr_matrix = readRDS("~/simulation/correlation_non_modify.rds")
parameter<-readRDS("~/simulation/lognormal.urt.rds")
beta0<-parameter$beta0
sigma2<-parameter$sigma2

m<-500 
n<-nnum
gamma<-gammain
mu=mumu

set.seed(q)
X0 <- matrix(exp(rnorm(m * n, beta0, sqrt((sigma2)))), nrow = m)
set.seed(q) 
taxon_frac = runif(m, min = 0.1, max = 1)
pi0 <- t(t(X0) / colSums(X0))
pi0.ave <- rowMeans(pi0)
if (any(pi0.ave == 0)) {
  ind <- which(pi0.ave == 0)
  pi0.ave[ind] <- min(pi0.ave[-ind]) / 10
  pi0.ave <- pi0.ave / sum(pi0.ave)
}
ratio=0.01
tmp <- (pi0.ave > ratio)
mu <- 1* mu * (n <= 50) + 1*mu * (n > 50)
mu.1 <- log(mu * tmp + mu * (ratio / pi0.ave) ^ (1 / 4) * (1 - tmp))

set.seed(q)
H <- rbinom(m, 1, gamma)
alpha <- mu.1 * H

#### setting A3 for strong compositional effects
# H=rep(0,m)
# abun <- order(pi0.ave, decreasing = TRUE)[1 : floor(m * 0.25)]
# set.seed(q)
# H[sample(abun, m * gamma)] <- 1
# alpha <- mu.1 * H

# ######## setting A1-A3&A5-A7 #########
# set.seed(q) 
# u <- rbinom(n, 1, 0.5) #A1-A3&A5
# # u <- rbinom(n, 1, 0.8) #A6
# # u = runif(n,-1,1) #A7
# confoun = as.matrix(rep(1,n))
# Z <- as.matrix(u)
# beta <- alpha
# colnames(Z)<-'u'




############## setting A4 ########
set.seed(q)  
z1 <- rbinom(n, 1, 0.5)
z1[which(z1 == 0)] <- -1
set.seed(q) 
z2 <- rnorm(n, 0, 1)
set.seed(q) 
beta1 <- rnorm(m, 1, 1)
set.seed(q) 
beta2 <- rnorm(m, 2, 1)
set.seed(q) 
u <- rbinom(n, 1, 1 / ( 1+exp(-0.5 * z1 - 0.5 * z2)))
confoun<-cbind(z1,z2)
Z <- cbind(u, confoun)
beta <- cbind(alpha, beta1, beta2)








set.seed(q) ##############################################
tmp <- beta %*% t(Z) + beta0
X <- matrix(exp(rnorm(m * n, tmp, rep(sqrt((sigma2)), n))), nrow = m)

# # ########### setting A5 for correlated taxa ################
# tmp <- beta %*% t(Z) + beta0
# eigen_corr=eigen(corr_matrix)
# vl = eigen_corr$values
# Sig <- eigen_corr$vectors %*% diag(sqrt(vl)) %*% t(eigen_corr$vectors)
# set.seed(q)
# X <- exp(Sig %*% matrix(rnorm(m * n), nrow = m) + tmp)
# ##########################################################

X=X*taxon_frac
pii <- t(t(X) / colSums(X))
set.seed(q) 
N <- 1*rnbinom(n, size = 6.5, mu = 8000)
# set.seed(q) 
# N[which(u==1)]=1*rnbinom(length(which(u==1)), size = 6.5, mu = 48000) # setting A2 different library size 


set.seed(q) 
Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pii[, s]))
Z <- as.data.frame(Z)
Y <- as.data.frame(Y)
colnames(Y) <- rownames(Z) <- paste0('sample', 1 : n)
rownames(Y) <- paste0('taxon', 1 : m)
