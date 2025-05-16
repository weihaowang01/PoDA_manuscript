load("/Users/wendywang/Downloads/throat.meta (1).RData")
####################################################################
####### estimate beta0 by the method proposed by Zhou et al. #######
####################################################################
para.fun <- function(otu.tab, m, D = NULL, model = c('log-normal', 'Gamma')) {
  otu.tab <- as.matrix(otu.tab)
  colnames(otu.tab) = NULL
  rownames(otu.tab) = NULL
  
  has.read <- rowSums(otu.tab > 0)
  ind.taxa <- sort(order(has.read, decreasing = TRUE)[1 : m])
  otu.tab.sel <- otu.tab[ind.taxa, ]
  
  if(model == 'log-normal'){
    W <- D %*% log(otu.tab.sel+0.1)
    n <- ncol(otu.tab.sel)
    
    beta.diff.est <- rowMeans(W)
    sigma2.add.est <- rowSums((W - beta.diff.est) ^ 2) / (n - 1)
    fit <- cv.glmnet(D, beta.diff.est, alpha = 1, standardize = FALSE, intercept = FALSE)
    beta.est <- coef(fit, s = 'lambda.min')[-1]
    
    s <- sum(sigma2.add.est) / (m - 1)
    ind <- t(t(matrix(Matrix::which(D != 0), nrow = m - 1)) - 
               (0 : (m - 1)) * m * (m - 1) /  2)
    tmp <- matrix(sigma2.add.est[ind], nrow = m - 1)
    sigma2.est <- (colSums(tmp) - s) / (m - 2)
    out <- list(otu.tab.sel = otu.tab.sel, beta0 = beta.est, sigma2 = sigma2.est)
  } else if(model == 'Gamma') {
    res <- dirmult(t(otu.tab.sel))
    pi0 <- res$pi
    theta0 <- res$theta
    eta0 <- res$gamma
    out <- list(otu.tab.sel = otu.tab.sel, pi0 = pi0, theta0 = theta0, eta0 = eta0)
  }
  return(out)
}

D.fun <- function (m) {
  D <- foreach(k = 1 : (m - 1), .combine = 'rbind') %do% {
    
    i <- rep(1 : (m - k), 2)
    j <- c(rep(k, m - k), (k + 1) : m)
    x <- c(rep(1, m - k), rep(-1, m - k))
    
    sparseMatrix(i, j, x = x, dims = c(m - k, m))
    
  }
  return(D)
}

####################################################################
####### estimate beta0 by the method proposed by Cao et al. ########
####################################################################

coat <- function(x, nFoler = 5, soft = 1){
  startTime <- proc.time()
  p <- ncol(x)
  clrX <- log(x) - rowSums(log(x)) %*%matrix(1,1,p) / p
  coatPred <- adaptThresoldCov(clrX, soft = soft)
  sigma <- coatPred$sigma
  corr <- coatPred$corr
  exeTimeClass <- proc.time() - startTime
  exeTime <- as.numeric(exeTimeClass[3])
  return(list(sigma = sigma, corr = corr, time = exeTime))
}



adaptThresoldCov <- function(x, nFolder = 5, soft = 1){
  n <- nrow(x)
  p <- ncol(x)
  nGrid <- 100
  gridInfo <- adaptThresholdRange(x)
  grid <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nGrid)/nGrid
  part <- 1 + sample(c(1:n))%%nFolder
  error <- matrix(0, nFolder, nGrid)
  for (i in 1:nFolder){
    xTest <- x[which(part == i),]
    xTrain <- x[which(part != i),]
    gridInfoTrain <- adaptThresholdRange(xTrain)
    covTest <- cov(xTest)*(n-1)/n
    for (j in 1:nGrid){
      sigmaTrain <- adaptThreshold(gridInfoTrain$cov,gridInfoTrain$theta,grid[j],soft)
      error[i,j] <- (norm(sigmaTrain-covTest, "F"))
    }
  }
  errorSum <- colSums(error)
  lambda <- grid[which(errorSum == min(errorSum))][1]
  sigma <- adaptThreshold(gridInfo$cov,gridInfo$theta,lambda,soft)
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))
  return(list(sigma = sigma, corr = corr))
}


adaptThresholdRange <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  cov <- cov(x)*(n-1)/n
  centered.x <- scale(x, scale = FALSE)
  theta <- (t(centered.x)^2)%*%(centered.x^2)/n - cov^2
  delta <- cov/(theta^0.5)
  delta <- abs(delta - diag(diag(delta)))
  upper <- max(delta)
  lower <- min(delta[which(delta != 0)])
  return(list(upper = upper, lower = lower, theta = theta, cov = cov))
}

adaptThreshold <- function(cov,theta,lambda,soft){
  covOffDiag <- cov - diag(diag(cov))
  thetaOffDiag <- theta - diag(diag(theta))
  sigmaTmp <- abs(covOffDiag) - lambda*thetaOffDiag^0.5
  sigmaTmp[which(sigmaTmp < 0)] <- 0
  if (soft == 1){
    sigma <- diag(diag(cov)) + sigmaTmp*sign(covOffDiag)
  }else{
    sigma <- cov
    sigma[which(sigmaTmp < 1e-10)] <- 0
    sigma <- sigma + diag(diag(cov))
  }
  return(sigma)
}


data(throat.otu.table, package = "LOCOM")
otu.tab=t(throat.otu.table)
otu.tab <- as.matrix(otu.tab)
otu.tab=otu.tab[,which(throat.meta$SmokingStatus=="NonSmoker")]
colnames(otu.tab) = NULL
rownames(otu.tab) = NULL


m <- 500
D <- D.fun(m)

res <- para.fun(otu.tab, m, D, model = 'log-normal') 
log.normal.para <- list(beta0 = res[[2]], sigma2 = res[[3]])



data(throat.otu.table, package = "LOCOM")
otu.tab=t(throat.otu.table)
otu.tab <- as.matrix(otu.tab)
otu.tab=otu.tab[,which(throat.meta$SmokingStatus=="NonSmoker")]
colnames(otu.tab) = NULL
rownames(otu.tab) = NULL



has.read <- rowSums(otu.tab > 0)
m=500
ind.taxa <- sort(order(has.read, decreasing = TRUE)[1 : m])
otu.tab.sel <- otu.tab[ind.taxa, ]
otu.tab.sel[which(otu.tab.sel==0)]=0.5 
otu.tab.sel.t = t(otu.tab.sel)
p <- ncol(otu.tab.sel.t )
otu.tab.sel.t  <- otu.tab.sel.t /(rowSums(otu.tab.sel.t )%*%matrix(1,1,p))
tmpsigma = coat(otu.tab.sel.t, nFoler = 5, soft = 1)
tmpsigma$sigma[1:6,1:6]
which(tmpsigma$sigma!=0)


eigen(tmpsigma$sigma)$value
which(tmpsigma$sigma!=0)
eigen(tmpsigma$sigma)$value



###### modify the eigenvalues #######
eigen_result <- eigen(tmpsigma$sigma)
eigen_result$values[eigen_result$values < 0] <- 1e-10
corr_modified <- eigen_result$vectors %*% diag(eigen_result$values) %*% t(eigen_result$vectors)
eigen(corr_modified)$value




beta0 = log.normal.para$beta0
sigma2 = diag(tmpsigma$sigma)
para=list(beta0=beta0,sigma2=sigma2)

saveRDS(para, "~/simulation/lognormal.urt2.rds")
saveRDS(tmpsigma$sigma,"~/simulation/correlation_non.rds")
saveRDS(corr_modified,"~/simulation/correlation_non_modify.rds")
