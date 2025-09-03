# In this modified function, we have turned off all of Linda's automatic zero-filling functions.

library(modeest)
lindaa=function (otu.tab, meta, formula, type = "count", adaptive = TRUE, 
                 imputation = FALSE, pseudo.cnt =0, corr.cut = 0.1, p.adj.method = "BH", 
                 alpha = 0.05, prev.cut = 0, lib.cut = 1, winsor.quan = NULL, 
                 n.cores = 1) 
{
  if (any(is.na(otu.tab))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- all.vars(as.formula(formula))
  Z <- as.data.frame(meta[, allvars])
  keep.sam <- which(colSums(otu.tab) >= lib.cut & rowSums(is.na(Z)) == 
                      0)
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  names(Z) <- allvars
  n <- ncol(Y)
  keep.tax <- which(rowSums(Y > 0)/n >= prev.cut)
  Y <- Y[keep.tax, ]
  m <- nrow(Y)
  if (any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }
  ind <- sapply(1:ncol(Z), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])
  if (!is.null(winsor.quan)) {
    Y <- winsor.fun(Y, winsor.quan)
  }
  if (grepl("\\(", formula)) {
    random.effect <- TRUE
  }
  else {
    random.effect <- FALSE
  }
  if (is.null(rownames(otu.tab))) {
    taxa.name <- (1:nrow(otu.tab))[keep.tax]
  }
  else {
    taxa.name <- rownames(otu.tab)[keep.tax]
  }
  if (is.null(rownames(meta))) {
    samp.name <- (1:nrow(meta))[keep.sam]
  }
  else {
    samp.name <- rownames(meta)[keep.sam]
  }
  if (type == "count") {
    if (any(Y == 0)) {
      N <- colSums(Y)
      if (adaptive) {
        logN <- log(N)
        if (random.effect) {
          tmp <- lmer(as.formula(paste0("logN", formula)), 
                      Z)
        }
        else {
          tmp <- lm(as.formula(paste0("logN", formula)), 
                    Z)
        }
        corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
        if (any(corr.pval <= corr.cut)) {
          cat("Imputation approach is used.\n")
          imputation <- FALSE
        }
        else {
          cat("Pseudo-count approach is used.\n")
          imputation <- FALSE
        }
      }
      if (imputation) {
        N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
        N.mat[Y > 0] <- 0
        tmp <- N[max.col(N.mat)]
        Y <- Y + N.mat/tmp
      }
      else {
        Y <- Y + pseudo.cnt
      }
    }
  }
  if (type == "proportion") {
    if (any(Y == 0)) {
      Y <- t(apply(Y, 1, function(x) {
        x[x == 0] <- 0.5 * min(x[x != 0])
        return(x)
      }))
    }
  }
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)
  oldw <- getOption("warn")
  options(warn = -1)
  if (!random.effect) {
    suppressMessages(fit <- lm(as.formula(paste0("W", formula)), 
                               Z))
    res <- do.call(rbind, coef(summary(fit)))
    d <- ncol(model.matrix(fit))
    df <- rep(n - d, m)
    tmp <- vcov(fit)
    res.cov <- foreach(i = 1:m) %do% {
      tmp[((i - 1) * d + 1):(i * d), ((i - 1) * d + 1):(i * 
                                                          d)]
    }
    wald <- list(beta = coef(fit), sig = sigma(fit), X = model.matrix(fit))
    res.cov <- do.call(rbind, res.cov)
    rownames(res.cov) <- rownames(res)
    colnames(res.cov) <- rownames(res)[1:d]
  }
  else {
    fun <- function(i) {
      w <- W[, i]
      fit <- lmer(as.formula(paste0("w", formula)), Z)
      a <- as_lmerModLmerTest(fit)
      rand.cov <- a@vcov_varpar
      Jac.beta.cov.rand <- a@Jac_list
      list(coef(summary(fit)), vcov(fit), rand.cov, Jac.beta.cov.rand)
    }
    if (n.cores > 1) {
      tmp <- mclapply(c(1:m), function(i) fun(i), mc.cores = n.cores)
    }
    else {
      suppressMessages(tmp <- foreach(i = 1:m) %do% fun(i))
    }
    res <- do.call(rbind, lapply(tmp, `[[`, 1))
    res.cov <- do.call(rbind, lapply(tmp, `[[`, 2))
    wald <- list(beta = do.call(cbind, lapply(lapply(tmp, 
                                                     `[[`, 1), function(x) x[, 1])), beta.cov = lapply(tmp, 
                                                                                                       `[[`, 2), rand.cov = lapply(tmp, `[[`, 3), Jac.beta.cov.rand = lapply(tmp, 
                                                                                                                                                                             `[[`, 4))
  }
  options(warn = oldw)
  res.intc <- res[which(rownames(res) == "(Intercept)"), ]
  rownames(res.intc) <- NULL
  options(warn = -1)
  suppressMessages(bias.intc <- mlv(sqrt(n) * res.intc[, 1], 
                                    method = "meanshift", kernel = "gaussian")/sqrt(n))
  options(warn = oldw)
  baseMean <- 2^(res.intc[, 1] - bias.intc)
  baseMean <- baseMean/sum(baseMean) * 1e+06
  output.fun <- function(x) {
    res.voi <- res[which(rownames(res) == x), ]
    rownames(res.voi) <- NULL
    if (random.effect) {
      df <- res.voi[, 3]
    }
    log2FoldChange <- res.voi[, 1]
    lfcSE <- res.voi[, 2]
    oldw <- getOption("warn")
    options(warn = -1)
    suppressMessages(bias <- mlv(sqrt(n) * log2FoldChange, 
                                 method = "meanshift", kernel = "gaussian")/sqrt(n))
    options(warn = oldw)
    log2FoldChange <- log2FoldChange - bias
    stat <- log2FoldChange/lfcSE
    pvalue <- 2 * pt(-abs(stat), df)
    padj <- p.adjust(pvalue, method = p.adj.method)
    reject <- padj <= alpha
    output <- cbind.data.frame(baseMean, log2FoldChange, 
                               lfcSE, stat, pvalue, padj, reject, df)
    rownames(output) <- taxa.name
    return(list(bias = bias, output = output))
  }
  cov.fun <- function(x) {
    tmp <- (1:ncol(res.cov))[-c(1, which(colnames(res.cov) == 
                                           x))]
    covariance <- as.data.frame(as.matrix(res.cov[which(rownames(res.cov) == 
                                                          x), tmp]))
    rownames(covariance) <- taxa.name
    colnames(covariance) <- colnames(res.cov)[tmp]
    return(covariance)
  }
  variables <- unique(rownames(res))[-1]
  variables.n <- length(variables)
  bias <- rep(NA, variables.n)
  output <- list()
  if (variables.n == 1) {
    covariance <- NULL
  }
  else {
    covariance <- list()
  }
  for (i in 1:variables.n) {
    tmp <- output.fun(variables[i])
    output[[i]] <- tmp[[2]]
    bias[i] <- tmp[[1]]
    if (variables.n > 1) {
      covariance[[i]] <- cov.fun(variables[i])
    }
  }
  names(output) <- variables
  if (variables.n > 1) {
    names(covariance) <- variables
  }
  rownames(Y) <- taxa.name
  colnames(Y) <- samp.name
  rownames(Z) <- samp.name
  wald[["bias"]] <- c(bias.intc, bias)
  return(list(variables = variables, bias = bias, output = output, 
              covariance = covariance, otu.tab.use = Y, meta.use = Z, 
              wald = wald))
}