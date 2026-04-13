
ancombc1.fun <- function(m,otu.tab, meta, formula, alpha,p_adj_method) {
  m=nrow(otu.tab)
  otu_mat=as.matrix(otu.tab)
  rownames(otu_mat) = paste0("taxon", 1:nrow(otu_mat)) 
  colnames(otu_mat) = paste0("sample", 1:ncol(otu_mat))
  assays = SimpleList(counts = otu_mat)
  
  
  smd = data.frame(meta,
                   row.names = paste0("sample", 1:ncol(otu_mat)),
                   stringsAsFactors = FALSE)
  smd = DataFrame(smd)
  colnames(smd)<-colnames(meta)
  
  aaaaa<-c(1:m)
  tax_tab = matrix(aaaaa, ncol = 1)
  rownames(tax_tab) = rownames(otu_mat)
  tax_tab = DataFrame(tax_tab)
  tse = TreeSummarizedExperiment(assays = assays,
                                 colData = smd,
                                 rowData = tax_tab)
  out = ancombc(data = tse, assay_name = "counts",
                tax_level = NULL, phyloseq = NULL,
                formula = formula,
                p_adj_method = p_adj_method, prv_cut = 0, lib_cut = 0,
                group = "u", struc_zero = FALSE, neg_lb = FALSE,
                tol = 1e-5, max_iter = 100, conserve = TRUE,
                alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
  
  res_prim = out$res
  rej<-which(res_prim$diff_abn[,3]==TRUE)
  
  indind<-match(res_prim$diff_abn[rej,1],rownames(otu_mat))
  return(indind)
}


ancombc2.fun <- function(m,otu.tab, meta, formula, alpha,p_adj_method) {
  meta<-as.data.frame(meta)
  otu_mat=as.matrix(otu.tab)
  rownames(otu_mat)=rownames(otu.tab)
  colnames(otu_mat) = colnames(otu.tab)
  assays = SimpleList(counts = otu_mat)
  smd = data.frame(meta,
                   row.names = colnames(otu.tab),
                   stringsAsFactors = FALSE)
  smd = DataFrame(smd)
  print(colnames(meta))
  colnames(smd)<-colnames(meta)
  aaaaa<-c(1:m)
  tax_tab = matrix(aaaaa, ncol = 1)
  rownames(tax_tab) = rownames(otu_mat)
  tax_tab = DataFrame(tax_tab)
  assays = S4Vectors::SimpleList(counts = otu_mat)
  smd = S4Vectors::DataFrame(smd)
  tse = TreeSummarizedExperiment(assays = assays,
                                 colData = smd)
  out = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                 fix_formula = formula,
                 rand_formula = NULL,
                 p_adj_method = p_adj_method, pseudo_sens = FALSE,
                 prv_cut = 0, lib_cut = 1, s0_perc = 0.05,
                 group = "u", struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 1),
                 trend_control = list(contrast =
                                        list(matrix(c(1, 0, -1, 1),
                                                    nrow = 2,
                                                    byrow = TRUE)),
                                      node = list(2),
                                      solver = "ECOS",
                                      B = 1))
  rej<-which(out$res$diff_u1=="TRUE")
  indind<-match(out$res$taxon[rej],rownames(otu_mat))
  return(indind)
}


ancombc2.fun.ss <- function(m,otu.tab, meta, formula, alpha,p_adj_method) {
  meta<-as.data.frame(meta)
  otu_mat=as.matrix(otu.tab)
  rownames(otu_mat)=rownames(otu.tab)
  colnames(otu_mat) = colnames(otu.tab)
  assays = SimpleList(counts = otu_mat)
  smd = data.frame(meta,
                   row.names = colnames(otu.tab),
                   stringsAsFactors = FALSE)
  smd = DataFrame(smd)
  print(colnames(meta))
  colnames(smd)<-colnames(meta)
  aaaaa<-c(1:m)
  tax_tab = matrix(aaaaa, ncol = 1)
  rownames(tax_tab) = rownames(otu_mat)
  tax_tab = DataFrame(tax_tab)
  assays = S4Vectors::SimpleList(counts = otu_mat)
  smd = S4Vectors::DataFrame(smd)
  tse = TreeSummarizedExperiment(assays = assays,
                                 colData = smd)
  out = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                 fix_formula = formula,
                 rand_formula = NULL,
                 p_adj_method = p_adj_method, pseudo_sens = TRUE,
                 prv_cut = 0, lib_cut = 1, s0_perc = 0.05,
                 group = "u", struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 1),
                 trend_control = list(contrast =
                                        list(matrix(c(1, 0, -1, 1),
                                                    nrow = 2,
                                                    byrow = TRUE)),
                                      node = list(2),
                                      solver = "ECOS",
                                      B = 1))
  rej<-which(out$res$diff_u1=="TRUE")
  indind<-match(out$res$taxon[rej],rownames(otu_mat))
  return(indind)
}




ancombc2.fun.cont <- function(m,otu.tab, meta, formula, alpha,p_adj_method) {
  meta<-as.data.frame(meta)
  otu_mat=as.matrix(otu.tab)
  rownames(otu_mat)=rownames(otu.tab)
  colnames(otu_mat) = colnames(otu.tab)
  assays = SimpleList(counts = otu_mat)

  
  smd = data.frame(meta,
                   row.names = colnames(otu.tab),
                   stringsAsFactors = FALSE)
  smd = DataFrame(smd)
  print(colnames(meta))
  colnames(smd)<-colnames(meta)
  
  aaaaa<-c(1:m)
  tax_tab = matrix(aaaaa, ncol = 1)
  rownames(tax_tab) = rownames(otu_mat)
  tax_tab = DataFrame(tax_tab)
  
  assays = S4Vectors::SimpleList(counts = otu_mat)
  smd = S4Vectors::DataFrame(smd)
  tse = TreeSummarizedExperiment(assays = assays,
                                 colData = smd)

  out = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                 fix_formula = formula,
                 rand_formula = NULL,
                 p_adj_method = p_adj_method, pseudo_sens = FALSE,
                 prv_cut = 0, lib_cut = 1, s0_perc = 0.05,
                 group = NULL, struc_zero = FALSE, neg_lb = FALSE, #############
                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 1),
                 trend_control = list(contrast =
                                        list(matrix(c(1, 0, -1, 1),
                                                    nrow = 2,
                                                    byrow = TRUE)),
                                      node = list(2),
                                      solver = "ECOS",
                                      B = 1))
  
  rej<-which(out$res$diff_u=="TRUE")
  indind<-match(out$res$taxon[rej],rownames(otu_mat))
  return(indind)
}



ancombc2.fun.contss <- function(m,otu.tab, meta, formula, alpha,p_adj_method) {
  meta<-as.data.frame(meta)
  otu_mat=as.matrix(otu.tab)
  rownames(otu_mat)=rownames(otu.tab)
  colnames(otu_mat) = colnames(otu.tab)
  assays = SimpleList(counts = otu_mat)
  
  
  smd = data.frame(meta,
                   row.names = colnames(otu.tab),
                   stringsAsFactors = FALSE)
  smd = DataFrame(smd)
  print(colnames(meta))
  colnames(smd)<-colnames(meta)
  
  aaaaa<-c(1:m)
  tax_tab = matrix(aaaaa, ncol = 1)
  rownames(tax_tab) = rownames(otu_mat)
  tax_tab = DataFrame(tax_tab)
  
  assays = S4Vectors::SimpleList(counts = otu_mat)
  smd = S4Vectors::DataFrame(smd)
  tse = TreeSummarizedExperiment(assays = assays,
                                 colData = smd)
  
  out = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                 fix_formula = formula,
                 rand_formula = NULL,
                 p_adj_method = p_adj_method, pseudo_sens = TRUE,
                 prv_cut = 0, lib_cut = 1, s0_perc = 0.05,
                 group = NULL, struc_zero = FALSE, neg_lb = FALSE, #############
                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 1),
                 trend_control = list(contrast =
                                        list(matrix(c(1, 0, -1, 1),
                                                    nrow = 2,
                                                    byrow = TRUE)),
                                      node = list(2),
                                      solver = "ECOS",
                                      B = 1))
  
  rej<-which(out$res$diff_u=="TRUE")
  indind<-match(out$res$taxon[rej],rownames(otu_mat))
  return(indind)
}