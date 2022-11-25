## This file implements the differential expression functions to be passed to
## the SummarizedBenchmark object

#' ZINBwaVe
compute_weights <- function(countdata, coldat, ...) {

  ## alternative code to use Zinbwave?
  # https://github.com/mikelove/zinbwave-deseq2/blob/master/zinbwave-deseq2.Rmd
  # zinb <- sim[keep,]
  # zinb$condition <- factor(zinb$Group)
  # # we need to reorganize the assays in the SumExp from splatter
  # nms <- c("counts", setdiff(assayNames(zinb), "counts"))
  # assays(zinb) <- assays(zinb)[nms]
  # assay(zinb) <- as.matrix(assay(zinb))
  # zinb <- zinbwave(zinb, K=0, observationalWeights=TRUE,
  #                  BPPARAM=SerialParam(), epsilon=1e12)
  # dds <- DESeqDataSet(zinb, design=~condition)
  set.seed(12345)
  computeExactWeights <- function (model, x) {
    set.seed(12345)
    library(zinbwave)

    mu <- getMu(model)
    pi <- getPi(model)
    theta <- getTheta(model)
    theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
    nb_part <- dnbinom(t(x), size = theta, mu = mu)
    zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
    zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
    zinbwg <- t(zinbwg)
    zinbwg[x > 0] <- 1
    zinbwg[zinbwg < 1e-15] <- 1e-15
    zinbwg
  }

  library(zinbwave)

  zinbmodel <- zinbFit(Y = countdata,
                       X = model.matrix(~ coldat),
                       K = 0,
                       # epsilon = nrow(countdata), #1e10
                       commondispersion = TRUE,
                       verbose = FALSE,
                       ...) #BPPARAM = BiocParallel::SerialParam()

  computeExactWeights(model = zinbmodel,
                      x = countdata)
}

#' DESeq2 defaults
deseq2_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  # message("Start DESeq2 default")
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds,
                         contrast = contrast,
                         test = "Wald",
                         independentFiltering = F)
  # or to shrink log fold changes association with condition:
  # lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

  runtime <- tictoc::toc(log = F, quiet = T)

  # message("End DESeq2 default")

  list(res = res,
       runtime = runtime$toc - runtime$tic)

}

#' 2. use LRT
deseq2lrt_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)
  dds <- DESeq2::DESeq(dds,
                       test = "LRT",
                       reduced = ~ 1)

  res <- DESeq2::results(dds, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 3. use betaPrior
deseq2bp_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)
  dds <- DESeq2::DESeq(dds,
                       betaPrior = T)
  res <- DESeq2::results(dds,
                         contrast = contrast,
                         test = "Wald",
                         independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 4. use recommended parameters for single-cell data
## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
## Use test="LRT"
## useT=TRUE
## minmu=1e-6 The default setting of minmu was benchmarked on bulk RNA-seq and
##            is not appropriate for single cell data when the expected count is
##            often much less than 1
## minReplicatesForReplace=Inf.
## setting sizeFactors from scran::computeSumFactors. The default size factors
##         are not optimal for single cell count matrices
## set fitType = "glmGamPoi"
deseq2zi_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)

  dds <- DESeq2::DESeq(dds,
                       # quiet = TRUE,
                       sfType = "poscounts",
                       useT = TRUE,
                       minmu = 1e-6,
                       minReplicatesForReplace = Inf,
                       test = "LRT",
                       reduced = ~ 1)
  # DESeq2::results(dds, contrast = contrast)
  res <- DESeq2::results(dds, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 5. use recommended parameters for single-cell data and the genefilter::shorth
## in estimating  the size factors (parameter 'locfunc')
## ?estimateSizeFactors
## locfunc: a function to compute a location for a sample. By default, the
##          median is used. However, especially for low counts, the shorth
##          function from the genefilter package may give better results.
deseq2lc_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)

  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts",
                                     locfunc = function(x){
                                       tryCatch(
                                         {
                                           genefilter::shorth(x)
                                         },
                                         error = function(cond) {
                                           warning(paste("The 'shorth' function failed with the following message:\n",
                                                         cond,
                                                         "Will try to use 'half.range.mode' instead"))
                                           genefilter::half.range.mode(x)}
                                       ) })

  dds <- DESeq2::DESeq(dds,
                       # quiet = TRUE,
                       # sfType = "poscounts",
                       useT = TRUE,
                       minmu = 1e-6,
                       minReplicatesForReplace = Inf,
                       test = "LRT",
                       reduced = ~ 1)
  # results(dds, contrast = contrast)
  res <- DESeq2::results(dds, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 6. use recommended parameters for single-cell data and the scran::computeSumFactors
## to estimate the size factors
## 'setting sizeFactors from scran::computeSumFactors. The default size factors
##  are not optimal for single cell count matrices'
deseq2sc_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)

  # "sizeFactors"(dds) <-
  #   DESeq2::sizeFactors(scran::computeSumFactors(dds,
  #                                                clusters = SummarizedExperiment::colData(dds)$condition))
  dds <- scran::computeSumFactors(dds,
                                  clusters = SummarizedExperiment::colData(dds)$condition)

  dds <- DESeq2::DESeq(dds,
                       # quiet = TRUE,
                       useT = TRUE,
                       minmu = 1e-6,
                       minReplicatesForReplace = Inf,
                       test = "LRT",
                       reduced = ~ 1)
  # DESeq2::results(dds, contrast = contrast)
  res <- DESeq2::results(dds, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 7. By using the argument fitType="glmGamPoi", one can leverage the faster NB GLM
## engine written by Constantin Ahlmann-Eltze. Note that glmGamPoi’s interface in
## DESeq2 requires use of test="LRT" and specification of a reduced design.
## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#speed-up-and-parallelization-thoughts
deseq2gpLrt_run <- function(countData, colData, design, contrast) {
  set.seed(12345)
  library(glmGamPoi)

  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)

  dds <- DESeq2::DESeq(dds,
                       # quiet = TRUE,
                       sfType = "poscounts",
                       # useT = TRUE,
                       minmu = 1e-6,
                       minReplicatesForReplace = Inf,
                       fitType = "glmGamPoi",
                       test = "LRT",
                       reduced = ~ 1)
  # DESeq2::results(dds, contrast = contrast)
  res <- DESeq2::results(dds, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 8. use ZinbWave weights
deseq2zw_run <- function(countData, colData, design, contrast, weights) {
  set.seed(12345)
  tictoc::tic()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData = colData, design = design)

  # if(is.null(weights)){
  #   message("Weights expected but NULL was given. DESeq2 will not use weights")
  # }else{
  weights[which(weights < 1e-6)] <- 1e-06
  SummarizedExperiment::assays(dds, withDimnames = F, "weights")[["weights"]] <- weights
  # }

  dds <- DESeq2::DESeq(dds,
                       # quiet = TRUE,
                       sfType = "poscounts",
                       # useT = TRUE,
                       minmu = 1e-6,
                       minReplicatesForReplace = Inf,
                       test = "LRT",
                       reduced = ~ 1)
  # DESeq2::results(dds, contrast = contrast)
  res <- DESeq2::results(dds, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## post functions
deseq2_pv <- function(x) {
  res <- x$res$pvalue
  # res[is.na(res)] <- 1
  res
}

deseq2_apv <- function(x) {
  res <- x$res$padj
  # res[is.na(res)] <- 1
  res
}

deseq2_lfc <- function(x) { x$res$log2FoldChange }

deseq2_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### edgeR family


## edgeR family

## 1. default edgeR
edgeR_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y) # method = "TMMwsp" for zero counts?
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- edgeR::estimateDisp(y, des) # min.row.sum = 5 default. Change?
  fit <- edgeR::glmFit(y, des)
  res <- edgeR::glmLRT(fit, coef=2)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 2. robust
edgeRrbst_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y) # method = "TMMwsp" for zero counts?
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- edgeR::estimateGLMRobustDisp(y, des) #maxit = 6
  fit <- edgeR::glmQLFit(y = y,
                         dispersion = y$tagwise.dispersion,
                         robust = TRUE,
                         design = des)
  res <- edgeR::glmLRT(fit, coef = 2)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 3. 50 degrees of freedom
## https://support.bioconductor.org/p/84338/
edgeRrbst50df_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y) # method = "TMMwsp" for zero counts?
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- edgeR::estimateGLMRobustDisp(y = y, des, prior.df = 50) #maxit = 6
  fit <- edgeR::glmQLFit(y = y,
                         dispersion = y$tagwise.dispersion,
                         robust = TRUE,
                         design = des)
  res <- edgeR::glmLRT(fit, coef = 2)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 4. auto estimation of df
edgeRrbstEdf_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y) # method = "TMMwsp" for zero counts?
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- edgeR::estimateDisp(y, des) #maxit = 6
  fit <- edgeR::glmQLFit(y = y,
                         dispersion = y$tagwise.dispersion,
                         robust = TRUE,
                         design = des)
  res <- edgeR::glmLRT(fit, coef = 2)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## 5. ZinbWave weights
edgeRzw_run <- function(countData, group, design, weights) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y)
  y$weights <- weights
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  # y <- edgeR::estimateGLMRobustDisp(y, des)
  y <- edgeR::estimateDisp(y, des)
  fit <- edgeR::glmFit(y, des)
  # res <- edgeR::glmLRT(fit, coef = 2)
  res <- zinbwave::glmWeightedF(fit, coef = 2, independentFiltering = F)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)

}

## 6. Quasi-likelihood dispersion estimate and empirical Bayes quasi-likelihood F-tests
## glmQLFit gives special attention to handling of zero counts
## https://support.bioconductor.org/p/84338/
## and
## https://support.bioconductor.org/p/84291/#84292
## "There are some situations where the QL F-test doesn't work well - for example, [...]
## where the dispersions are very large and the counts are very small, whereby
## some of the approximations in the QL framework seem to fail. In such cases,
## I usually switch to the LRT rather than using the exact test, for the reasons
## of experimental flexibility that I mentioned above."
edgeRql_run <- function(countData, group, design, weights) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y)
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- edgeR::estimateDisp(y = y, design = des, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design = des, robust = TRUE)
  res <- edgeR::glmQLFTest(fit, coef = 2)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

#### edgeR TMMwsp robust
edgeRrbstTMMswp_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y, method = "TMMwsp") ## "TMMwsp" for zero counts
  des <- model.matrix(~group, data = y$samples)
  y <- edgeR::estimateGLMRobustDisp(y, des)
  fit <- edgeR::glmQLFit(y = y,
                         dispersion = y$tagwise.dispersion,
                         robust = TRUE,
                         design = des)
  res <- edgeR::glmLRT(fit, coef = 2)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## post functions
edgeR_pv <- function(x) {
  res <- x$res$table$PValue
  # res[is.na(res)] <- 1
  res
}

edgeR_apv <- function(x) {
  res <- p.adjust(p = x$res$table$PValue, method = "BH")
  # res <- topTags(x, number = Inf, sort.by = "none")$table$FDR
  # res[is.na(res)] <- 1
  res
}

edgeR_lfc <- function(x) { x$res$table$logFC }

edgeR_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


# ## post functions
# edgeR_pv <- function(x) {
#   res <- x$res$table$PValue
#   res[is.na(res)] <- 1
#   res
# }
#
# edgeR_apv <- function(x) {
#   res <- p.adjust(p = x$res$table$PValue, method = "BH")
#   # res <- topTags(x, number = Inf, sort.by = "none")$table$FDR
#   res[is.na(res)] <- 1
#   res
# }
#
# edgeR_lfc <- function(x) { x$res$table$logFC }
#
# edgeR_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### limma-voom family


## limma-voom family
voom_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y)
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- limma::voom(y, des)
  res <- limma::eBayes(limma::lmFit(y, des))

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## estimation of df.prior and var.prior be robustified against outlier sample variances
voomRbst_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  y <- edgeR::DGEList(countData, group = group)
  y <- edgeR::calcNormFactors(y) # method = "TMMwsp" for zero counts?
  # des <- model.matrix(design)
  des <- model.matrix(~group, data = y$samples)
  y <- limma::voom(y, des)
  res <- limma::eBayes(limma::lmFit(y, des),
                       robust = T)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## voom quantile normalization
voomQn_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  # des <- model.matrix(design)
  y <- edgeR::DGEList(countData, group = group)
  des <- model.matrix(~group, data = y$samples)
  voom.data <- limma::voom(countData,
                           design = des,
                           normalize.method = 'quantile')
  #"none", "scale", "quantile" or "cyclicloess
  res <- limma::eBayes(limma::lmFit(voom.data, design = des))

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## limma-voom code as in the vignette (no prior edgeR; mind the robust = T param)
voomSimple_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  # des <- model.matrix(design)
  y <- edgeR::DGEList(countData, group = group)
  des <- model.matrix(~group, data = y$samples)
  res <- limma::eBayes(limma::lmFit(limma::voom(counts = countData,
                                                design = des),
                                    des),
                       robust = TRUE)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}


## voomLmFit is more robust to zero counts than calling voom,
## duplicateCorrelation and lmFit separately and provides more
## rigorous error rate control.
## ?limma::voom
## Note that edgeR::voomLmFit is now recommended over voom for
## sparse counts with a medium to high proportion of zeros.
voomLmFit_run <- function(countData, group, design) {
  set.seed(12345)
  tictoc::tic()

  # des <- model.matrix(design)
  y <- edgeR::DGEList(countData, group = group)
  des <- model.matrix(~group, data = y$samples)
  ## Empirical sample quality weights will be estimated if sample.weights=TRUE
  ## or if var.design or var.group are non-NULL. In that case, voomLmFit is
  ## analogous to running voomWithQualityWeights followed by lmFit.
  ## voomLmFit is usually followed by running eBayes on the fitted model object.
  res <- limma::eBayes(edgeR::voomLmFit(counts = countData,
                                        design = des,
                                        # block = NULL,
                                        # prior.weights = NULL,
                                        sample.weights = TRUE,
                                        # var.design = NULL,
                                        # var.group = NULL,
                                        # lib.size = NULL,
                                        # normalize.method = "none",
                                        # span = 0.5,
                                        # plot = FALSE,
                                        save.plot = FALSE))

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

# limma-voom with ZINBWaVE weights
voomZw_run <- function(countData, group, design, weights) {
  set.seed(12345)
  tictoc::tic()

  # des <- model.matrix(design)
  y <- edgeR::DGEList(countData, group = group)
  des <- model.matrix(~group, data = y$samples)
  v <- limma::voom(counts = countData,
                   design = des,
                   plot = FALSE,
                   weights = weights)
  # v$weights <- v$weights * weights
  fit <- limma::lmFit(v, design = des, weights = v$weights)
  # fit$df.residual <- rowSums(weights) - ncol(design)
  res <- limma::eBayes(fit, robust = T)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

# voomWithQualityWeights

## post functions
voom_pv <- function(x) {
  res <- x$res$p.value[, 2] # res <- x$F.p.value	?
  # res[is.na(res)] <- 1
  res
}

voom_apv <- function(x) {
  # res <- p.adjust(p = x$p.value[, 2], method = "BH")
  res <- limma::topTable(x$res, number = Inf, sort.by = "none")$adj.P.Val
  # res[is.na(res)] <- 1
  res
}

voom_lfc <- function(x) {
  # x$coefficients[, 2]
  limma::topTable(x$res, number = Inf, sort.by = "none")$logFC
}

voom_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### circMeta


## circMeta
cMeta_run <- function(countData, group, sf = TRUE) {
  set.seed(12345)
  tictoc::tic()

  groupVar <- group
  names(groupVar) <- colnames(countData)
  countData <- apply(countData, 2, function(x) {storage.mode(x) <- 'integer'; x})
  m = rowMeans(countData)

  sfs = colSums(countData)
  sfs = sfs / min(sfs)
  if(sf) countData = sweep(countData, 2, sfs, FUN='/')
  n0 = sum(groupVar == levels(groupVar)[1])
  n1 = sum(groupVar == levels(groupVar)[2]) # assume 2 levels factor groups
  m0 = rowMeans(countData[, groupVar == levels(groupVar)[1]])
  m1 = rowMeans(countData[, groupVar == levels(groupVar)[2]])
  n  = nrow(countData)
  pval = rep(1, n)
  for(i in 1:n){
    z = (m1[i] - m0[i]) / sqrt(m1[i] / n1 + m0[i] / n0)
    pval[i] = 2 * pnorm(-abs(z))
  }

  # fdr = p.adjust(pval, method = 'fdr')
  lfc <- log( (m1 + 1) / (m0 + 1), 2)

  names(pval) <- rownames(countData)
  res <- data.frame(p.value = pval,
                    logFC = lfc)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

cMetaLc_run <- function(countData, group, sf = TRUE) {
  set.seed(12345)
  tictoc::tic()

  groupVar <- group
  names(groupVar) <- colnames(countData)
  countData <- apply(countData, 2, function(x) {storage.mode(x) <- 'integer'; x})
  m = rowMeans(countData)

  sfs = colSums(countData)
  sfs = sfs / min(sfs)
  if(sf) countData = sweep(countData, 2, sfs, FUN='/')
  n0 = sum(groupVar == levels(groupVar)[1])
  n1 = sum(groupVar == levels(groupVar)[2]) # assume 2 levels factor groups
  m0 = rowMeans(countData[, groupVar == levels(groupVar)[1]])
  m1 = rowMeans(countData[, groupVar == levels(groupVar)[2]])
  n  = nrow(countData)
  pval = rep(1, n)
  for(i in 1:n){
    z = ( sqrt(m1[i]) - sqrt(m0[i]) ) / (1 / 2 * sqrt(1 / n1 + 1 / n0))
    pval[i] = 2 * pnorm(-abs(z))
  }

  # fdr = p.adjust(pval, method = 'fdr')
  lfc <- log( (m1 + 1) / (m0 + 1), 2)

  names(pval) <- rownames(countData)
  res <- data.frame(p.value = pval,
                    logFC = lfc)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## post functions
cMeta_pv <- function(x) {
  res <- x$res$p.value
  # res[is.na(res)] <- 1
  res
}

cMeta_apv <- function(x) {
  res <- p.adjust(x$res$p.value, method = "BH")
  # res[is.na(res)] <- 1
  res
}

cMeta_lfc <- function(x) { x$res$logFC }

cMeta_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### lncDIFF


lncdiff_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  dgeList <- edgeR::DGEList(counts = countData)
  # calculate normalization factors
  dgeList <- edgeR::calcNormFactors(object = dgeList, method = "TMM")
  ## normalize counts
  norm_data <- edgeR::cpm(dgeList,
                          normalized.lib.sizes = TRUE,
                          log = F)

  res <- lncDIFF::lncDIFF(edata = norm_data,
                          group = as.character(group),
                          covariate = NULL,
                          link.function = 'log',
                          CompareGroups = unique(as.character(group)),
                          simulated.pvalue = FALSE,
                          permutation = 100)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)

}

lncdiff_pv <- function(x) {
  res <- x$res$DE.results$Pvalue
  # res[is.na(res)] <- 1
  res
}

lncdiff_apv <- function(x) {
  res <- p.adjust(x$res$DE.results$Pvalue, method = "BH")
  # res[is.na(res)] <- 1
  res
}

lncdiff_lfc <- function(x) {
  x$res$DE.results$Log2.Fold.Change
}

lncdiff_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$DE.results$Pvalue)) }


### SAMseq


samseq_run <- function(countData, group) {
  set.seed(12345)
  library(samr)

  tictoc::tic()

  condition.12 <- rep(1, length(group))
  condition.12[which(group == levels(factor(group))[2])] <- 2

  SAMseq.test <- samr::SAMseq(countData,
                              condition.12,
                              resp.type = 'Two class unpaired',
                              geneid = rownames(countData),
                              genenames = rownames(countData),
                              nperms = 100,
                              nresamp = 20,
                              fdr.output = 1)
  pv <- samr.pvalues.from.perms(tt = SAMseq.test$samr.obj$tt,
                                ttstar = SAMseq.test$samr.obj$ttstar)

  SAMseq.result <- rbind(SAMseq.test$siggenes.table$genes.up,
                         SAMseq.test$siggenes.table$genes.lo)

  SAMseq.statistic <- rep(0, nrow(countData))

  SAMseq.statistic[match(SAMseq.result[, 1],
                         rownames(countData))] <- as.numeric(SAMseq.result[, 3])

  SAMseq.FDR <- rep(1, nrow(countData))

  SAMseq.FDR[match(SAMseq.result[, 1],
                   rownames(countData))] <- as.numeric(SAMseq.result[, 5]) / 100

  SAMseq.FC <- rep(NA, nrow(countData))

  SAMseq.FC[match(SAMseq.result[, 1],
                  rownames(countData))] <- as.numeric(SAMseq.result[, 4])

  SAMseq.score <- 1 - SAMseq.FDR

  res <- data.frame('statistic' = SAMseq.statistic,
                    'Pvalue' = pv,
                    'FDR' = SAMseq.FDR,
                    'FC' = SAMseq.FC,
                    'score' = SAMseq.score)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

samseq_pv <- function(x) {
  res <- x$res$Pvalue
  # res[is.na(res)] <- 1
  res
}


samseq_apv <- function(x) {
  # p.adjust(x$res$Pvalue, "BH")
  res <- x$res$FDR
  # res[is.na(res)] <- 1
  res
}

samseq_lfc <- function(x) {
  log2(x$res$FC)
}

samseq_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### NBID


nbid_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  res <- NBID::DEUsingNBID(countData, group, ncore = 1)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

## use NBID with scran size factors
nbid_sc_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  # sf <- scran::calculateSumFactors(countData, clusters = group)
  sf <- scuttle::pooledSizeFactors(countData, clusters = group)
  # sf <- scuttle::pooledSizeFactors(countData) # this reduces to librarySizeFactors(countData)
  res <- NBID::DEUsingNBID(countData, group, ncore = 1, countPerCell = sf)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

nbid_pv <- function(x) {
  res <- x$res[, "pvalue"]
  # res[is.na(res)] <- 1
  res
}

nbid_apv <- function(x) {
  res <- p.adjust(x$res[, "pvalue"], method = "BH")
  # res[is.na(res)] <- 1
  res
}

nbid_lfc <- function(x) {
  x$res[, "log2FC2"]
}

nbid_time <- function(x) {
  rep(as.numeric(x$runtime), length(x$res[, "pvalue"]))
}


### NOISeqBIO


noiseqbio_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  res <- NOISeq::noiseqbio(NOISeq::readData(countData, group),
                           k = 0.5,
                           norm = "tmm",
                           factor = "condition",
                           lc = 0,
                           adj = 1.5,
                           plot = FALSE,
                           a0per = 0.9,
                           random.seed = 12345,
                           filter = 0) #disabled filter

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

noiseqbio_pv <- function(x) {

  ## NB: NOISeq does not provide P-values, but only "adjusted p-value"-like probabilities
  res <- 1 - x$res@results[[1]]$prob
  # res[is.na(res)] <- 1
  res
}

noiseqbio_apv <- function(x) {

  # res@results[[1]]$prob gives the estimated probability of differential
  # expression for each feature. Note that when using NOISeq, these probabilities
  # are not equivalent to p-values. The higher the probability, the more likely
  # that the difference in expression is due to the change in the experimental
  # condition and not to chance.
  # when using NOISeqBIO, the probability of differential expression would be
  # equivalent to 1 − FDR, where FDR can be considered as an adjusted p-value.

  res <- 1 - x$res@results[[1]]$prob
  # res[is.na(res)] <- 1
  res
}

noiseqbio_lfc <- function(x) {
  x$res@results[[1]]$log2FC
}

noiseqbio_time <- function(x) {
  rep(as.numeric(x$runtime), length(x$res@results[[1]]$log2FC))
}


### PoissonSeq


poissonSeq_run <- function(countData, group) {
  set.seed(12345)
  # library(PoissonSeq)
  library(splines)
  tictoc::tic()

  dat <- list(n = countData,
              y = as.numeric(group),
              pair = FALSE,
              type = 'twoclass',
              gname = rownames(countData))

  ## prevent concurrent parallel access to the pow.txt file
  pow_filename <- tempfile(pattern = "poissonSeq_pow.txt", tmpdir = tempdir())

  res <- PoissonSeq::PS.Main(dat = dat,
                             para = list(ct.sum = 2,
                                         ct.mean = 0.1,
                                         pow.file = pow_filename))

  res <- res[rownames(countData), ]

  runtime <- tictoc::toc(log = F, quiet = T)

  if (file.exists(pow_filename)) {
    #Delete file if it exists
    file.remove(pow_filename)
  }

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

poissonSeq_pv <- function(x) {
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

poissonSeq_apv <- function(x) {
  res <- x$res$fdr
  # res[is.na(res)] <- 1
  res
}

poissonSeq_lfc <- function(x) {
  x$res$log.fc
}

poissonSeq_time <- function(x) {
  rep(as.numeric(x$runtime), length(x$res$pval))
}


### Limma + sctransform::vst


limmaStVst_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  sampleTab <- data.frame(condition = setNames(group, colnames(countData)))
  des <- model.matrix(~ condition, data = sampleTab)

  set.seed(12345)
  vst_out <- sctransform::vst(countData,
                              cell_attr = sampleTab,
                              latent_var = c("condition"), #c("log_umi"),
                              method = "glmGamPoi",
                              # method = "glmGamPoi_offset",
                              # method = "offset",
                              return_gene_attr = TRUE,
                              return_cell_attr = TRUE,
                              min_cells = 3,
                              n_genes = NULL,
                              verbosity = 1)

  fit <- limma::eBayes(limma::lmFit(vst_out$y, des))

  res <- limma::topTable(fit,
                         sort.by = "none",   # number = nrow(countData)
                         number = Inf)[rownames(countData), ]

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

limmaStVst_pv <- function(x) {
  res <- x$res$P.Value
  # res[is.na(res)] <- 1
  res
}

limmaStVst_apv <- function(x) {
  # res <- p.adjust(p = x$p.value[, 2], method = "BH")
  res <- x$res$adj.P.Val
  # res[is.na(res)] <- 1
  res
}

limmaStVst_lfc <- function(x) {
  x$res$logFC
}

limmaStVst_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### ROTS CPM


# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_ROTScpm.R
ROTScpm_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  grp <- as.character(group)
  dge <- edgeR::DGEList(counts = countData, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  cpms <- edgeR::cpm(dge)
  rots <- ROTS::ROTS(data = cpms,
                     groups = grp,
                     B = 1000,
                     K = 1000,
                     log = FALSE,
                     seed = 123)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = rots,
       runtime = runtime$toc - runtime$tic)
}

ROTScpm_pv <- function(x) {
  res <- x$res$pvalue
  # res[is.na(res)] <- 1
  res
}

ROTScpm_apv <- function(x) {
  # res <- p.adjust(p = x$res$pvalue, method = "BH")
  res <- x$res$FDR
  # res[is.na(res)] <- 1
  res
}

ROTScpm_lfc <- function(x) {
  x$res$logfc
}

ROTScpm_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res$data)) }


### DESingle


DEsingle_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  results <- DEsingle::DEsingle(counts = countData, group = group)
  results <- data.frame(results, stringsAsFactors = FALSE)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = results,
       runtime = runtime$toc - runtime$tic)

}

DEsingle_pv <- function(x) {
  res <- x$res$pvalue
  # res[is.na(res)] <- 1
  res
}

DEsingle_apv <- function(x) {
  # res <- p.adjust(p = x$res$pvalue, method = "BH")
  res <- x$res$pvalue.adj.FDR
  # res[is.na(res)] <- 1
  res
}

DEsingle_lfc <- function(x) {
  log2(x$res$foldChange)
}

DEsingle_time <- function(x) { rep(as.numeric(x$runtime), nrow(x$res)) }


### metagenomeSeq


# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_metagenomeSeq.R
metagenomeSeq_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  sampleTable <- data.frame(condition = group,
                            row.names = colnames(countData))

  obj <- metagenomeSeq::newMRexperiment(countData,
                                        phenoData = new("AnnotatedDataFrame",
                                                        data = sampleTable))
  p <- metagenomeSeq::cumNormStatFast(obj)
  obj <- metagenomeSeq::cumNorm(obj, p = p)
  mod <- model.matrix(~ condition, data = sampleTable) # data = Biobase::pData(obj)
  res <- metagenomeSeq::fitFeatureModel(obj = obj, mod = mod, coef = 2)
  tbl <- metagenomeSeq::MRtable(obj = res, number = Inf, by = 2,
                                adjustMethod = "BH")[rownames(countData), ]

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = tbl,
       runtime = runtime$toc - runtime$tic)
}

metagenomeSeq_pv <- function(x) {
  res <- x$res$pvalues
  # res[is.na(res)] <- 1
  res
}

metagenomeSeq_apv <- function(x) {
  res <- x$res$adjPvalues
  # res[is.na(res)] <- 1
  res
}

metagenomeSeq_lfc <- function(x) {
  x$res$logFC
}

metagenomeSeq_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$pvalues)) }


### Seurat family


## Seurat Bimod
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_SeuratTobit.R
seuratBimod_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  colnames(countData) <- paste0(group, "__", 1:ncol(countData))
  ##CreateSeuratObject::CreateSeuratObject()
  seur <- Seurat::CreateSeuratObject(counts = countData,
                                     project = "scrnaseq",
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = "__")
  res <- Seurat::FindMarkers(seur,
                             # min.pct = 0,
                             verbose = F,
                             ident.1 = levels(factor(group))[1],
                             ident.2 = levels(factor(group))[2],
                             test.use = "bimod") #"wilcox" "MAST"

  ## FindMarkers may return less features than the input
  ## the following fixes the missing features
  missed_genes <- rownames(countData)[which(!sub("_", "-", rownames(countData)) %in% rownames(res))]

  nadf <- data.frame("p_val" = rep(NA, length(missed_genes)),
                     "avg_log2FC" = rep(NA, length(missed_genes)),
                     "pct.1" = rep(NA, length(missed_genes)),
                     "pct.2" = rep(NA, length(missed_genes)),
                     "p_val_adj" = rep(NA, length(missed_genes)),
                     row.names = sub("_", "-", missed_genes))

  ## attach missed rows and order according to initial names
  res <- rbind(res, nadf)[sub("_", "-", rownames(countData)), ]
  rownames(res) <- rownames(countData)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}



### Seurat Wilcoxon
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_SeuratTobit.R
seuratWilcox_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  colnames(countData) <- paste0(group, "__", 1:ncol(countData))
  ##CreateSeuratObject::CreateSeuratObject()
  seur <- Seurat::CreateSeuratObject(counts = countData,
                                     project = "scrnaseq",
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = "__")
  res <- Seurat::FindMarkers(seur,
                             # min.pct = 0,
                             verbose = F,
                             ident.1 = levels(factor(group))[1],
                             ident.2 = levels(factor(group))[2],
                             test.use = "wilcox")

  ## FindMarkers may return less features than the input
  ## the following fixes the missing features
  missed_genes <- rownames(countData)[which(!sub("_", "-", rownames(countData)) %in% rownames(res))]

  nadf <- data.frame("p_val" = rep(NA, length(missed_genes)),
                     "avg_log2FC" = rep(NA, length(missed_genes)),
                     "pct.1" = rep(NA, length(missed_genes)),
                     "pct.2" = rep(NA, length(missed_genes)),
                     "p_val_adj" = rep(NA, length(missed_genes)),
                     row.names = sub("_", "-", missed_genes))

  ## attach missed rows and order according to initial names
  res <- rbind(res, nadf)[sub("_", "-", rownames(countData)), ]
  rownames(res) <- rownames(countData)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       runtime = runtime$toc - runtime$tic)
}

seurat_pv <- function(x) {
  res <- x$res$p_val
  # res[is.na(res)] <- 1
  res
}

seurat_apv <- function(x) {
  res <- x$res$p_val_adj
  # res[is.na(res)] <- 1
  res
}

seurat_lfc <- function(x) {
  x$res$avg_log2FC
}

seurat_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$p_val)) }


### MAST


# MAST (on CPM counts)
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpm.R
MASTcpm_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  grp <- group
  dge <- edgeR::DGEList(counts = countData, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  lgcpms <- edgeR::cpm(dge, log = T, prior.count = 1)
  sca <- MAST::FromMatrix(exprsArray = lgcpms,
                          cData = data.frame(wellKey = colnames(lgcpms),
                                             grp = grp))

  zlmdata <- MAST::zlm(formula = ~ grp, sca = sca)
  mast <- MAST::lrTest(zlmdata, "grp")

  lfc <- MAST::getLogFC(zlmdata)

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = list(mast = mast, lfc = lfc),
       runtime = runtime$toc - runtime$tic)
}

MASTcpm_pv <- function(x) {
  res <- x$res$mast[, "hurdle", "Pr(>Chisq)"]
  # res[is.na(res)] <- 1
  res
}

MASTcpm_apv <- function(x) {
  res <- p.adjust(x$res$mast[, "hurdle", "Pr(>Chisq)"], method = "BH")
  # res[is.na(res)] <- 1
  res
}

MASTcpm_lfc <- function(x) {
  x$res$lfc$logFC
}

MASTcpm_time <- function(x) { rep(as.numeric(x$runtime), dim(x$res$mast)[1]) }


### Monocle


# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_monoclecount.R
monocle_run <- function(countData, group) {
  set.seed(12345)
  ## this is just to compute the LFCs. Do not count this computation time
  dge <- edgeR::DGEList(counts = countData, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  lgcpms <- edgeR::cpm(dge, log = T, prior.count = 1)

  tictoc::tic()

  mon <- monocle::newCellDataSet(as.matrix(countData),
                                 phenoData = new("AnnotatedDataFrame",
                                                 data = data.frame(condition = group,
                                                                   row.names = colnames(countData))),
                                 lowerDetectionLimit = 0.1, #0.5,
                                 expressionFamily = VGAM::negbinomial())
  mon <- BiocGenerics::estimateSizeFactors(mon)
  mon <- BiocGenerics::estimateDispersions(mon)
  monres <- monocle::differentialGeneTest(cds = mon,
                                          fullModelFormulaStr = " ~ condition",
                                          reducedModelFormulaStr = "~1")

  lfc <- rowMeans(lgcpms[, group[group == levels(group)[1]]]) -
    rowMeans(lgcpms[, group[group == levels(group)[2]]])

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = monres,
       lfc = lfc,
       runtime = runtime$toc - runtime$tic)
}

# monocle3_run <- function(countData, group) {
#
#     ## this is just to compute the LFCs. Do not count this computation time
#     dge <- edgeR::DGEList(counts = countData, group = group)
#     dge <- edgeR::calcNormFactors(dge, method = "TMM")
#     lgcpms <- edgeR::cpm(dge, log = T, prior.count = 1)
#
#     tictoc::tic()
#
#     cds <- monocle3::new_cell_data_set(expression_data = countData,
#                                        cell_metadata = data.frame(condition = group,
#                                                                   row.names = colnames(countData)))
#
#     # Perform differential expression analysis with regression:
#     gene_fits <- monocle3::fit_models(cds = cds,
#                                       # expression_family = "quasipoisson", #the default
#                                       expression_family = "negbinomial",
#                                       # expression_family = "poisson",
#                                       # expression_family = "binomial",
#                                       # expression_family = "gaussian",
#                                       # expression_family = "zipoisson",
#                                       # expression_family = "zinegbinomial",
#                                       model_formula_str = "~ condition")
#
#     fit_coefs <- monocle3::coefficient_table(gene_fits)
#     cond_terms <- fit_coefs %>% dplyr::filter(term == "condition2")
#     cond_terms <- cond_terms %>% dplyr::mutate(q_value = p.adjust(p_value))
#
#     # hist(cond_terms$p_value, 50)
#
#     res <- data.frame(pval = cond_terms$p_value,
#                       qval = cond_terms$q_value)
#     rownames(res) <- cond_terms$gene_id
#
#     lfc <- rowMeans(lgcpms[, group[group == levels(group)[1]]]) -
#         rowMeans(lgcpms[, group[group == levels(group)[2]]])
#
#     runtime <- tictoc::toc(log = F, quiet = T)
#
#     list(res = res,
#          lfc = lfc,
#          runtime = runtime$toc - runtime$tic)
# }

monocle_pv <- function(x) {
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

monocle_apv <- function(x) {
  res <- x$res$qval
  # res[is.na(res)] <- 1
  res
}

monocle_lfc <- function(x) {
  x$lfc
}

monocle_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$pval)) }


### Wilcoxon
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_Wilcoxon.R

wilcoxon_run <- function(countData, group) {
  set.seed(12345)
  tictoc::tic()

  dge <- edgeR::DGEList(counts = countData, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  cpms <- edgeR::cpm(dge)
  idx <- 1:nrow(cpms)
  names(idx) <- rownames(cpms)
  wilcox_p <- sapply(idx, function(i) {
    wilcox.test(countData[i, ] ~ group, exact = F)$p.value
  })

  # hist(wilcox_p, 50)

  res <- data.frame(pval = wilcox_p,
                    row.names = names(wilcox_p))

  lfc <- log2(rowMeans(cpms[, which(group == levels(group)[1])] + 1) /
                rowMeans(cpms[, which(group == levels(group)[2])] + 1))

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       lfc = lfc,
       runtime = runtime$toc - runtime$tic)
}

# wilcoxonLog_run <- function(countData, group) {
#
#     tictoc::tic()
#
#     dge <- edgeR::DGEList(counts = countData, group = group)
#     dge <- edgeR::calcNormFactors(dge, method = "TMM")
#     cpms <- edgeR::cpm(dge, log = T, prior.count = 1) ## log scale
#     idx <- 1:nrow(cpms)
#     names(idx) <- rownames(cpms)
#     wilcox_p <- sapply(idx, function(i) {
#         wilcox.test(countData[i, ] ~ group, exact = F)$p.value
#     })
#
#     # hist(wilcox_p, 50)
#
#     res <- data.frame(pval = wilcox_p,
#                       row.names = names(wilcox_p))
#
#     lfc <- rowMeans(cpms[, which(group == levels(group)[1])] + 1) -
#         rowMeans(cpms[, which(group == levels(group)[2])] + 1)
#
#     runtime <- tictoc::toc(log = F, quiet = T)
#
#     list(res = res,
#          lfc = lfc,
#          runtime = runtime$toc - runtime$tic)
# }

wilcoxon_pv <- function(x) {
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

wilcoxon_apv <- function(x) {
  res <- p.adjust(x$res$pval, method = "BH")
  # res[is.na(res)] <- 1
  res
}

wilcoxon_lfc <- function(x) {
  x$lfc
}

wilcoxon_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$pval)) }


scCODE_run <- function(countData, group) {
  ## scCODE installation
  # necessary1 <- c('doParallel', 'samr','doSNOW','pls')
  # installed <- necessary1 %in% installed.packages()[, 'Package']
  # if (length(necessary1[!installed]) >=1){
  #   install.packages(necessary1[!installed])
  # }
  #
  # necessary2<-c('DESeq2', 'DEsingle',
  #               'edgeR', 'limma', 'MAST',
  #               'S4Vectors', 'scDD', 'scmap',
  #               'SingleCellExperiment', 'SummarizedExperiment')
  # # installed <- necessary2 %in% installed.packages()[, 'Package']
  # #
  # # if (length(necessary2[!installed]) >=1){
  # #   if (!requireNamespace("BiocManager", quietly = TRUE))
  # #     install.packages("BiocManager")
  # #   library(BiocManager)
  # #   BiocManager::install(necessary2[!installed])
  # # }
  # ## exploit renv
  # install.packages(paste0("bioc::", necessary2))
  # install.packages(c("nghiavtr/BPSC"))
  # OGFSC_file <- "OGFSC_0.2.3.tar.gz"
  # download.file(url = "https://github.com/XZouProjects/OGFSC-R/raw/master/OGFSC_0.2.3.tar.gz",
  #               destfile = OGFSC_file)
  # install.packages(OGFSC_file, repos = NULL, type="source")
  # scCODE_file <- "scCODE_1.2.0.0.tar.gz"
  # download.file(url = "https://github.com/XZouProjects/scCODE/releases/download/1.2.0.0/scCODE_1.2.0.0.tar.gz",
  #               destfile = scCODE_file)
  # install.packages(scCODE_file, repos = NULL, type="source")
  library(scCODE)

  set.seed(12345)
  tictoc::tic()

  pvals <- setNames(rep_len(x = NA, length.out = nrow(countData)), rownames(countData))

  data1 <- countData[, which(group == levels(group)[1])]
  data2 <- countData[, which(group == levels(group)[2])]
  # There are two parameters selectable. "light", True or False,
  # default as True, run scCODE in a light version which saves time.
  # "top_ranked", the number of top-ranked strategies selected (5-10), default as 5.
  results <- scCODE(data1, data2, light = TRUE, top_ranked = 5)

  ## scCODE gives only the significant genes
  pvals[as.character(results$DE_results$Gene_name)] <- as.numeric(results$DE_results$P_adjust)

  res <- data.frame(pval = pvals[rownames(countData)],
                    row.names = rownames(countData))

  lfc <- setNames(rep_len(x = 0, length.out = nrow(countData)), rownames(countData))
  lfc[as.character(results$DE_results$Gene_name)] <- results$DE_results$logFC

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       lfc = lfc,
       runtime = runtime$toc - runtime$tic)
}

scCODE_pv <- function(x) {
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

scCODE_apv <- function(x) {
  ## scCODE only gives the P-adjusted
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

scCODE_lfc <- function(x) {
  x$lfc
}

scCODE_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$pval)) }

scDEA_run <- function(countData, group) {
  ## install
  # deps <- c("aggregation",
  #           "nghiavtr/BPSC",
  #           "statOmics/zingeR",
  #           "Zhangxf-ccnu/scDEA",
  #           paste0("bioc::",
  #                  c("sctransform", "DEsingle", "DESeq2", "edgeR", "MAST", "monocle", "scDD",
  #                    "limma", "Seurat", "SingleCellExperiment", "scater")))
  # install.packages(deps)

  library(scDEA)

  set.seed(12345)
  tictoc::tic()

  Pvals <- scDEA_individual_methods(raw.count = countData,
                                    cell.label = group,
                                    BPSC.parallel = T,
                                    DEsingle.parallel = T,
                                    DESeq2.parallel = T,
                                    MAST.parallel = T)
  combination.Pvals <- lancaster.combination(Pvals = Pvals, weight = TRUE, trimmed = 0.2)
  adjusted.Pvals <- scDEA.p.adjust(combination.Pvals, adjusted.method = "BH")

  res <- data.frame(pval = combination.Pvals[rownames(countData)],
                    padj = adjusted.Pvals[rownames(countData)],
                    row.names = rownames(countData))

  lfc <- setNames(rep_len(x = 0, length.out = nrow(countData)), rownames(countData))
  # lfc[as.character(results$DE_results$Gene_name)] <- results$DE_results$logFC

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       lfc = lfc,
       runtime = runtime$toc - runtime$tic)

}

scDEA_pv <- function(x) {
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

scDEA_apv <- function(x) {
  ## scCODE only gives the P-adjusted
  res <- x$res$padj
  # res[is.na(res)] <- 1
  res
}

scDEA_lfc <- function(x) {
  x$lfc
}

scDEA_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$pval)) }

hytest_run <- function(countData, group) {
  # download.file(url = "https://github.com/gianluca-sottile/A-Novel-Statistical-Test-For-Differential-Expression-Analysis/raw/main/0.hy.test.R",
  #               destfile = "hy_test.R")
  #https://github.com/gianluca-sottile/A-Novel-Statistical-Test-For-Differential-Expression-Analysis

  ##### hy-test #####
  source("../utils/hy_test.R")
  # n <- ncol(dati.geni.log) / 2
  # tissue <- rep(1:2, each = n)

  set.seed(12345)
  tictoc::tic()

  dge <- edgeR::DGEList(counts = countData, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  lcpms <- edgeR::cpm(dge, log = T, prior.count = 1)

  dati.geni.log <- lcpms
  tissue <- group

  dati.geni.log.tolist <- as.list(data.frame(t(dati.geni.log)))

  # out <- pbmclapply(dati.geni.log.tolist,
  #                   function(gene, tissue) hy.test(gene ~ tissue, data = data.frame(gene, tissue)),
  #                   tissue = tissue, mc.cores = 4L)
  # out <- BiocParallel::bplapply(dati.geni.log.tolist,
  #               function(gene, tissue) hy.test(gene ~ tissue,
  #                                              data = data.frame(gene, tissue)),
  #               tissue = tissue,
  #               BPPARAM = BiocParallel::MulticoreParam(BiocParallel::multicoreWorkers()))
  out <- lapply(dati.geni.log.tolist,
                function(gene, tissue) hy.test(gene ~ tissue,
                                               data = data.frame(gene, tissue)),
                tissue = tissue)

  pvals <- sapply(out, function(x) x$p.value)

  res <- data.frame(pval = pvals,
                    padj = p.adjust(pvals, method = "BH"),
                    row.names = rownames(countData))

  ## mock LFC
  lfc <- setNames(rep_len(x = 0,
                          length.out = nrow(countData)),
                  rownames(countData))

  runtime <- tictoc::toc(log = F, quiet = T)

  list(res = res,
       lfc = lfc,
       runtime = runtime$toc - runtime$tic)
}

hytest_pv <- function(x) {
  res <- x$res$pval
  # res[is.na(res)] <- 1
  res
}

hytest_apv <- function(x) {
  res <- x$res$padj
  # res[is.na(res)] <- 1
  res
}

hytest_lfc <- function(x) {
  x$lfc
}

hytest_time <- function(x) { rep(as.numeric(x$runtime), length(x$res$pval)) }
