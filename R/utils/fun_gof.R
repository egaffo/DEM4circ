suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pscl))
suppressPackageStartupMessages(library(MASS)) ## for the glm.nb
suppressPackageStartupMessages(library(ggpointdensity))

### Extraction of ZINB coefs
computeExp <- function(zinbModel){
  (1 - t(getPi(zinbModel))) * t(getMu(zinbModel))
}
computeVar <- function(zinbModel){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  (1 - pi) * mu * (1 + mu*(phi + pi))
}

computeP0 <- function(zinbModel){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  pi + (1 - pi) * (1 + phi * mu) ^ (-1/phi)
}

zinb.loglik.matrix <- function(model, x) {
  mu <- getMu(model)
  theta <- getTheta(model)
  theta_mat <- matrix(rep(theta, each = nrow(x)), ncol = ncol(x))
  pi <- getPi(model)
  lik <- pi * (x == 0) + (1 - pi) * dnbinom(x, size = theta_mat, mu = mu)
  lik[lik == 0] <- min(lik[lik != 0]) #to avoid log lik to be infinite
  log(lik)
}

# Negative Binomial fitting with TMM norm factors
fitNB_TMM <- function(counts, design){
  normFacts <- edgeR::calcNormFactors(counts, method = "TMM")
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- normFacts
  disp <- estimateDisp(y = dge,design = design,tagwise = TRUE)
  fit <- glmFit(dge$counts,dispersion = disp$tagwise.dispersion, design = design)
  list(fitted = fit$fitted.values, disp = disp$tagwise.dispersion)
}

#' Scale counts and compute mean number of zero counts per row
#'
#' @param counts
#' @param scale _NULL_ (default) no scaling; _median_ median scaling;
#' _default_ is CPM
#'
#' @return
#' @export
#'
#' @examples
prepareObserved <- function(counts,
                            scale = NULL){
  if(!is.null(scale)){

    if(scale == "median"){

      counts <- counts * stats::median(colSums(counts)) / colSums(counts)

    } else if(scale == "default"){

      counts <- counts * 1e6 / colSums(counts)

    } else stop("When specified, 'scale' must be 'median' or 'default'")
  }

  Y <- log1p(rowMeans(counts))
  Y0 <- rowMeans(counts == 0)

  return(data.frame("Y" = Y, "Y0" = Y0))
}

meanDifferences <- function(estimated,
                            observed){
  if(nrow(estimated) != nrow(observed))
    stop("Estimated and Observed data.frames have different number of rows")
  MD <- estimated$Y - observed$Y
  ZPD <- estimated$Y0 - observed$Y0
  return(data.frame("MD" = MD, "ZPD" = ZPD))
}

RMSE <- function(differences){
  sqrt(mean(differences^2,na.rm = TRUE))
}

# Negative Binomial fitting
fitNB <- function(dge){

  message("Fitting NB ...")

  fit <- edgeR::glmFit(dge$counts,
                       dispersion = dge$tagwise.dispersion)
  nbloglik <- rowSums(dnbinom(x = dge$counts,
                              size = 1 / dge$tagwise.dispersion,
                              mu = rowMeans(fit$fitted.values),
                              log = TRUE))
  # Fitted values extraction
  # Return the log(average fitted values + 1) to

  Y <- log1p(rowMeans(fit$fitted.values))

  Y0 <- rowMeans((1 + fit$fitted.values * dge$tagwise.dispersion)^(-1 / dge$tagwise.dispersion))

  return(data.frame("Y" = Y,
                    "Y0" = Y0,
                    "nbloglik" = nbloglik,
                    row.names = names(Y)))
}

fit_zeroinfl <-
  function(ithrow, coefs, disps, sf, logmean){

    zinb.mean <- exp(logmean)
    zinb.prop <- -Inf
    zinb.disp <- disps
    zinb.loglik <- NA

    if(any(ithrow == 0)){
      ## if at least a count is zero, fit a ZINB
      try({
        zfit <- zeroinfl(formula = ithrow ~ 1 | 1,
                         dist = "negbin",
                         offset = log(sf),
                         EM = F,
                         start = list(count = coefs,
                                      zero = -3,
                                      theta = 1 / disps))

        zinb.mean <- mean(exp(zfit$coefficients$count))
        zinb.prop <- zfit$coefficients$zero
        zinb.disp <- 1 / zfit$theta
        zinb.loglik <- zfit$loglik
      })
    } else {
      ## fit a NB if no zero counts
      zfit <- MASS::glm.nb(formula = c(ithrow) ~ 1)

      zinb.mean <- mean(exp(zfit$coefficients))
      zinb.prop <- -Inf
      zinb.disp <- 1 / zfit$theta
      zinb.loglik <- zfit$twologlik / 2
    }

    list("zinb.mean" = zinb.mean,
         "zinb.prop" = zinb.prop,
         "zinb.disp" = zinb.disp,
         "zinb.loglik" = zinb.loglik)
  }

fitZINB <- function(dge) {

  message("Fitting ZINB ...")

  sf <- dge$samples[, "norm.factors"]
  coefs <- coef(glmFit(dge$counts,
                       dispersion = dge$tagwise.dispersion))

  centered.off <- getOffset(dge)
  centered.off <- centered.off - mean(centered.off)
  logmeans <- mglmOneGroup(dge$counts,
                           offset = centered.off,
                           dispersion = dge$tagwise.dispersion)

  library(BiocParallel)
  fitted <-
    rbindlist(bpmapply(FUN = fit_zeroinfl,
                       ithrow = asplit(dge$counts, 1),
                       coefs = coefs,
                       disps = dge$tagwise.dispersion,
                       MoreArgs = list(sf = sf),
                       logmean = logmeans,
                       SIMPLIFY = F),
              idcol = "circ_id")

  zinb.prop <- exp(fitted$zinb.prop) / (1 + exp(fitted$zinb.prop))

  EY <- log1p((1 - zinb.prop) * fitted$zinb.mean)

  EY0 = zinb.prop + (1 - zinb.prop) * (1 + fitted$zinb.disp * fitted$zinb.mean) ^ (-1 / fitted$zinb.disp)

  return(data.frame("Y" = EY,
                    "Y0" = EY0,
                    "zinb.loglik" = fitted$zinb.loglik,
                    row.names = fitted$circ_id))
}

#' Fit NB and/or ZINB models
#'
#' @param counts The count matrix: genes in rows, samples as columns
#' @param models A character specifying which models to fit ("NB" or "ZINB") or
#' the vector c("NB", "ZINB") if both models has to be fitted
#' @param coldata As default, a design matrix with only an intercept is used to
#' compute the dispersion. Set this parameter with a data.frame reporting one
#' _condition_ column to consider experimental design.
#' @param ... additional parameters for internal functions.
#' F.i. set scale = "median"
#'
#' @return
#' @export
#'
#' @examples
fitModels <- function(counts,
                      models = c("NB", "ZINB"),
                      coldata = NULL,
                      ...) {

  nCovariates <- 0

  if(is.null(coldata)){

    message("NULL coldata. A model with only the intercept will be used")
    ## generate a single-condition model
    coldata <-
      data.frame(condition = rep("A", dim(counts)[2]))
    rownames(coldata) <- colnames(counts)

    design <- model.matrix( ~ 1, data = coldata)

  }else{

    design <- model.matrix( ~ condition, data = coldata)

  }

  if (!is.null(design)) {

    covariates <- as.matrix(design)

    p <- dim(design)[2]

    if(p == 1 && length(unique(design)) == 1){

      # only intercept
      message("Only the intercept in design")
      nCovariates <- 1

    } else {

      message(paste(p, "covariates in design"))
      nCovariates <- p + 1
    }

    message(paste("Number of covariates:", nCovariates))
  }

  observed <- prepareObserved(counts, ...)

  dge <- DGEList(counts = counts,
                 group = coldata[colnames(counts), "condition"])
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, tagwise = T, design = design)

  nObs <- ncol(dge$counts)

  fittedModels <- list()

  if("NB" %in% models){

    fitted <- fitNB(dge)

    MD <- meanDifferences(estimated = fitted,
                          observed = observed)

    aicNB <- -2 * fitted$nbloglik + 2 * (nCovariates + 1)

    bicNB <- -2 * fitted$nbloglik + (log(nObs) * (nCovariates + 1))

    fittedModels$NB <- data.frame(observed, ## 'observed' is a named list
                                  MD, ## 'MD' is a named list
                                  EY = fitted$Y,
                                  EY0 = fitted$Y0,
                                  loglik = fitted$nbloglik,
                                  aic = aicNB,
                                  bic = bicNB,
                                  circ_id = rownames(fitted))
  }

  if("ZINB" %in% models){

    fitted <- fitZINB(dge)

    MD <- meanDifferences(estimated = fitted,
                          observed = observed)

    aicZINB <- -2 * fitted$zinb.loglik + 2 * (nCovariates + 1)

    bicNB <- -2 * fitted$zinb.loglik + (log(nObs) * (nCovariates + 1))

    fittedModels$ZINB <- data.frame(observed, ## 'observed' is a named list
                                    MD, ## 'MD' is a named list
                                    EY = fitted$Y,
                                    EY0 = fitted$Y0,#fitted$zinb.prop,
                                    loglik = fitted$zinb.loglik,
                                    aic = aicZINB,
                                    bic = bicNB,
                                    circ_id = rownames(fitted))
  }

  return(fittedModels)
}

extract_values <- function(site,
                           distributions = c("ZINB", "NB"),
                           varnames = c("EY", "EY0")){

  cbind(ldply(lapply(site[distributions],
                     function(dist){
                       as.data.frame(dist[varnames])
                     }),
              .id = "Models"),
        ldply(lapply(site["OBSERVED"],
                     function(dist){
                       observed_df <- as.data.frame(dist)
                       return(observed_df)
                     })))
}

# From model list to data.frame
model_to_data.frame <- function(model_list){
  # df_list <- lapply(model_list,extract_values)
  # df <- ldply(df_list)
  df_list = lapply(model_list, function(x) {
    temp = rbindlist(x, idcol = "Models")
    temp$circ_id = c(rownames(x$NB), rownames(x$ZINB))
    temp
  })
  df = ldply(df_list)
  colnames(df) <- c("Method", "Model",
                    "Y", "Y0", "MD", "ZPD", "EY", "EY0",
                    "loglik", "AIC", "circ_id")
  return(df)
}


# GOF indexes
# Root mean square errors
compute_RMSE <- function(df){
  ddply(df,.variables = ~ Models + Methods, function(x){
    MD <- sqrt(mean(x$mean_diff^2,na.rm = TRUE))
    # MD <- sqrt(x$mean_diff^2)
    MD_sd <- sd(MD)
    ZPD <- sqrt(mean(x$zero_prob^2,na.rm = TRUE))
    # ZPD <- sqrt(abs(x$zero_prob)^2)
    ZPD_sd <- sd(ZPD)
    return(data.frame("MD" = MD, "ZPD" = ZPD))
  })
}

compute_MAE <- function(df){
  ddply(df,.variables = ~ Models + Methods, function(x){
    MD <- mean(x$mean_diff,na.rm = TRUE)
    # MD <- sqrt(x$mean_diff^2)
    MD_sd <- sd(MD)
    ZPD <- mean(x$zero_prob,na.rm = TRUE)
    # ZPD <- sqrt(abs(x$zero_prob)^2)
    ZPD_sd <- sd(ZPD)
    return(data.frame("MD" = MD, "ZPD" = ZPD))
  })
}

RMSE.func <- function(differences){
  sqrt(mean(differences^2, na.rm = TRUE))
}


MDPlot <- function(data, difference = NULL, split = TRUE){

  if(difference == "MD"){
    if("Model" %in% colnames(data) &
       "Method" %in% colnames(data) &
       "Y" %in% colnames(data) &
       "MD" %in% colnames(data)){

      RMSE_MD <- plyr::ddply(.data = data,
                             .variables = ~ Method + Model,
                             .fun = function(m) cbind("RMSE" = RMSE.func(m$MD)))

      gobj <- ggplot(data = data, aes(x = Y, y = MD)) + #, color = Model
        # scale_color_manual(values = method_cols) +
        ggtitle(label = "Mean Differences plot",
                subtitle = paste0("Observed = log(mean(counts*)+1)",
                                  "\n",
                                  "Estimated = log(mean(fitted*)+1)"))

      if(split){
        gobj <- gobj +
          geom_text(data = RMSE_MD, color = "black",
                    aes(x = mean(data$Y) + 2,
                        y = max(data$MD,na.rm = TRUE),
                        label = paste0("RMSE:", round(RMSE, 3))), size = 3)
      }

    } else {
      stop(paste0("'data' should contains 'Method', Model', 'Y', and 'MD' ",
                  "columns for the method that generated the counts, ",
                  "the model name, observed values, and mean difference ",
                  "values, respectively."))}

  } else if(difference == "ZPD"){

    if("Model" %in% colnames(data) &
       "Method" %in% colnames(data) &
       "Y0" %in% colnames(data) &
       "ZPD" %in% colnames(data)){

      RMSE_ZPD <- ddply(.data = data,
                        .variables = ~ Method + Model,
                        .fun = function(m) cbind("RMSE" = RMSE.func(m$ZPD)))

      gobj <- ggplot(data = data, aes(x = Y0,  y = ZPD)) + #,  color = Model
        # scale_color_manual(values = method_cols) +
        ggtitle(label = "Zero Probability Differences plot", subtitle =
                  "Observed = mean(counts=0)\nEstimated = mean(P(Y=0))")

      if(split){
        gobj <- gobj +
          geom_text(data = RMSE_ZPD, color = "black",
                    aes(x = 0.25,
                        y = max(data$ZPD,na.rm = TRUE),
                        label = paste0("RMSE:", round(RMSE, 2))), size = 3)
      }

    } else {
      stop(paste0("'data' should contains 'Method', Model', 'Y0', and 'ZPD' ",
                  "columns for the method that generated the counts, ",
                  "the model name, zero rate observed values, and zero probability ",
                  "values, respectively."))
    }
  } else stop("Difference must be 'MD' or 'ZPD'")

  gobj <- gobj +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    theme(legend.position = "bottom") +
    xlab("Observed") +
    ylab("Estimated-Observed")

  if(length(unique(data$Model)) > 1){
    if(split){

      gobj <- gobj +
        facet_grid(Model ~ Method, labeller = label_wrap_gen(width = 10)) +
        # geom_point(alpha = .3) +
        geom_pointdensity(size = .5) +
        geom_smooth(color = "gray", method = "lm") +
        theme_classic() +
        theme(legend.position = "bottom",
              plot.margin = unit(c(0, 0, 0, 0), "cm"))
    } else {
      gobj <- gobj +
        geom_smooth()
    }
  }

  return(gobj)
}


RMSEPlot <- function(data, difference = NULL, method_cols){
  if(difference == "MD"){
    if("Model" %in% colnames(data) & "Y" %in% colnames(data) &
       "MD" %in% colnames(data)){

      RMSE <- plyr::ddply(.data = data,
                          .variables = ~ Model + Method,
                          .fun = function(m) cbind("RMSE" = RMSE.func(m$MD)))

      gobj <- ggplot(data = RMSE, aes(x = Model, y = RMSE, fill = Model)) +
        scale_fill_manual(values=method_cols) +
        geom_col() +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE,2)),
                   fill = "white") +
        ggtitle(label = "RMSE",
                subtitle = "Mean differences")

    } else {stop("data should contains 'Model', 'Y', and 'MD' columns for model
              name, observed values and mean difference values respectively.")}
  } else if(difference == "ZPD"){
    if("Model" %in% colnames(data) &
       "Y0" %in% colnames(data) &
       "ZPD" %in% colnames(data)){

      RMSE <- plyr::ddply(.data = data,
                          .variables = ~ Model + Method,
                          .fun = function(m) cbind("RMSE" = RMSE.func(m$ZPD)))

      gobj <- ggplot(data = RMSE, aes(x = Model, y = RMSE, fill = Model)) +
        geom_col() +
        scale_fill_manual(values=method_cols) +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE,4)),
                   fill = "white") +
        ggtitle(label = "RMSE",
                subtitle = "Zero probability difference")

    } else {stop("df should contains 'Model', 'Y0', and 'ZPD' columns for model
              name, zero rate observed values and zero probability difference
              values respectively.")}

  } else stop("Difference must be 'MD' or 'ZPD'")

  gobj <- gobj +
    facet_grid( ~ Method, labeller = labeller(.cols = label_both)) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete() +
    theme_classic() +
    theme(legend.position = "bottom",
          strip.text.x = element_text(face = "bold.italic", size = 5),
          strip.text.y = element_text(face = "bold.italic", size = 8), #angle = 75),
          strip.background = element_rect(#fill = "lightblue",
            colour = "grey", size = 1),
          axis.text.x = element_text(hjust = 1, angle = 45),
          # axis.text.x = element_blank(),
          # axis.text.y = element_blank(),
          # panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0,0,0,0), "cm"))

  return(gobj)
}
