#' Calculate statistics for pairwise comparison of data sets
#'
#' Calculate a range of statistics and p-values for comparison of two data sets.
#'
#' @param df The input data frame. Must contain at least a column named
#'   'dataset' and an additional column with values to use as the basis for the
#'   comparison.
#' @param ds1,ds2 The names of the two data sets to be compared.
#' @param column The name of the column(s) of \code{df} to be used as the basis
#'   for the comparison.
#' @param subsampleSize The number of observations for which certain
#'   time-consuming statistics will be calculated. The observations will be
#'   selected randomly among the rows of \code{df}.
#' @param permute Whether to permute the dataset column of \code{df} before
#'   calculating the statistics.
#' @param kmin,kfrac For statistics that require the extraction of k nearest
#'   neighbors of a given point, the number of neighbors will be max(kmin, kfrac
#'   * nrow(df)).
#' @param xmin,xmax Smallest and largest value of \code{column}, used to
#'   normalize the x-axis when calculating the area between the eCDFs.
#'
#' @return A vector with statistics and p-values
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
#' @importFrom stats ks.test ecdf chisq.test
#'
calculateStats <- function(df, ds1, ds2, column, subsampleSize,
                           permute = FALSE, kmin = 5, kfrac = 0.01, xmin, xmax) {

    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(stats))

    ## Remove rows with NA values in column(s) of interest
    df <- df[rowSums(is.na(df[, column, drop = FALSE])) == 0, ]

    ## Permute dataset column if requested
    if (permute) {
        df$dataset <- sample(df$dataset, nrow(df))
    }

    if (length(column) == 1) {
        ## Kolmogorov-Smirnov test
        kst <- stats::ks.test(df[, column][df$dataset == ds1],
                              df[, column][df$dataset == ds2])

        ## Area between ecdfs
        e1 <- stats::ecdf(df[, column][df$dataset == ds1])
        e2 <- stats::ecdf(df[, column][df$dataset == ds2])
        xv <- sort(unique(df[, column]))
        ediff <- abs(e1(xv) - e2(xv))
        earea <- caTools::trapz((xv - xmin)/(xmax - xmin), ediff)
    }

    ## Silhouette width and nearest neighbour distribution for subsample of
    ## observations
    ## Calculate overall proportion of values from the two data sets
    overallprop <- table(factor(df$dataset, levels = c(ds1, ds2)))
    overallprop <- overallprop/sum(overallprop)
    nrep <- 1
    chisilh <- t(vapply(seq_len(nrep), function(m) {
        ## Select subset of observations for calculation of silhouette widths and
        ## nearest neighbor distributions
        idx <- sample(seq_len(nrow(df)), min(nrow(df), subsampleSize))

        ## For each selected observation, calculate statistics
        tmp <- t(vapply(idx, function(j) {

            ## Distances to all observations
            dists <- sqrt(rowSums((sweep(df[, column, drop = FALSE],
                                         2, unlist(df[j, column,
                                                      drop = FALSE]), "-")) ^ 2))
            dists_this <- dists[df$dataset == df$dataset[j]]
            dists_other <- dists[df$dataset != df$dataset[j]]

            ## Silhouette width
            silh <- (mean(dists_other, na.rm = TRUE) -
                         mean(dists_this, na.rm = TRUE)) /
                max(mean(dists_other, na.rm = TRUE), mean(dists_this, na.rm = TRUE))

            ## Local silhouette width
            dists_this_local <-
                sort(dists_this)[seq_len(max(kmin, kfrac * nrow(df)) *
                                             overallprop[df$dataset[j]])]
            dists_other_local <-
                sort(dists_other)[seq_len(max(kmin, kfrac * nrow(df)) *
                                              (1 - overallprop[df$dataset[j]]))]
            silh_local <-
                (mean(dists_other_local, na.rm = TRUE) -
                     mean(dists_this_local, na.rm = TRUE)) /
                max(mean(dists_other_local, na.rm = TRUE),
                    mean(dists_this_local, na.rm = TRUE))
            if (all(c(dists_this_local, dists_other_local) == 0)) {
                silh_local <- 0
            }

            ## Chi-square test comparing distribution of data set labels among
            ## nearest neighbors to the overall distribution
            x <- df$dataset[order(dists)][seq_len(max(kmin, kfrac * nrow(df)))]
            chisqp <- stats::chisq.test(x = table(factor(x, levels = c(ds1, ds2))),
                                        p = overallprop)$p.value
            ## Return values
            c(silh = silh, chisqp = chisqp, silh_local = silh_local)
        }, c(silh = NA_real_, chisqp = NA_real_, silh_local = NA_real_)))
        c(silh = mean(tmp[, "silh"], na.rm = TRUE),
          silhlocal = mean(tmp[, "silh_local"], na.rm = TRUE),
          chisqp = mean(tmp[, "chisqp"] <= 0.05, na.rm = TRUE))
    }, c(silh = NA_real_, silhlocal = NA_real_, chisqp = NA_real_)))

    ## Runs test
    if (length(column) == 1) {
        # df <- df %>% arrange_(column)
        df <- df %>% arrange(column)
        runs_res <- randtests::runs.test(as.numeric(df$dataset == ds1),
                                         threshold = 0.5,
                                         alternative = "left.sided", plot = FALSE)
    }

    if (length(column) == 1)
        c(ksstatistic = signif(kst$statistic, 3),
          kspvalue = signif(kst$p.value, 3),
          diffarea = signif(earea, 3),
          runsstatistic = signif(runs_res$statistic, 3),
          runspvalue = signif(runs_res$p.value, 3),
          NNmismatch = signif(mean(chisilh[, "chisqp"]), 3),
          avesilh = signif(mean(chisilh[, "silh"]), 3),
          avesilhlocal = signif(mean(chisilh[, "silhlocal"]), 3))
    else
        c(NNmismatch = signif(mean(chisilh[, "chisqp"]), 3),
          avesilh = signif(mean(chisilh[, "silh"]), 3),
          avesilhlocal = signif(mean(chisilh[, "silhlocal"]), 3))
}

#' Construct data frame with pairwise statistics
#'
#' Construct a data frame containing statistics and p-values for pairwise
#' comparison of data sets.
#'
#' @param df The input data frame. Must contain at least a column named
#'   'dataset' and an additional column with values
#' @param column The name of the column(s) of \code{df} to be used as the basis
#'   for the comparison
#' @param permutationPvalues Whether or not to calculate p-values of statistics
#'   via permutation
#' @param nPermutations The number of permutations (only used if
#'   permutationPvalues = TRUE)
#' @param subsampleSize The number of observations for which certain
#'   (time-consuming) statistics will be calculated. The observations will be
#'   selected randomly among the rows of \code{df}
#' @param kmin,kfrac For statistics that require the extraction of k nearest
#'   neighbors of a given point, the number of neighbors will be max(kmin, kfrac
#'   * nrow(df))
#'
#' @return A data table with statistics and p-values for pairwise comparisons of
#'   data sets, based on the provided \code{column}
#'
#' @author Charlotte Soneson
#' @keywords internal
#'
#' @import dplyr
#' @importFrom utils combn
#'
makeDF <- function(df, column, permutationPvalues, nPermutations,
                   subsampleSize, kmin, kfrac,
                   bpparam = NULL) {

    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(utils))
    suppressPackageStartupMessages(library(BiocParallel))

    if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()

    ## Generate data frame with comparisons of all pairs of data sets
    if (length(unique(df$dataset)) > 1) {
        ## Initialize data frame by populating it with all data set pairs
        ds_df <- data.frame(t(utils::combn(unique(df$dataset), 2)),
                            stringsAsFactors = FALSE) %>%
            dplyr::rename(dataset1 = X1, dataset2 = X2)

        ## Calculate statistics for each pair of data sets
        ds_res <- as.data.frame(t(vapply(seq_len(nrow(ds_df)), function(i) {
            dftmp <- df %>% dplyr::filter(dataset %in%
                                              c(ds_df$dataset1[i], ds_df$dataset2[i]))

            obs_stats <- calculateStats(df = dftmp, ds1 = ds_df$dataset1[i],
                                        ds2 = ds_df$dataset2[i], column = column,
                                        subsampleSize = subsampleSize,
                                        permute = FALSE, kmin = kmin, kfrac = kfrac,
                                        xmax = max(df[, column], na.rm = TRUE),
                                        xmin = min(df[, column], na.rm = TRUE))
            if (permutationPvalues) {
                # perm_stats <- t(vapply(seq_len(nPermutations),
                #                        function(s) {
                #                            calculateStats(df = dftmp,
                #                                           ds1 = ds_df$dataset1[i],
                #                                           ds2 = ds_df$dataset2[i],
                #                                           column = column,
                #                                           subsampleSize = subsampleSize,
                #                                           permute = TRUE,
                #                                           kmin = kmin,
                #                                           kfrac = kfrac,
                #                                           xmax = max(df[, column], na.rm = TRUE),
                #                                           xmin = min(df[, column], na.rm = TRUE))
                #                        },
                #                        defaultStats(length(column),
                #                                     withP = FALSE)))
                ## doing in parallel is very dangerous as it will fill up quickly all the RAM
                ## Using bpparam = BiocParallel::SerialParam() should be safe
                perm_stats_l <-
                    bplapply(seq_len(nPermutations),
                             function(s, calculateStats) {
                                 calculateStats(df = dftmp, ds1 = ds_df$dataset1[i],
                                                ds2 = ds_df$dataset2[i], column = column,
                                                subsampleSize = subsampleSize,
                                                permute = TRUE, kmin = kmin, kfrac = kfrac,
                                                xmax = max(df[, column], na.rm = TRUE),
                                                xmin = min(df[, column], na.rm = TRUE))
                             },
                             calculateStats = calculateStats,
                             BPPARAM = bpparam)
                perm_stats <- data.frame(rbindlist(lapply(perm_stats_l, function(x)data.table(t(x)))))
                ## Permutation p-values
                obs_stats["NNmismatchP"] <-
                    signif(mean(c(obs_stats["NNmismatch"], perm_stats[, "NNmismatch"]) >=
                                    obs_stats["NNmismatch"]), 3)
                obs_stats["avesilhP"] <-
                    signif(mean(c(obs_stats["avesilh"], perm_stats[, "avesilh"]) >=
                                    obs_stats["avesilh"]), 3)
                obs_stats["avesilhlocalP"] <-
                    signif(mean(c(obs_stats["avesilhlocal"],
                                  perm_stats[, "avesilhlocal"]) >=
                                    obs_stats["avesilhlocal"]), 3)
                if (length(column) == 1)
                    obs_stats["diffareaP"] <-
                    signif(mean(c(obs_stats["diffarea"], perm_stats[, "diffarea"]) >=
                                    obs_stats["diffarea"]), 3)
            }
            obs_stats
        }, defaultStats(length(column), withP = permutationPvalues))))

        # ## Populate the output table
        # if (length(column) == 1) {
        #     ds_df$`K-S statistic` <- ds_res$ksstatistic
        #     ds_df$`K-S p-value` <- ds_res$kspvalue
        #     ds_df$`Scaled area between eCDFs` <- ds_res$diffarea
        #     if (permutationPvalues) {
        #         ds_df$`Scaled area between eCDFs perm p-value` <- ds_res$diffareaP
        #     }
        #     ds_df$`Runs statistic` <- ds_res$runsstatistic
        #     ds_df$`Runs p-value` <- ds_res$runspvalue
        # }
        # ds_df$`NN rejection fraction` <- ds_res$NNmismatch
        # if (permutationPvalues) {
        #     ds_df$`NN rejection fraction perm p-value` <- ds_res$NNmismatchP
        # }
        # ds_df$`Average silhouette width` <- ds_res$avesilh
        # if (permutationPvalues) {
        #     ds_df$`Average silhouette width perm p-value` <- ds_res$avesilhP
        # }
        # ds_df$`Average local silhouette width` <- ds_res$avesilhlocal
        # if (permutationPvalues) {
        #     ds_df$`Average local silhouette width perm p-value` <-
        #         ds_res$avesilhlocalP
        # }

        # ## Return output table
        # DT::datatable(ds_df, options = list(
        #     scrollX = TRUE,
        #     columnDefs = list(list(className = "dt-center", targets = 3:ncol(ds_df)))
        # ))

        ds_res
    }
}

#' Return a vector of NA scores
#'
#' @param n Number of columns to use for the comparison
#' @param withP Whether or not to include p-value columns
#'
#' @return A vector with NA values for all applicable statistics
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
defaultStats <- function(n, withP = FALSE) {
    if (n == 1) {
        if (withP) {
            c(ksstatistic = NA_real_, kspvalue = NA_real_, diffarea = NA_real_,
              runsstatistic = NA_real_, runspvalue = NA_real_, NNmismatch = NA_real_,
              avesilh = NA_real_, avesilhlocal = NA_real_, NNmismatchP = NA_real_,
              avesilhP = NA_real_, avesilhlocalP = NA_real_, diffareaP = NA_real_)
        } else {
            c(ksstatistic = NA_real_, kspvalue = NA_real_, diffarea = NA_real_,
              runsstatistic = NA_real_, runspvalue = NA_real_, NNmismatch = NA_real_,
              avesilh = NA_real_, avesilhlocal = NA_real_)
        }
    } else {
        if (withP) {
            c(NNmismatch = NA_real_, avesilh = NA_real_, avesilhlocal = NA_real_,
              NNmismatchP = NA_real_, avesilhP = NA_real_, avesilhlocalP = NA_real_)
        } else {
            c(NNmismatch = NA_real_, avesilh = NA_real_, avesilhlocal = NA_real_)
        }
    }
}

##
parallel_calcStats <- function(df, column, permutationPvalues, nPermutations,
                               subsampleSize, kmin, kfrac,
                               ref_ds,
                               bpparam = NULL) {

    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(utils))
    suppressPackageStartupMessages(library(BiocParallel))

    bpparam_ds <- bpparam
    cpus4permut = 1

    if(is.null(bpparam)){
        bpparam <- BiocParallel::SerialParam()
    } else {

        cpus4permut <- tryCatch({
            bpparam$resources$ncpus},
            error = function(x){1}
        )
    }

    df_l <- split(df, df$dataset)

    ## Calculate statistics for each target:source pair of data sets
    ds_res_l <- bplapply(df_l, function(x, ref, calculateStats,
                                        column,
                                        subsampleSize,
                                        kmin,
                                        kfrac,
                                        df,
                                        permutationPvalues,
                                        nPermutations,
                                        cpus4permut) {

        dftmp <- data.frame(data.table::rbindlist(list(x, ref), use.names = T, fill = T))

        obs_stats <- calculateStats(df = dftmp,
                                    ds1 = unique(x$dataset),
                                    ds2 = unique(ref$dataset),
                                    column = column,
                                    subsampleSize = min(subsampleSize, nrow(dftmp)),
                                    permute = FALSE,
                                    kmin = kmin,
                                    kfrac = kfrac,
                                    xmax = max(df[, column], na.rm = TRUE),
                                    xmin = min(df[, column], na.rm = TRUE))

        if (permutationPvalues) {
            ## doing in parallel is very dangerous as it will fill up quickly all the RAM
            ## Using bpparam = BiocParallel::SerialParam() should be safe
            perm_stats_l <-
                BiocParallel::bplapply(seq_len(nPermutations),
                                       function(s, dftmp, ref, calculateStats,
                                                column,
                                                subsampleSize,
                                                kmin,
                                                kfrac, df, cpus4permut) {

                                           calculateStats(df = dftmp,
                                                          ds1 = unique(x$dataset),
                                                          ds2 = unique(ref$dataset),
                                                          column = column,
                                                          subsampleSize = subsampleSize,
                                                          permute = TRUE,
                                                          kmin = kmin,
                                                          kfrac = kfrac,
                                                          xmax = max(df[, column], na.rm = TRUE),
                                                          xmin = min(df[, column], na.rm = TRUE))
                                       },
                                       ref = ref,
                                       calculateStats = calculateStats,
                                       column = column,
                                       kmin = kmin,
                                       kfrac = kfrac,
                                       subsampleSize = subsampleSize, df = df, dftmp = dftmp,
                                       BPPARAM = BiocParallel::MulticoreParam(workers = cpus4permut))

            perm_stats <- data.frame(data.table::rbindlist(lapply(perm_stats_l, function(x)data.table::data.table(t(x)))))

            ## Permutation p-values
            obs_stats["NNmismatchP"] <-
                signif(mean(c(obs_stats["NNmismatch"], perm_stats[, "NNmismatch"]) >=
                                obs_stats["NNmismatch"]), 3)
            obs_stats["avesilhP"] <-
                signif(mean(c(obs_stats["avesilh"], perm_stats[, "avesilh"]) >=
                                obs_stats["avesilh"]), 3)
            obs_stats["avesilhlocalP"] <-
                signif(mean(c(obs_stats["avesilhlocal"],
                              perm_stats[, "avesilhlocal"]) >=
                                obs_stats["avesilhlocal"]), 3)
            if (length(column) == 1)
                obs_stats["diffareaP"] <-
                signif(mean(c(obs_stats["diffarea"], perm_stats[, "diffarea"]) >=
                                obs_stats["diffarea"]), 3)
        }
        obs_stats
    },
    ref = ref_ds,
    calculateStats = calculateStats,
    column = column,
    kmin = kmin,
    kfrac = kfrac,
    subsampleSize = subsampleSize,
    df = df,
    permutationPvalues = permutationPvalues,
    nPermutations = nPermutations,
    cpus4permut = cpus4permut,
    BPPARAM = bpparam)

    ds_res <- data.frame(data.table::rbindlist(lapply(ds_res_l, function(x)data.table::data.table(t(x))), idcol = "dataset"))
    ds_res

}
