#' SPsimSeq simulation
#'
SPsimSeq_simulate <-
  function(params_list, ...) {

    library(SPsimSeq)
    seed.num <- 2021
    set.seed(seed.num)

    counts_mt <- params_list$counts_mt
    reps_per_group <- params_list$reps_per_group
    Ncircrnas <- params_list$Ncircrnas
    fracDECs <- params_list$fracDECs
    group_df <- params_list$group_df
    model_zeros <- params_list$model_zeros
    n.sims <- params_list$n.sims
    cand_DE_genes <- NULL
    if(!is.null(params_list$cand_DE_genes)) {
      cand_DE_genes <- list(nonnull.genes = params_list$cand_DE_genes,
                            null.genes = rownames(counts_mt)[!rownames(counts_mt) %in%
                                                               params_list$cand_DE_genes])
    }
    w <- params_list$w

    sim_pars <- paste(n.sims, "data sets of",
                      Ncircrnas, "circRNAs per",
                      2 * reps_per_group, "samples, modelling zero probability",
                      model_zeros)

    start_time <- Sys.time()
    message(paste(start_time, "start SPsimSeq to simulate",
                  sim_pars))

    simdatasets <-
      SPsimSeq::SPsimSeq(n.sim = n.sims,
                         s.data = counts_mt,
                         n.genes = Ncircrnas,
                         tot.samples = 2 * reps_per_group,
                         pDE = fracDECs,
                         cand.DE.genes = cand_DE_genes,
                         group = as.integer(factor(group_df[colnames(counts_mt), ])),
                         group.config = c(.5, .5), # fixed to equal group sizes
                         model.zero.prob = model_zeros,
                         batch = rep(1, ncol(counts_mt)),
                         batch.config = 1,
                         w = w,
                         result.format = "list",
                         return.details = TRUE,
                         ...)

    runtime <- difftime(Sys.time(), start_time, units = "s")
    message(paste("SPsimSeq took", runtime, "seconds"))

    list(Datasets = simdatasets,
         runtime = runtime)
  }
