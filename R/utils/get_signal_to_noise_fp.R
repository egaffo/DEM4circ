get_signal_to_noise <- function(x) {

  mat <- x@BenchDesign@data@data$cntdat

  cpms <- edgeR::cpm(mat, lib.size = colSums(mat) * edgeR::calcNormFactors(mat))

  stats_dt <-
    data.table(gene_id = rownames(cpms),
               fraczeros = apply(mat, 1, function(x) mean(x == 0)),
               avecpm = apply(cpms, 1, mean),
               varcpm = apply(cpms, 1, var),
               cvcpm = apply(cpms, 1, sd) / apply(cpms, 1, mean))

  ## get DE of each method
  alpha_thr <- 0.05

  adj_pv <- assays(x)[["adj_pv"]]
  rownames(adj_pv) <- rownames(mat)

  meth_stats_dt <-
    merge(melt(data.table(adj_pv, keep.rownames = "gene_id"),
               id.vars = "gene_id",
               value.name = "padj",
               variable.name = "Method")[, is_de := padj <= alpha_thr],
          stats_dt, by = "gene_id")

  ## use only methods that detected >= 5 DECs
  methods <- as.character(meth_stats_dt[is_de == T, .N, by = .(Method)][N >= 5, Method])

  if(length(methods) > 0){
    ## signal-to-noise statistic
    # (mean(DECs) - mean(notDECs)) / (sd(DECs) + sd(notDECs))

    dcast(dcast(melt(meth_stats_dt[Method %in% methods,
                                   .(`M_Fraction zeros` = mean(fraczeros),
                                     # `M_log2(average CPM)` = log2(mean(avecpm)),
                                     # `M_log2(variance CPM)` = log2(mean(varcpm)),
                                     `M_average CPM` = mean(avecpm),
                                     `M_variance CPM` = mean(varcpm),
                                     `M_CV(CPM)` = mean(cvcpm),
                                     `SD_Fraction zeros` = sd(fraczeros),
                                     # `SD_log2(average CPM)` = log2(sd(avecpm)),
                                     # `SD_log2(variance CPM)` = log2(sd(varcpm)),
                                     `SD_average CPM` = sd(avecpm),
                                     `SD_variance CPM` = sd(varcpm),
                                     `SD_CV(CPM)` = sd(cvcpm)),
                                   by = .(Method, is_de)],
                     id.vars = c("Method", "is_de"),
                     value.name = "Stat_val",
                     variable.name = "Stat"),
                formula = Method + Stat ~ is_de,
                value.var = "Stat_val")[, c("Stat", "Val") := tstrsplit(Stat, "_")],
          Method + Val ~ Stat,
          value.var = c("FALSE", "TRUE"))[, .(StN = (TRUE_M - FALSE_M) / (TRUE_SD + FALSE_SD)),
                                          by = .(Method, Val)]
  }
}
