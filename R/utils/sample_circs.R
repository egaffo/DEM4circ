#' function to sample the desired fraction of circRNAs from each quantile
#'
sample_circs <- function(x, N, fracDE) {
  ## x: the set to sample from
  ## N: the number of elements to sample from the set
  ## fracDE: the fraction of sampled elements deemed to be differentially expressed

  if(N > length(x)){
    warning(paste("The number of elements to be sampled",
                  "is larger than the set:", N, ">", length(x), "\n",
                  "Will sample", length(x), "elements."))
    N <- length(x)
  }

  circs_to_sim <- sample(x, size = N)
  nonnull_circs <- sample(circs_to_sim, size = ceiling(length(circs_to_sim) * fracDE))
  null_circs <- circs_to_sim[! circs_to_sim %in% nonnull_circs]

  list(circs_to_sim = circs_to_sim,
       nonnull_circs = nonnull_circs,
       null_circs = null_circs)
}
