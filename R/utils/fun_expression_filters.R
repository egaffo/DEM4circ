#' Fraction of samples filter
#'
#' Each circRNA must be detected in at least the specified fraction of samples
#'
#' @param x the matrix to be filtered. It must have a \code{circ_id} column
#' @param rthr the minimum BJR count for a circRNA to be kept
#' @param frac the fraction of the samples in which each circRNA must be
#'   detected. Default is that each circRNA must be detected in al least half
#'   the samples
#'
#' @return the matrix bearing circRNAs detected in at least a fraction of the
#'   samples
#' @export
#'
#' @examples
#' @import data.table
fraction_the_samples_filter <-
  function(x, rthr = 0, frac = .5) {
    xmat <- as.matrix(data.frame(x, row.names = "circ_id"))
    x[circ_id %in%
        rownames(xmat)[rowSums(xmat >= rthr) >=
                         ceiling(length(colnames(xmat)) * frac)]]
  }

#' Smallest condition size filter
#'
#' Each circRNA must be detected in at least as many samples as in the smallest
#' condition group
#'
#' @param x the matrix to be filtered. It must have a \code{circ_id} column
#' @param rthr the minimum BJR count for a circRNA to be kept
#' @param cond the sample metadata table. Must have \code{sample_id} and
#'   \code{condition} columns
#'
#' @return a matrix bearing circRNAs expressed in as many samples as the size of
#'   the smallest condition group
#' @export
#'
#' @examples
#' @import data.table
smallest_group_filter <-
  function(x, rthr, cond) {

    smallest_group_n <- cond[, .N, by = condition][, min(N)]

    circ_ids_to_keep <-
      melt(x,
           id.vars = "circ_id",
           variable.name = "sample_id",
           value.name = "BJR")[BJR >= rthr,
                               .N,
                               by = .(circ_id)][N >= smallest_group_n,
                                                unique(circ_id)]
    x[circ_id %in% circ_ids_to_keep]
  }

#' Condition consistency filter
#'
#' CircRNAs must be detected in all the samples of at least one condition
#'
#' @param x the matrix to be filtered. It must have a \code{circ_id} column
#' @param rthr the minimum BJR count for a circRNA to be kept
#' @param cond the sample metadata table. Must have \code{sample_id} and
#'   \code{condition} columns
#'
#' @return the matrix bearing circRNAs expressed in all samples of at least one
#'   condition
#' @export
#'
#' @examples
#' @import data.table
condition_consistent_filter <-
  function(x, rthr, cond) {

    circ_ids_to_keep <-
      merge(cond[, .(samples_per_cond = .N),
                 by = condition],
            merge(melt(x,
                       id.vars = "circ_id",
                       variable.name = "sample_id",
                       value.name = "BJR"),
                  cond,
                  by = "sample_id")[BJR >= rthr
                  ][, .N,
                    by = .(condition,
                           circ_id)],
            by = "condition")[samples_per_cond == N,
                              unique(circ_id)]

    x[circ_id %in% circ_ids_to_keep]
  }


#' Keep circRNAs expressed in >= N samples
#'
#' @param x
#' @param rthr
#' @param Nsamples
#'
#' @return
#' @export
#'
#' @examples
n_samples_filter <-
  function (x, rthr, Nsamples = 1) {
    xmat <- as.matrix(data.frame(x, row.names = "circ_id"))
    x[circ_id %in% rownames(xmat)[rowSums(xmat >= rthr) >= Nsamples]]
  }
