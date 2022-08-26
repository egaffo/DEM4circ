args <- commandArgs(trailingOnly = T)
input <- args[1] #"PMID35078526_rawCIRI2output.tsv"
outdir <- args[2]

library(SimSeq)
library(data.table)
library(qs)
library(BiocParallel)
library(edgeR)
library(fdrtool)

nsims <- 30

outdir <- file.path(outdir, "simdata")
dir.create(outdir, showWarnings = F, recursive = T)

simdata_file <- file.path(outdir, "nonp_simdata_PC.qs")
# simulated_datasets_qs <- file.path(params$outdir, "PC_trimmed_simulated_datasets.qs")
source_dataset_qs <- file.path(outdir, "PC_trimmed_source_dataset.qs")

# setwd(utils::getSrcDirectory()[1])
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

counts_dt <- fread(input,
                   data.table = T,
                   showProgress = F,
                   drop = c("chr", "start", "end", "strand", "gene_id", "circRNA_type"))

## filter the original dataset
counts_dt_m <- melt(counts_dt,
                    id.vars = "circRNA_ID",
                    variable.name = "sample_id",
                    value.name = "BJR")
counts_dt_m[, c("ID", "Group") := tstrsplit(sample_id, "_"), by = sample_id]

## keep only circRNAs expressed in > 2 samples per condition
keep <- counts_dt_m[BJR > 0, .N,
                    by = .(circRNA_ID, Group)][N > 2 &
                                                   Group != "MPC", unique(circRNA_ID)]

## read sample groups
sampleTable <- data.frame(unique(counts_dt_m[Group != "MPC",
                                             .(sample_id, Group)][order(Group, sample_id)]),
                          row.names = "sample_id")

counts <-
    as.matrix(data.frame(counts_dt, check.names = F,
                         row.names = "circRNA_ID"))[keep, rownames(sampleTable)]
# dim(counts)

## check possible biases
dge <- DGEList(counts = counts[, rownames(sampleTable)], group = sampleTable$Group)
dge <- calcNormFactors(dge, method = "TMM")
cpms <- edgeR::cpm(dge, log = T, prior.count = 1)

pcs <- prcomp(x = t(cpms), scale. = F, center = T)
df <- cbind(sampleTable, pcs$x)

# library(ggplot2)
# ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
#     geom_point()
#
# plot(pcs$x[, 2])

df$Batch <- ifelse(df$PC2 > 0, "B1", "B2")
sampleTable$Batch <- df[rownames(sampleTable), "Batch"]
## get batch-group sample sizes
table(sampleTable[, c("Batch", "Group")])

## check group separation after correcting by PC2
# library(limma)
# mat <- removeBatchEffect(x = cpms,
#                          batch = df$Batch,
#                          # covariates = df$PC2,
#                          design = model.matrix(~ Group, data = sampleTable))
#
# pcs <- prcomp(x = t(mat), scale. = F, center = T)
# df <- cbind(sampleTable, pcs$x)
# ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
#     geom_point()

## keep only samples from the largest batch
st <- sampleTable[sampleTable$Batch == "B2", "Group", drop = F]
ct <- counts[, rownames(st)]
## keep circRNAs detected in at least 3 samples
ct <- ct[rowSums(ct > 0) > 2, ]
# ## keep circRNAs detected in at least half samples
# ct <- ct[rowSums(ct > 0) >= floor(dim(ct)[2] / 2), ]
dim(ct)

## save source dataset
orig_mat_dge <- edgeR::DGEList(counts = ct,
                               group = st[colnames(ct), "condition"])
qsave(x = orig_mat_dge,
      file = source_dataset_qs,
      nthreads = multicoreWorkers(),
      preset = "fast")

## parameters for the simulations
# n_genes <- nrow(ct) # simulate as many circRNAs as in the source data set
# n_diff <- floor(n_genes * .1) # 10% DECs

nf <- calcNormFactors(ct, method = "TMM")
sort.list <- SortData(counts = ct,
                      treatment = st$Group,
                      replic = NULL,
                      sort.method = "unpaired",
                      norm.factors = nf)

sort.list$probs <-
    CalcPvalWilcox(counts = sort.list$counts,
                   treatment = sort.list$treatment,
                   sort.method = "unpaired",
                   sorted = TRUE,
                   norm.factors = sort.list$norm.factors,
                   exact = FALSE)

sort.list$weights <- 1 - fdrtool(sort.list$probs,
                                 statistic = "pvalue",
                                 plot = FALSE,
                                 verbose = FALSE)$lfdr

#' 'set_name' must be on the form 'NSS_ID', where SS is the number of samples
#' of each group and ID is the simulation identifier. F.i. N03_01 will be the
#' name of the simulation 01 having 3 samples per group
simData <- function(set_name, sort.list) {

    samplesPerGroup <- as.integer(sub("N", "", unlist(strsplit(set_name, "_"))[1]))

    simdata <-
        SimSeq::SimData(counts = sort.list$counts,
                        treatment = sort.list$treatment,
                        replic = sort.list$replic,
                        sort.method = "unpaired",
                        k.ind = samplesPerGroup,
                        n.genes = nrow(sort.list$counts), # simulate as many circRNAs as in the source data set,
                        n.diff = ifelse(is.null(sort.list$n.diff),
                                        floor(nrow(sort.list$counts) * .1), # 10% DECs
                                        sort.list$n.diff),
                        norm.factors = sort.list$norm.factors,
                        # samp.independent = FALSE,
                        # genes.select = NULL,
                        # genes.diff = NULL,
                        # switch.trt = F,
                        # probs = NULL,
                        weights = sort.list$weights,
                        # exact = FALSE,
                        power = 1)

    simdata$set_name <- set_name

    ## add fields to make it compatible with SPsimSeq objects
    colnames(simdata$counts) <- paste0("Sample_", 1:ncol(simdata$counts))

    simdata$rowData <- data.frame(DE.ind = simdata$DE.ind,
                                  source.ID = rownames(simdata$counts),
                                  row.names = rownames(simdata$counts))

    simdata$coldata <- data.frame(Batch = 1,
                                  Group = factor(as.integer(simdata$treatment)+1), # will relevel factor from {0,1} into {1,2}
                                  sim.Lib.size = colSums(simdata$counts))
    rownames(simdata$coldata) <- colnames(simdata$counts)

    simdata$colData <- simdata$coldata

    simdata
}

## simulate data sets
set.seed(123)

# samplesPerGroup <- 5 # number of samples in each group to simulate

set_names <- c(paste0(paste0("N", formatC(3, width = 2, flag = "0"), "_"),
                      formatC(seq_len(nsims), width = 2, flag = "0")),
               paste0(paste0("N", formatC(5, width = 2, flag = "0"), "_"),
                      formatC(seq_len(nsims), width = 2, flag = "0")),
               paste0(paste0("N", formatC(10, width = 2, flag = "0"), "_"),
                      formatC(seq_len(nsims), width = 2, flag = "0")))

## DE
sort.list$n.diff <- NULL
de_simdataL <- sapply(X = set_names,
                      FUN = simData,
                      sort.list = sort.list,
                      USE.NAMES = T,
                      simplify = F)

names(de_simdataL) <- paste0(names(de_simdataL), "_bulk")

## not DE
sort.list$n.diff <- 0
mock_simdataL <- sapply(X = set_names,
                        FUN = simData,
                        sort.list = sort.list,
                        USE.NAMES = T,
                        simplify = F)

names(mock_simdataL) <- paste0(names(mock_simdataL), "_bulk_mock")

sim_ds_list <- list("N03_bulk" = list(Datasets = list(sim.data.list = de_simdataL[grepl("N03", names(de_simdataL))]), runtime = NA),
                    "N05_bulk" = list(Datasets = list(sim.data.list = de_simdataL[grepl("N05", names(de_simdataL))]), runtime = NA),
                    "N10_bulk" = list(Datasets = list(sim.data.list = de_simdataL[grepl("N10", names(de_simdataL))]), runtime = NA),
                    "N03_bulk_mock" = list(Datasets = list(sim.data.list = mock_simdataL[grepl("N03", names(mock_simdataL))]), runtime = NA),
                    "N05_bulk_mock" = list(Datasets = list(sim.data.list = mock_simdataL[grepl("N05", names(mock_simdataL))]), runtime = NA),
                    "N10_bulk_mock" = list(Datasets = list(sim.data.list = mock_simdataL[grepl("N10", names(mock_simdataL))]), runtime = NA))

## save simulations
qsave(x = sim_ds_list, file = simdata_file, nthreads = multicoreWorkers())
