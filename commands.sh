#!/bin/bash

## N.B. run the "Rscript"s from the main directory to use renv

#---- analyse the circRNA expression of the source data sets ----#
Rscript -e "rmarkdown::render(input = 'R/circexp_characteristics/spliced_read_fractions.Rmd', \
output_dir = 'R/circexp_characteristics', params = list(datadir = '../../data/PRJCA000751', \
outdir = 'spliced_read_fractions'))"

Rscript -e "rmarkdown::render(input = 'R/circexp_characteristics/source_ds_characteristics.Rmd', \
output_dir = 'R/circexp_characteristics')"

Rscript -e "rmarkdown::render(input = 'R/circexp_characteristics/bjr_filters.Rmd', \
output_dir = 'R/circexp_characteristics', params = list(datadir = '../../data', \
outdir = 'bjr_filters'))"

Rscript -e "rmarkdown::render(input = 'R/circexp_characteristics/gof.Rmd', \
output_dir = 'R/circexp_characteristics', params = list(outdir = 'gof', datadir = '../../data', \
filtdatadir = 'bjr_filters'))"

#---- simulate the data sets ----#
Rscript -e "rmarkdown::render(input = 'R/simulations/semiparametric/datasets_simulation_DM1.Rmd', \
params = list(n_sims = 30, max_circrrnas = 10000, dec_fraction = 0.1), \
output_dir = 'R/simulations/semiparametric')"

Rscript -e "rmarkdown::render(input = 'R/simulations/semiparametric/datasets_simulation_IDC.Rmd', \
params = list(n_sims = 30, max_circrrnas = 10000, dec_fraction = 0.1), \
output_dir = 'R/simulations/semiparametric')"

Rscript -e "rmarkdown::render(input = 'R/simulations/semiparametric/datasets_simulation_IPF.Rmd', \
params = list(n_sims = 30, max_circrrnas = 10000, dec_fraction = 0.1), \
output_dir = 'R/simulations/semiparametric')"

Rscript -e "rmarkdown::render(input = 'R/simulations/semiparametric/datasets_simulation_MS.Rmd',  \
params = list(n_sims = 30, max_circrrnas = 10000, dec_fraction = 0.1), \
output_dir = 'R/simulations/semiparametric')"

---- assess the quality of the simulated data sets ----#
Rscript -e "rmarkdown::render(input = 'R/simulations/simdata_qual/assess_simdata_qual.Rmd', \
params = list(simulated_datasets_qs = '../semiparametric/simdata/DM1_trimmed_simulated_datasets.qs', \
trimmed_source_dataset.qs = '../semiparametric/simdata/DM1_trimmed_source_dataset.qs', \
de_ds_stats_outfile = 'simdata_qual/DM1_de_ds_stats.csv', \
null_ds_stats_outfile = 'simdata_qual/DM1_null_ds_stats.csv'), \
output_file = 'R/simulations/simdata_qual/assess_simdata_qual_DM1.html', output_dir = 'R/simulations/simdata_qual')"

Rscript -e "rmarkdown::render(input = 'R/simulations/simdata_qual/assess_simdata_qual.Rmd', \
params = list(simulated_datasets_qs = '../semiparametric/simdata/IDC_trimmed_simulated_datasets.qs', \
trimmed_source_dataset.qs = '../semiparametric/simdata/IDC_trimmed_source_dataset.qs', \
de_ds_stats_outfile = 'simdata_qual/IDC_de_ds_stats.csv', \
null_ds_stats_outfile = 'simdata_qual/IDC_null_ds_stats.csv'), \
output_file = 'R/simulations/simdata_qual/assess_simdata_qual_IDC.html', output_dir = 'R/simulations/simdata_qual')"

Rscript -e "rmarkdown::render(input = 'R/simulations/simdata_qual/assess_simdata_qual.Rmd', \
params = list(simulated_datasets_qs = '../semiparametric/simdata/IPF_trimmed_simulated_datasets.qs', \
trimmed_source_dataset.qs = '../semiparametric/simdata/IPF_trimmed_source_dataset.qs', \
de_ds_stats_outfile = 'simdata_qual/IPF_de_ds_stats.csv', \
null_ds_stats_outfile = 'simdata_qual/IPF_null_ds_stats.csv'), \
output_file = 'R/simulations/simdata_qual/assess_simdata_qual_IPF.html', output_dir = 'R/simulations/simdata_qual')"

Rscript -e "rmarkdown::render(input = 'R/simulations/simdata_qual/assess_simdata_qual.Rmd', \
params = list(simulated_datasets_qs = '../semiparametric/simdata/MS_trimmed_simulated_datasets.qs',  \
trimmed_source_dataset.qs = '../semiparametric/simdata/MS_trimmed_source_dataset.qs', \
de_ds_stats_outfile = 'simdata_qual/MS_de_ds_stats.csv', \
null_ds_stats_outfile = 'simdata_qual/MS_null_ds_stats.csv'), \
output_file = 'R/simulations/simdata_qual/assess_simdata_qual_MS.html',  output_dir = 'R/simulations/simdata_qual')"

#---- filter the simulated data sets ----#
Rscript -e "rmarkdown::render(input = 'R/benchmark/filter_simdata.Rmd', \
params = list(outdir = 'filtered_simdata/DM1', \
simulated_datasets_qs = '../simulations/semiparametric/simdata/DM1_trimmed_simulated_datasets.qs'), \
output_file = 'R/benchmark/filter_simdata_DM1.html', output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/filter_simdata.Rmd', \
params = list(outdir = 'filtered_simdata/IDC', \
simulated_datasets_qs = '../simulations/semiparametric/simdata/IDC_trimmed_simulated_datasets.qs'), \
output_file = 'R/benchmark/filter_simdata_IDC.html', output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/filter_simdata.Rmd', \
params = list(outdir = 'filtered_simdata/IPF', \
simulated_datasets_qs = '../simulations/semiparametric/simdata/IPF_trimmed_simulated_datasets.qs'), \
output_file = 'R/benchmark/filter_simdata_IPF.html', output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/filter_simdata.Rmd', \
params = list(outdir = 'filtered_simdata/MS',  \
simulated_datasets_qs = '../simulations/semiparametric/simdata/MS_trimmed_simulated_datasets.qs'), \
output_file = 'R/benchmark/filter_simdata_MS.html',  output_dir = 'R/benchmark')"

# ---- run differential expression methods and compute performance metrics ----#
Rscript -e "rmarkdown::render(input = 'R/benchmark/benchmark.Rmd', \
params = list(outdir = 'benchmark_results/DM1', \
inputData = 'filtered_simdata/DM1/datasetList.qs'), \
output_file = 'R/benchmark/benchmark_DM1.html', output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/benchmark.Rmd', \
params = list(outdir = 'benchmark_results/IDC', \
inputData = 'filtered_simdata/IDC/datasetList.qs'), \
output_file = 'R/benchmark/benchmark_IDC.html', output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/benchmark.Rmd', \
params = list(outdir = 'benchmark_results/IPF', \
inputData = 'filtered_simdata/IPF/datasetList.qs'), \
output_file = 'R/benchmark/benchmark_IPF.html', output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/benchmark.Rmd', \
params = list(outdir = 'benchmark_results/MS', \
inputData = 'filtered_simdata/MS/datasetList.qs'), \
output_file = 'R/benchmark/benchmark_MS.html',  output_dir = 'R/benchmark')"

#---- evaluate performance metrics ----#
rm -r R/benchmark/overall_evaluations
Rscript -e "rmarkdown::render(input = 'R/benchmark/overall_evaluations.Rmd', output_dir = 'R/benchmark')"
Rscript -e "rmarkdown::render(input = 'R/benchmark/CAT.Rmd', output_dir = 'R/benchmark')"
Rscript -e "rmarkdown::render(input = 'R/benchmark/jaccard.Rmd', output_dir = 'R/benchmark')"

---- simulate the data sets with nonparametric approach ----#
Rscript R/simulations/nonparametric/simSeq_PMID35078526.R "data/PMID35078526_rawCIRI2output.tsv" "R/simulations/nonparametric"

#---- assess the quality of the simulated data sets ----#
Rscript -e "rmarkdown::render(input = 'R/simulations/simdata_qual/assess_simdata_qual.Rmd', \
params = list(simulated_datasets_qs = '../nonparametric/simdata/nonp_simdata_PC.qs', \
trimmed_source_dataset.qs = '../nonparametric/simdata/PC_trimmed_source_dataset.qs', \
de_ds_stats_outfile = 'simdata_qual/PC_de_ds_stats.csv', \
null_ds_stats_outfile = 'simdata_qual/PC_null_ds_stats.csv'), \
output_file = 'R/simulations/simdata_qual/assess_simdata_qual_PC.html', output_dir = 'R/simulations/simdata_qual')"

#---- filter the simulated data sets ----#
Rscript -e "rmarkdown::render(input = 'R/benchmark/filter_simdata.Rmd', \
params = list(outdir = 'filtered_simdata/PC', \
simulated_datasets_qs = '../simulations/nonparametric/simdata/nonp_simdata_PC.qs'), \
output_file = 'R/benchmark/filter_simdata_PC.html', output_dir = 'R/benchmark')"

#---- run differential expression methods and compute performance metrics ----#
Rscript -e "rmarkdown::render(input = 'R/benchmark/benchmark.Rmd', \
params = list(outdir = 'benchmark_results/PC', \
inputData = 'filtered_simdata/PC/datasetList.qs'), \
output_file = 'R/benchmark/benchmark_PC.html', output_dir = 'R/benchmark')"

#---- evaluate performance metrics ----#
rm -r R/benchmark/PC_overall_evaluations
Rscript -e "rmarkdown::render(input = 'R/benchmark/PC_overall_evaluations.Rmd', \
output_dir = 'R/benchmark')"

Rscript -e "rmarkdown::render(input = 'R/benchmark/sp_vs_np_evaluations.Rmd', \
output_dir = 'R/benchmark')"
