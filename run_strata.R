#!/usr/bin/env Rscript
## Usage:
## run_strata.R AZ_Prescott 1 trend 2 0 2000 1
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) {
    data_fn <- args[1]
    season <- as.numeric(args[2])
    trend <- args[3]
    config_model <- as.numeric(args[4])
    config_holdout <- as.numeric(args[5])
    niter <- as.numeric(args[6])
    replicate <- as.numeric(args[7])
} else cat("No config values found.\n")

if(FALSE) {
    ## run_strata.R AZ_Prescott 1 trend 2 0 2000 1
    data_fn <- 'CA_Quincy'
    season <- 5
    trend <- 'trend'
    config_model <- 3
    config_holdout <- 0
    niter <- 20000
    replicate <- 1
}

source(paste0('config-model-', config_model, '.R'))
source(paste0('config-holdout-', config_holdout, '.R'))
source('config-dirs.R')

cat("Configuring strata run with:\n")
cat("  Dataset: ", data_fn, "\n")
cat("  Season: ", season, "\n")
cat("  Trend: ", trend, "\n")
cat("  Model config: ", config_model, "\n")
cat("  Holdout config: ", config_holdout, "\n")
cat("  niter: ", niter, "\n")
cat("  replicate(seed): ", replicate, "\n")
cat("  output dir: ", output_dir, "\n")

source('fit.R')
