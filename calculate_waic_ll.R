#!/usr/bin/env Rscript
## Usage:
## calculate_waic.R AZ_Prescott 1 trend 2 1 2000 
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) {
    data_fn <- args[1]
    season <- as.numeric(args[2])
    trend <- args[3]
    config_model <- as.numeric(args[4])
    config_holdout <- as.numeric(args[5])
    niter <- as.numeric(args[6])
} else cat("No config values found.\n")

if(FALSE) {
    data_fn <- 'CA_Berkeley'
    season <- 1
    trend <- 'trend'
    config_model <- 5
    config_holdout <- 0
    niter <- 20000
}

library(nimble)

source('functions.R')
source(paste0('config-model-', config_model, '.R'))
source('config-dirs.R')

source('nimble_functions.R')

output_pattern <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, sep = '-')
## Load in basic objects saved in output file using first output file.
load(file.path(output_dir, paste0(output_pattern, '-1.Rda')))

nburnin <- (niter/thin)*burnin_waic_frac

rain_model <- nimbleModel(rain_code, rain_constants, rain_data, rain_inits, dimensions = rain_dims)

rain_model_compiled <- compileNimble(rain_model)

samples <- combine_chains(data_fn, season, trend, config_model, config_holdout, niter, nburnin = nburnin)$samples 

if(nrow(samples) != 7500)
    warning('expecting 7500 samples')

waic <- calculateWAIC(samples, rain_model_compiled)
print(waic)

ll <- combine_ll(data_fn, season, trend, config_model, config_holdout, niter, nburnin = nburnin)

output_fn <- paste('waic', data_fn, season,
                   ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
output_fn <- file.path(output_dir, paste0(output_fn, '.Rda'))

save(waic, ll, file = output_fn)


