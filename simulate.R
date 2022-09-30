#!/usr/bin/env Rscript
## Usage:
## simulate.R AZ_Prescott 1 trend 2 1 2000 
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
    data_fn <- 'AZ_Litchfield'
    season <- 2
    trend <- 'trend'
    config_model <- 5
    config_holdout <- 0
    niter <- 20000
}

## CAUTION: this uses about 88 GB memory for 5x20000 cases I have been doing.


library(nimble)

source("functions.R")

out <- combine_chains(data_fn, season, trend, config_model, config_holdout, niter, nburnin = 0) # 10000/10)
samples <- out$samples

source(paste0('config-model-', config_model, '.R'))
source('config-dirs.R')

source("load_data.R")
load(file.path(data_dir, 'basis_matrices.Rda'))  # for X1, X2

X1 <- X1[[season]]
X2 <- X2[[season]]

pD <- list(
    get_values(samples, 'lpD', 'beta1_pD', 'beta2_pD', X1, X2, int_element = 1, beta_element = NULL, doy = doy),
    get_values(samples, 'lpD', 'beta1_pD', 'beta2_pD', X1, X2, int_element = 2, beta_element = NULL, doy = doy),
    get_values(samples, 'lpD', 'beta1_pD', 'beta2_pD', X1, X2, int_element = 3, beta_element = NULL, doy = doy)
)
pDW <- get_values(samples, 'lpDW', 'beta1_pDW', 'beta2_pDW', X1, X2, int_element = NULL, beta_element = NULL, doy = doy)
pW <- list(
    get_values(samples, 'lpW', 'beta1_pW', 'beta2_pW', X1, X2, int_element = 1, beta_element = 1, doy = doy),
    get_values(samples, 'lpW', 'beta1_pW', 'beta2_pW', X1, X2, int_element = 2, beta_element = 2, doy = doy)
)
pWW <- list(
    get_values(samples, 'lpWW', 'beta1_pWW', 'beta2_pWW', X1, X2, int_element = 1, beta_element = 1, doy = doy),
    get_values(samples, 'lpWW', 'beta1_pWW', 'beta2_pWW', X1, X2, int_element = 2, beta_element = 2, doy = doy)
)
pWD <- list(
    get_values(samples, 'lpWD', 'beta1_pWD', 'beta2_pWD', X1, X2, int_element = 1, beta_element = 1, doy = doy),
    get_values(samples, 'lpWD', 'beta1_pWD', 'beta2_pWD', X1, X2, int_element = 2, beta_element = 2, doy = doy)
)


Pmat <- create_Pmat(pD, pDW, pW, pWW, pWD)

wh <- which(Pmat < 0)

if(length(wh)) {
    warning("some Pmat less than zero: ", min(Pmat[wh]))
## why necessary?
    Pmat[wh]  <- 0
}

pi <- list(
    get_values(samples, 'eta', 'beta1_eta', 'beta2_eta', X1, X2, int_element = 1, beta_element = 1, doy = doy),
    get_values(samples, 'eta', 'beta1_eta', 'beta2_eta', X1, X2, int_element = 2, beta_element = 2, doy = doy),
    get_values(samples, 'eta', 'beta1_eta', 'beta2_eta', X1, X2, int_element = 3, beta_element = 3, doy = doy)
)
scale <- list(
    get_values(samples, 'alpha', 'beta1_alpha', 'beta2_alpha', X1,  X2, int_element = 1, beta_element = 1, link = 'exp', doy = doy),
    get_values(samples, 'alpha', 'beta1_alpha', 'beta2_alpha', X1,  X2, int_element = 2, beta_element = 2, link = 'exp', doy = doy),
    get_values(samples, 'alpha', 'beta1_alpha', 'beta2_alpha', X1,  X2, int_element = 3, beta_element = 3, link = 'exp', doy = doy)
)
if(dens_type == 2) link <- 'exp' else link <- 'identity'
shape <- list(
    get_values(samples, 'gamma', 'beta1_gamma', 'beta2_gamma', X1, X2, int_element = 1, beta_element = 1, link = link, doy = doy),
    get_values(samples, 'gamma', 'beta1_gamma', 'beta2_gamma', X1, X2, int_element = 2, beta_element = 2, link = link, doy = doy),
    get_values(samples, 'gamma', 'beta1_gamma', 'beta2_gamma', X1, X2, int_element = 3, beta_element = 3, link = link, doy = doy)
)

P_zero <- samples[ , grep("P_zero", colnames(samples))]

if(dens_type == 3) {
    kappa <- colMeans(samples[ , grep("kappa", colnames(samples))])
} else kappa <- rep(1, W+1)

rSim <- simulate_vec_rep(P_zero, Pmat, 
                 pi, scale, shape, 
                 kappa, nT = nT, n = n, round = rounded, seed = 1)

output_fn <- paste('sims', data_fn, season,
                   ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
output_fn <- file.path(output_dir, paste0(output_fn, '.Rda'))

save(rSim, r, file = output_fn)

if(FALSE) {
    ## think more about GPD tail
    rSim[rSim > max(r)]  <- max(r)
    
    ff <- function(x) mean(x == 0)
    image.plot(1:118,1:92,apply(rSim,c(1,2),ff),zlim=c(0,1))
}

## tail analysis code
if(FALSE) {
    season <- 4
    type <- 2
    data <- "PA_station" # "AZ_Prescott" # "PA_station"
    load(paste0('output/sims-', data, '-', season, '-trend-', type, '-1.Rda'))
    rSim[rSim > max(r)]  <- 3*max(r)
#    it <- 300
#    sum(rSim[,,it] > max(r))
    qqplot(r,rSim[,,])
    abline(0,1)
}
