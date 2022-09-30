sets <- 1:8

trend <- 'trend'
config_model <- 5
config_holdout <- 0
niter <- 20000

stns <- c(rep('AZ_Winslow', 4), rep('CA_Berkeley', 3), 'CA_Quincy')
seas <- c(1:4, 1,2,4,5)

source(paste0('config-model-', config_model, '.R'))
source(paste0('config-holdout-', config_holdout, '.R'))
source('config-dirs.R')
source('functions.R')

resultsDir <- file.path('..', 'paper1', 'results')

## Stan approach to Rhat - columns are chains, rows are iterations
## sims <- matrix(rnorm(500), nrow = 100, ncol = 5)
## Rhat(sims)
## ess_bulk(sims)
## ess_tail(sims)

library(rstan)
library(coda)

RhatMat <- function(x, nreps = 5) {
    Rhat(matrix(x, ncol = nreps))
}
essMat <- function(x, nreps = 5) {
    ess_bulk(matrix(x, ncol = nreps))
}
esstMat <- function(x, nreps = 5) {
    ess_tail(matrix(x, ncol = nreps))
}

RhatOldMat <- function(x, nreps = 5) {
    mlist <- as.mcmc.list(lapply(seq_len(nreps), function(i)
                                 as.mcmc(x[chains[,i]])))
    gelman.diag(mlist)$psrf[1,1]
}

essOldMat <- function(x, nreps = 5) {
    mlist <- as.mcmc.list(lapply(seq_len(nreps), function(i)
                                 as.mcmc(x[chains[,i]])))
    effectiveSize(mlist)
}

cmPerInch <- 2.54
dryCut <- 0.3  # cm below which a day is considered 'dry'

calc_rain <- function(x, cut = dryCut / cmPerInch, na.rm = FALSE) mean(x>cut, na.rm = na.rm)

## Load in for a single dataset to get some of the configuration values
data_fn <- 'CA_Berkeley'
season <- 1
fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
load(file.path(output_dir, paste0(fn, '.Rda')))

burnin <- burnin_waic_frac*niter/thin
niter_per_chain <- niter/thin

nreps <- dim(rSim)[[3]] / niter_per_chain

starts <- seq(1, (nreps-1)*niter_per_chain+1, by = niter_per_chain) + burnin
ends <- starts + niter_per_chain - burnin - 1

keep <- starts[1]:ends[1]
for(i in 2:length(starts))
    keep <- c(keep, starts[i]:ends[i])

chains <- matrix(1:length(keep), ncol = nreps)



yrs <- 1900:2021 
timeInts <- list(1920:2021, 1950:2021, 1980:2021)
rg <- match(timeInts[[1]], yrs)

years <- yrs[rg]

qs <- c(.5, .95, .99)

trans <- list()

for(i in sets) {
    print(c(i, date()))
    data_fn <- stns[i]
    season <- seas[i]

    fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
    load(file.path(output_dir, paste0(fn, '.Rda')))

    rSim <- rSim[rg, , keep]
    trans[[i]] <- apply(rSim, c(1,3), calc_trans, 0.3/2.54)
}

summ_senDry <- summ_senWet <- matrix(0, length(seas), 3)
summOld_senDry <- summOld_senWet <- matrix(0, length(seas), 2)

rhat_qu <- ess_qu <- array(0, c(length(seas), 12, 3))
rhatOld_qu <- essOld_qu <- array(0, c(length(seas), 12, 3))


for(i in sets) {
    data_fn <- stns[i]
    season <- seas[i]

    fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
    load(file.path(output_dir, paste0(fn, '.Rda')))

    rSim <- rSim[rg, , keep]
    nIts <- dim(rSim)[3]

    rainTime <- apply(rSim, c(1,3), calc_rain, 0.3/2.54)
    rainSeas <- apply(rSim, c(2,3), calc_rain, 0.3/2.54)
    rainTime1 <- apply(rSim, c(1,3), calc_rain, 1/2.54)
    rainSeas1 <- apply(rSim, c(2,3), calc_rain, 1/2.54)
    rainTime2 <- apply(rSim, c(1,3), calc_rain, 2/2.54)
    rainSeas2 <- apply(rSim, c(2,3), calc_rain, 2/2.54)

    rhat_rainTime <- apply(rainTime, 1, RhatMat)
    rhat_rainSeas <- apply(rainSeas, 1, RhatMat)
    rhat_rainTime1 <- apply(rainTime1, 1, RhatMat)
    rhat_rainSeas1 <- apply(rainSeas1, 1, RhatMat)
    rhat_rainTime2 <- apply(rainTime2, 1, RhatMat)
    rhat_rainSeas2 <- apply(rainSeas2, 1, RhatMat)

    ess_rainTime <- apply(rainTime, 1, essMat)
    ess_rainSeas <- apply(rainSeas, 1, essMat)
    ess_rainTime1 <- apply(rainTime1, 1, essMat)
    ess_rainSeas1 <- apply(rainSeas1, 1, essMat)
    ess_rainTime2 <- apply(rainTime2, 1, essMat)
    ess_rainSeas2 <- apply(rainSeas2, 1, essMat)

    rhatOld_rainTime <- apply(rainTime, 1, RhatOldMat)
    rhatOld_rainSeas <- apply(rainSeas, 1, RhatOldMat)
    rhatOld_rainTime1 <- apply(rainTime1, 1, RhatOldMat)
    rhatOld_rainSeas1 <- apply(rainSeas1, 1, RhatOldMat)
    rhatOld_rainTime2 <- apply(rainTime2, 1, RhatOldMat)
    rhatOld_rainSeas2 <- apply(rainSeas2, 1, RhatOldMat)

    essOld_rainTime <- apply(rainTime, 1, essOldMat)
    essOld_rainSeas <- apply(rainSeas, 1, essOldMat)
    essOld_rainTime1 <- apply(rainTime1, 1, essOldMat)
    essOld_rainSeas1 <- apply(rainSeas1, 1, essOldMat)
    essOld_rainTime2 <- apply(rainTime2, 1, essOldMat)
    essOld_rainSeas2 <- apply(rainSeas2, 1, essOldMat)
    
    ## transDW, transWW not useful as just one minus transDD and transWD

    ## NAs can occur when never in a given state (generally the wet state)
    ## set to one since most NAs are for transWD and transWW and we
    ## focus on transWD, which would generally be 1 if it had rained in cases where days with precip are rare
    trans[[i]][is.na(trans[[i]])] <- 1

    rhat_transDD <- apply(trans[[i]][1,,], 1, RhatMat)
    rhat_transDW <- apply(trans[[i]][3,,], 1, RhatMat)
    rhat_transWD <- apply(trans[[i]][2,,], 1, RhatMat)
    rhat_transWW <- apply(trans[[i]][4,,], 1, RhatMat)

    ess_transDD <- apply(trans[[i]][1,,], 1, essMat)
    ess_transDW <- apply(trans[[i]][3,,], 1, essMat)
    ess_transWD <- apply(trans[[i]][2,,], 1, essMat)
    ess_transWW <- apply(trans[[i]][4,,], 1, essMat)

    rhatOld_transDD <- apply(trans[[i]][1,,], 1, RhatOldMat)
    rhatOld_transDW <- apply(trans[[i]][3,,], 1, RhatOldMat)
    rhatOld_transWD <- apply(trans[[i]][2,,], 1, RhatOldMat)
    rhatOld_transWW <- apply(trans[[i]][4,,], 1, RhatOldMat)

    essOld_transDD <- apply(trans[[i]][1,,], 1, essOldMat)
    essOld_transDW <- apply(trans[[i]][3,,], 1, essOldMat)
    essOld_transWD <- apply(trans[[i]][2,,], 1, essOldMat)
    essOld_transWW <- apply(trans[[i]][4,,], 1, essOldMat)
    ## probably not useful as just 1 minus the other quantity

    ## dsl, wet spell magnitude
    fn <- paste('spells', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
    load(file.path(output_dir, paste0(fn, '.Rda')))

    dryLenSim <- dryLenSim[ , rg]  
    wetAvgSim <- wetAvgSim[ , rg]  

    rhat_dsl <- apply(dryLenSim, 2, RhatMat)
    rhat_int <- apply(wetAvgSim, 2, RhatMat)
   
    ess_dsl <- apply(dryLenSim, 2, essMat)
    ess_int <- apply(wetAvgSim, 2, essMat)

    rhatOld_dsl <- apply(dryLenSim, 2, RhatOldMat)
    rhatOld_int <- apply(wetAvgSim, 2, RhatOldMat)
   
    essOld_dsl <- apply(dryLenSim, 2, essOldMat)
    essOld_int <- apply(wetAvgSim, 2, essOldMat)

    ## summaries
    rhat_qu[i, 1, ] <-  quantile(rhat_rainTime, qs)
    rhat_qu[i, 2, ] <-  quantile(rhat_rainTime1, qs)
    rhat_qu[i, 3, ] <-  quantile(rhat_rainTime2, qs)
    rhat_qu[i, 4, ] <-  quantile(rhat_rainSeas, qs)
    rhat_qu[i, 5, ] <-  quantile(rhat_rainSeas1, qs)
    rhat_qu[i, 6, ] <-  quantile(rhat_rainSeas2, qs)
    rhat_qu[i, 7, ] <-  quantile(rhat_transDD, qs)
    rhat_qu[i, 8, ] <-  quantile(rhat_transDW, qs)
    rhat_qu[i, 9, ] <-  quantile(rhat_transWD, qs)
    rhat_qu[i, 10, ] <-  quantile(rhat_transWW, qs)
    rhat_qu[i, 11, ] <-  quantile(rhat_dsl, qs)
    rhat_qu[i, 12, ] <-  quantile(rhat_int, qs)
    
    ess_qu[i, 1, ] <-  quantile(ess_rainTime, 1-qs)
    ess_qu[i, 2, ] <-  quantile(ess_rainTime1, 1-qs)
    ess_qu[i, 3, ] <-  quantile(ess_rainTime2, 1-qs)
    ess_qu[i, 4, ] <-  quantile(ess_rainSeas, 1-qs)
    ess_qu[i, 5, ] <-  quantile(ess_rainSeas1, 1-qs)
    ess_qu[i, 6, ] <-  quantile(ess_rainSeas2, 1-qs)
    ess_qu[i, 7, ] <-  quantile(ess_transDD, 1-qs)
    ess_qu[i, 8, ] <-  quantile(ess_transDW, 1-qs)
    ess_qu[i, 9, ] <-  quantile(ess_transWD, 1-qs)
    ess_qu[i, 10, ] <-  quantile(ess_transWW, 1-qs)
    ess_qu[i, 11, ] <-  quantile(ess_dsl, 1-qs)
    ess_qu[i, 12, ] <-  quantile(ess_int, 1-qs)


    rhatOld_qu[i, 1, ] <-  quantile(rhatOld_rainTime, qs)
    rhatOld_qu[i, 2, ] <-  quantile(rhatOld_rainTime1, qs)
    rhatOld_qu[i, 3, ] <-  quantile(rhatOld_rainTime2, qs)
    rhatOld_qu[i, 4, ] <-  quantile(rhatOld_rainSeas, qs)
    rhatOld_qu[i, 5, ] <-  quantile(rhatOld_rainSeas1, qs)
    rhatOld_qu[i, 6, ] <-  quantile(rhatOld_rainSeas2, qs)
    rhatOld_qu[i, 7, ] <-  quantile(rhatOld_transDD, qs)
    rhatOld_qu[i, 8, ] <-  quantile(rhatOld_transDW, qs)
    rhatOld_qu[i, 9, ] <-  quantile(rhatOld_transWD, qs)
    rhatOld_qu[i, 10, ] <-  quantile(rhatOld_transWW, qs)
    rhatOld_qu[i, 11, ] <-  quantile(rhatOld_dsl, qs)
    rhatOld_qu[i, 12, ] <-  quantile(rhatOld_int, qs)
    
    essOld_qu[i, 1, ] <-  quantile(essOld_rainTime, 1-qs)
    essOld_qu[i, 2, ] <-  quantile(essOld_rainTime1, 1-qs)
    essOld_qu[i, 3, ] <-  quantile(essOld_rainTime2, 1-qs)
    essOld_qu[i, 4, ] <-  quantile(essOld_rainSeas, 1-qs)
    essOld_qu[i, 5, ] <-  quantile(essOld_rainSeas1, 1-qs)
    essOld_qu[i, 6, ] <-  quantile(essOld_rainSeas2, 1-qs)
    essOld_qu[i, 7, ] <-  quantile(essOld_transDD, 1-qs)
    essOld_qu[i, 8, ] <-  quantile(essOld_transDW, 1-qs)
    essOld_qu[i, 9, ] <-  quantile(essOld_transWD, 1-qs)
    essOld_qu[i, 10, ] <-  quantile(essOld_transWW, 1-qs)
    essOld_qu[i, 11, ] <-  quantile(essOld_dsl, 1-qs)
    essOld_qu[i, 12, ] <-  quantile(essOld_int, 1-qs)

    ## Sen's slope for dsl, intensity (using full 1920-2021 period)
    
    summ_senDry[i,1] <- RhatMat(trendsDrySim[[1]])
    summ_senDry[i,2] <- essMat(trendsDrySim[[1]])
    summ_senDry[i,3] <- esstMat(trendsDrySim[[1]])

    summ_senWet[i,1] <- RhatMat(trendsWetSim[[1]])
    summ_senWet[i,2] <- essMat(trendsWetSim[[1]])
    summ_senWet[i,3] <- esstMat(trendsWetSim[[1]])

    mlist <- as.mcmc.list(lapply(seq_along(starts), function(i)
                                 as.mcmc(trendsDrySim[[1]][chains[,i]])))
    summOld_senDry[i,1] <- gelman.diag(mlist)$psrf[1,1]
    summOld_senDry[i,2] <- effectiveSize(mlist)

    mlist <- as.mcmc.list(lapply(seq_along(starts), function(i)
                                 as.mcmc(trendsWetSim[[1]][chains[,i]])))
    summOld_senWet[i,1] <- gelman.diag(mlist)$psrf[1,1]
    summOld_senWet[i,2] <- effectiveSize(mlist)

    if(i == 3) { ## Winslow JJA
        summ_senDry_Winslow_JJA_1950_2021 <- RhatMat(trendsDrySim[[2]])
        summ_senDry_Winslow_JJA_1980_2021 <- RhatMat(trendsDrySim[[3]])
    }
}

save(summ_senDry, summ_senWet, rhat_qu, ess_qu,
     summOld_senDry, summOld_senWet, rhatOld_qu, essOld_qu,
     summ_senDry_Winslow_JJA_1950_2021, summ_senDry_Winslow_JJA_1980_2021,
     file = file.path(resultsDir, 'mixing.Rda'))

## Table 1: summ_senDry, summ_senWet
## Table 3 (Appendix): rhat_qu[,c(1:7,9,11,12),1], rhat_qu[,c(1:7,9,11,12),3]
