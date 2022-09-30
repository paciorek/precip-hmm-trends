#!/usr/bin/env Rscript
## Usage:
## calculate_spells.R AZ_Winslow 1 trend 5 0 20000
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
    data_fn <- 'CA_Quincy'
    season <- 5
    trend <- 'trend'
    config_model <- 5
    config_holdout <- 0
    niter <- 20000
}

## maxSimValRatio <- 2  # set sims greater than twice the maximum observed to that value

source(paste0('config-model-', config_model, '.R'))
source(paste0('config-holdout-', config_holdout, '.R'))
source('config-dirs.R')

source('load_data.R')
source('functions.R')

cmPerInch <- 2.54
dryCut <- 0.3  # cm below which a day is considered 'dry'

analyze_trend <- function(data, full = FALSE) {
    tmp <- modifiedmk::mkttest(data)
    if(full) return(tmp) else return(tmp[2])
}

## 1900:2017 for original AZ_Prescott, PA_station
## Zhang used 1925:2019, 1955:2019, 1976:2019
yrs <- 1900:2021 # 1920:2021
timeInts <- list(1920:2021, 1950:2021, 1980:2021)


output_pattern <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, sep = '-')
## Load in basic objects saved in output file using first output file.
load(file.path(output_dir, paste0(output_pattern, '-1.Rda')))

files <- list.files(output_dir, output_pattern)
 
nreps <- length(files)

burnin <- burnin_waic_frac*niter/thin
niter_per_chain <- niter/thin

starts <- seq(1, (nreps-1)*niter_per_chain+1, by = niter_per_chain) + burnin
ends <- starts + niter_per_chain - burnin - 1

keep <- starts[1]:ends[1]
for(i in 2:length(starts))
    keep <- c(keep, starts[i]:ends[i])


## Full posterior simulation-based analysis

## Get simulated data
fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
load(file.path(output_dir, paste0(fn, '.Rda')))

rSim <- rSim[ , , keep]

## rSim[rSim > maxSimValRatio * max(r)] <- maxSimValRatio * max(r)

nT <- dim(rSim)[1]

spellsSim <- list()
length(spellsSim) <- nT

## This uses about 4 GB
for(i in seq_len(nT)) 
    spellsSim[[i]] <- apply(rSim[i, , ], 2, get_spells, cutoff = dryCut/cmPerInch)

dryLenSim <- sapply(spellsSim, function(x) {
    sapply(x, function(y) mean(y$dryLen[ , 'dryLen']))
})

wetAvgSim <- sapply(spellsSim, function(x) {
    sapply(x, function(y) y$wetAvg)
})
wetAvgSim[is.na(wetAvgSim)] <- 0

if(FALSE) {
    wetMedSim <- sapply(spellsSim, function(x) {
        sapply(x, function(y) y$wetMed)
    })
    wetMedSim[is.na(wetMedSim)] <- 0
}

trendsDrySim <- list(); length(trendsDrySim) <- length(timeInts)
trendsWetSim <- list(); length(trendsWetSim) <- length(timeInts)

for(i in seq_along(timeInts)) {
    rg <- match(timeInts[[i]], yrs)
    trendsDrySim[[i]] <- apply(dryLenSim[ , rg], 1, analyze_trend)
    trendsWetSim[[i]] <- apply(wetAvgSim[ , rg], 1, analyze_trend)
}


## Imputation-based analysis

imputations <- combine_chains(data_fn, season, trend, config_model, config_holdout, niter, nburnin = 0)$imputations
imputations <- imputations[keep, ]

## imputations[imputations > maxSimValRatio * max(r)] <- maxSimValRatio * max(r)

## missingness was filled in by row
nIts <- nrow(imputations)
rImpute <- array(r, c(dim(r), nrow(imputations)))
for(i in 1:nIts) {
    tmp <- t(rImpute[ , , i])
    tmp[t(missing)] <- imputations[i, ]
    rImpute[ , , i] <- t(tmp)
}


nT <- dim(rImpute)[1]
spellsImp <- list(); length(spellsImp) <- nT

for(i in seq_len(nT)) 
    spellsImp[[i]] <- apply(rImpute[i, , ], 2, get_spells, cutoff = dryCut/cmPerInch)

dryLenImp <- sapply(spellsImp, function(x) {
    sapply(x, function(y) mean(y$dryLen[ , 'dryLen']))
})

wetAvgImp <- sapply(spellsImp, function(x) {
    sapply(x, function(y) y$wetAvg)
})
wetAvgImp[is.na(wetAvgImp)] <- 0

if(FALSE) {
    wetMedImp <- sapply(spellsImp, function(x) {
        sapply(x, function(y) y$wetMed)
    })
    wetMedImp[is.na(wetMedImp)] <- 0
}

trendsDryImp <- list(); length(trendsDryImp) <- length(timeInts)
pvalsDryImp <- list(); length(pvalsDryImp) <- length(timeInts)
statsDryImp <- list(); length(statsDryImp) <- length(timeInts)
trendsWetImp <- list(); length(trendsWetImp) <- length(timeInts)
pvalsWetImp <- list(); length(pvalsWetImp) <- length(timeInts)
statsWetImp <- list(); length(statsWetImp) <- length(timeInts)

analyze_trend_full <- function(data) {
   tmp <- modifiedmk::mkttest(data)
   return(tmp)
}

## Get average Sen's slope
for(i in seq_along(timeInts)) {
    rg <- match(timeInts[[i]], yrs)
    trendsDryImp[[i]] <- mean(apply(dryLenImp[ , rg], 1, analyze_trend))
    trendsWetImp[[i]] <- mean(apply(wetAvgImp[ , rg], 1, analyze_trend))
}

## Get p-value
for(i in seq_along(timeInts)) {
    rg <- match(timeInts[[i]], yrs)
    statsDryImp[[i]] <- apply(dryLenImp[ , rg], 1, analyze_trend, full = TRUE)
    S <- statsDryImp[[i]][3,]
    S[S>0] <- S[S>0]-1
    S[S<0] <- S[S<0]+1
    v <- mean(statsDryImp[[i]][4,]) + (1+1/nIts)*var(S)
    z <- mean(S)/sqrt(v)
    pvalsDryImp[[i]] <- 2*pnorm(-abs(z)) 

    statsWetImp[[i]] <- apply(wetAvgImp[ , rg], 1, analyze_trend, full = TRUE)
    S <- statsWetImp[[i]][3,]
    S[S>0] <- S[S>0]-1
    S[S<0] <- S[S<0]+1
    v <- mean(statsWetImp[[i]][4,]) + (1+1/nIts)*var(S)
    z <- mean(S)/sqrt(v)
    pvalsWetImp[[i]] <- 2*pnorm(-abs(z)) 
}

output_fn <- paste('spells', gsub(".csv", "", data_fn), season,
                   ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
output_fn <- file.path(output_dir, paste0(output_fn, '.Rda'))

save(spellsSim, dryLenSim, wetAvgSim, trendsDrySim, trendsWetSim,
     spellsImp, dryLenImp, wetAvgImp, trendsDryImp, trendsWetImp, pvalsDryImp, pvalsWetImp,
     statsDryImp, statsWetImp,
     file = output_fn)


