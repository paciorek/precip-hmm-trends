
## Extract holdout LL and WAIC values

resultsDir <- file.path('..', 'paper1', 'results')

source('config-dirs.R')

locn <- c(rep('AZ_Winslow', 3), rep('AZ_Winslow', 3), rep('CA_Berkeley', 3), rep('CA_Quincy',3)) 
seas <- rep(c(1,3,1,5), each = 3)
tr <- rep(c("notrend","trend","trend"), 4)
config_mod <- rep(c(5,5,4), 4)

output <- data.frame(locn = locn, season = seas, trend = tr, dist = config_mod,
                     ll_train = NA, ll_train_sd = NA,
                     ll_test1 = NA, ll_test1_sd = NA, ll_test2 = NA, ll_test2_sd = NA, waic = NA, pwaic = NA)

output_pattern <- paste('waic', locn, seas, ifelse(tr == 'trend', 'trend', 'notrend'), config_mod, 1, sep = '-')

for(i in seq_along(locn)) {
    
    load(file.path(output_dir, paste0(output_pattern[i], '.Rda')))
    
    output[i, 'waic'] <- waic$WAIC
    output[i, 'pwaic'] <- waic$pWAIC
    output[i, 'll_train'] <- mean(ll$ll_train, na.rm = TRUE)
    output[i, 'll_train_sd'] <- sd(ll$ll_train, na.rm = TRUE)
    output[i, 'll_test1'] <- mean(ll$ll_test1, na.rm = TRUE)
    output[i, 'll_test1_sd'] <- sd(ll$ll_test1, na.rm = TRUE)
    output[i, 'll_test2'] <- mean(ll$ll_test2, na.rm = TRUE)
    output[i, 'll_test2_sd'] <- sd(ll$ll_test2, na.rm = TRUE)
}


## Table 2 values
save(output, file = file.path(resultsDir, 'll_waic.Rda'))




## check Ga vs GPD tails

data_fn <- 'CA_Berkeley'
season <- 1
trend <- 'trend'
config_model <- 4
config_holdout <- 1
niter <- 20000

source(paste0('config-model-', config_model, '.R'))
source(paste0('config-holdout-', config_holdout, '.R'))
source('config-dirs.R')
source('functions.R')

source('load_data.R')

output_pattern <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, sep = '-')

files <- list.files(output_dir, output_pattern)
 
nreps <- length(files)

burnin <- burnin_waic_frac*niter/thin
niter_per_chain <- niter/thin

starts <- seq(1, (nreps-1)*niter_per_chain+1, by = niter_per_chain) + burnin
ends <- starts + niter_per_chain - burnin - 1

keep <- starts[1]:ends[1]
for(i in 2:length(starts))
    keep <- c(keep, starts[i]:ends[i])

yrs <- 1900:2021 # 1920:2021
timeInts <- list(1920:2021, 1950:2021, 1980:2021)

## Load in basic objects saved in output file using first output file.
load(file.path(output_dir, paste0(output_pattern, '-1.Rda')))

fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
load(file.path(output_dir, paste0(fn, '.Rda')))

rg <- match(timeInts[[1]], yrs)
nIts <- dim(rSim)[3]

rSimMissing <- rSim
for(i in seq_len(nIts)) {
    tmp <- rSimMissing[,,i]
    tmp[missing] <- NA
    rSimMissing[,,i] <- tmp
}

rSimMissing <- rSimMissing[rg, , keep]

rSave <- r
r[missing] <- NA

r <- r[rg, ]

nItsKeep <- dim(rSimMissing)[3]

qq1 <- sapply(seq_len(nItsKeep), function(i) {
    out <- qqplot(r, rSimMissing[,,i], plot.it = FALSE)
    return(out$y)
})

base1 <- qqplot(r, rSimMissing[,,1], plot.it = FALSE)
qus1 <- apply(qq1, 1, quantile, c(.025,.5, .975))

maxval <- max(r, na.rm = TRUE)

plot(base1$x, qus1[2,], pch=16, cex = .5, xlim = c(0,maxval), ylim = c(0,maxval))
lines(base1$x, qus1[1,], lty = 2)
lines(base1$x, qus1[3,], lty = 2)
abline(0,1)

