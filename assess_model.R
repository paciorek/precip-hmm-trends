set.seed(1)

sets <- 1:8

trend <- 'trend'
config_model <- 5
config_holdout <- 0
niter <- 20000

stns <- c(rep('AZ_Winslow', 4), rep('CA_Berkeley', 3), 'CA_Quincy')
seas <- c(1:4,1,2,4,5)

labs <- c("Winslow AZ, DJF", "Winslow AZ, MAM", "Winslow AZ, JJA", "Winslow AZ, SON",
          "Berkeley CA, DJF", "Berkeley CA, MAM", "Berkeley CA, SON", "Quincy CA, wet")

plotDir <- file.path('..', 'paper1', 'plots')
resultsDir <- file.path('..', 'paper1', 'results')

source(paste0('config-model-', config_model, '.R'))
source(paste0('config-holdout-', config_holdout, '.R'))
source('config-dirs.R')
source('functions.R')


cmPerInch <- 2.54
dryCut <- 0.3

## Not used in assessment, just for reporting outliers.
maxSimValRatio <- 2  # set sims greater than twice the maximum observed to that value

calc_rain <- function(x, cut = dryCut/cmPerInch, na.rm = FALSE) mean(x>cut, na.rm = na.rm)

calc_mean_rain  <- function(x, cut = dryCut/cmPerInch, na.rm = FALSE) mean(x[x>cut], na.rm = na.rm)

data_fn <- stns[1]
season <- seas[1]
output_pattern <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, sep = '-')

yrs <- 1900:2021 
timeInts <- list(1920:2021, 1950:2021, 1980:2021)
files <- list.files(output_dir, output_pattern)
 
nreps <- length(files)

burnin <- burnin_waic_frac*niter/thin
niter_per_chain <- niter/thin

starts <- seq(1, (nreps-1)*niter_per_chain+1, by = niter_per_chain) + burnin
ends <- starts + niter_per_chain - burnin - 1

keep <- starts[1]:ends[1]
for(i in 2:length(starts))
    keep <- c(keep, starts[i]:ends[i])

rg <- match(timeInts[[1]], yrs)

years <- yrs[rg]

rMissing <- list()
rSimFull <- rSimMissing <- list()
nDays <- list()
miss <- list()
spellsSimFull  <- dryLenSimFull <- list()

for(i in sets) {
    data_fn <- stns[i]
    season <- seas[i]
    source('load_data.R')

    nDays[[i]] <- dim(r)[[2]]
    
    r <- r[rg, ]
    missing <- missing[rg, ]
    miss[[i]] <- missing
    
    r[missing] <- NA
    rMissing[[i]] <- r

    fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
    load(file.path(output_dir, paste0(fn, '.Rda')))

    rSim <- rSim[rg, , keep]
    nIts <- dim(rSim)[3]

    ## hopefully not needed with new runs
    if(any(is.na(rSim)))
        print("NaNs found")
    if(any(rSim > maxSimValRatio*max(r))) {
        cat(sum(rSim > maxSimValRatio*max(r)), " extreme values found, with max of ",
            max(rSim, na.rm = TRUE), " compared to max in obs of ", max(r), "\n")
    }
    ## Not used after all.
    ## rSim[rSim > maxSimValRatio*max(r)] <- maxSimValRatio*max(r)

    rSimMissing[[i]] <- rSim
    rSimFull[[i]] <- rSim
    
    for(j in seq_len(nIts)) {
        tmp <- rSimMissing[[i]][,,j]
        tmp[missing] <- NA
        rSimMissing[[i]][,,j] <- tmp
    }

    fn <- paste('spells', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
    load(file.path(output_dir, paste0(fn, '.Rda')))

    spellsSimFull[[i]] <- spellsSim[rg]
    dryLenSimFull[[i]] <- dryLenSim[ , rg]  
}

## seasonal precip

pdf(file.path(plotDir, "diag_total.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {
    totalEmp <- apply(rMissing[[i]], 1, sum, na.rm = TRUE) * cmPerInch
    totalSim <- apply(rSimMissing[[i]], c(1,3), sum, na.rm = TRUE) * cmPerInch
    mn <- rowMeans(totalSim, na.rm = TRUE)
    qs <- apply(totalSim, 1, quantile, c(.05,.95), na.rm = TRUE)

    use <- rowSums(!miss[[i]])/nDays[[i]] > 0.75
    totalEmp[!use] <- NA
    mn[!use] <- NA
    qs[1, !use] <- NA
    qs[2, !use] <- NA

    plot(years, totalEmp, type = 'l', main = labs[i],
         xlab = 'year', ylab = 'seasonal precipitation (cm)')
    lines(years, mn, col = 'red')
    lines(years, qs[1, ], col = 'red', lty = 2)
    lines(years, qs[2, ], col = 'red', lty = 2)
}
dev.off()
    

## probability of rain by year

pdf(file.path(plotDir, "diag_rain_year.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {
    rainEmpYear <- apply(rMissing[[i]], 1, calc_rain, na.rm = TRUE)
    rainSimYear <- apply(rSimMissing[[i]], c(1,3), calc_rain, na.rm = TRUE)
    mn <- rowMeans(rainSimYear, na.rm = TRUE)
    qs <- apply(rainSimYear, 1, quantile, c(.05,.95), na.rm = TRUE)
    
    use <- rowSums(!miss[[i]])/nDays[[i]] > 0.75
    rainEmpYear[!use] <- NA
    mn[!use] <- NA
    qs[1, !use] <- NA
    qs[2, !use] <- NA

    plot(years, rainEmpYear, type = 'l', main = labs[i],
         xlab = 'year', ylab = 'probability of precipitation')
    lines(years, mn, col = 'red')
    lines(years, qs[1, ], col = 'red', lty = 2)
    lines(years, qs[2, ], col = 'red', lty = 2)
}
dev.off()

## probability of rain by season

pdf(file.path(plotDir, "diag_rain_seas.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {
    rainEmpSeas <- apply(rMissing[[i]], 2, calc_rain, na.rm = TRUE)
    rainSimSeas <- apply(rSimMissing[[i]], c(2,3), calc_rain, na.rm = TRUE)
    mn <- rowMeans(rainSimSeas, na.rm = TRUE)
    qs <- apply(rainSimSeas, 1, quantile, c(.05,.95), na.rm = TRUE)
    
    use <- colSums(!miss[[i]])/length(years) > 0.75
    rainEmpSeas[!use] <- NA
    mn[!use] <- NA
    qs[1, !use] <- NA
    qs[2, !use] <- NA

    plot(seq_len(nDays[[i]]), rainEmpSeas, type = 'l', main = labs[i],
         xlab = 'day of season', ylab = 'probability of precipitation',
         ylim = c(min(c(qs[1,],rainEmpSeas), na.rm = TRUE), max(c(qs[2,],rainEmpSeas), na.rm = TRUE)))
    lines(seq_len(nDays[[i]]), mn, col = 'red')
    lines(seq_len(nDays[[i]]), qs[1,], col = 'red', lty = 2)
    lines(seq_len(nDays[[i]]), qs[2,], col = 'red', lty = 2)
}
dev.off()

## mean daily rain by year

## can occasionally get sims with no precip greater than cutoff in entire season

pdf(file.path(plotDir, "diag_int_year.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {
    intEmpYear <- apply(rMissing[[i]], 1, calc_mean_rain, na.rm = TRUE) * cmPerInch
    intSimYear <- apply(rSimMissing[[i]], c(1,3), calc_mean_rain, na.rm = TRUE) * cmPerInch
    mn <- rowMeans(intSimYear, na.rm = TRUE)
    qs <- apply(intSimYear, 1, quantile, c(.05,.95), na.rm = TRUE)
    
    use <- rowSums(!miss[[i]])/nDays[[i]] > 0.75
    intEmpYear[!use] <- NA
    mn[!use] <- NA
    qs[1, !use] <- NA
    qs[2, !use] <- NA

    plot(years, intEmpYear, type = 'l', main = labs[i],
         xlab = 'year', ylab = 'precipitation intensity (cm)')
    lines(years, mn, col = 'red')
    lines(years, qs[1, ], col = 'red', lty = 2)
    lines(years, qs[2, ], col = 'red', lty = 2)
}
dev.off()


## mean daily rain by season

pdf(file.path(plotDir, "diag_int_seas.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {
    intEmpSeas <- apply(rMissing[[i]], 2, calc_mean_rain, na.rm = TRUE) * cmPerInch
    intSimSeas <- apply(rSimMissing[[i]], c(2,3), calc_mean_rain, na.rm = TRUE) * cmPerInch
    mn <- rowMeans(intSimSeas, na.rm = TRUE)
    qs <- apply(intSimSeas, 1, quantile, c(.05,.95), na.rm = TRUE)
    
    use <- colSums(!miss[[i]])/length(years) > 0.75
    intEmpSeas[!use] <- NA
    mn[!use] <- NA
    qs[1, !use] <- NA
    qs[2, !use] <- NA

    plot(seq_len(nDays[[i]]), intEmpSeas, type = 'l', main = labs[i],
         xlab = 'day of season', ylab = 'precipitation intensity (cm)',
         ylim = c(min(c(qs[1,],intEmpSeas), na.rm = TRUE), max(c(qs[2,],intEmpSeas), na.rm = TRUE)))
    lines(seq_len(nDays[[i]]), mn, col = 'red')
    lines(seq_len(nDays[[i]]), qs[1,], col = 'red', lty = 2)
    lines(seq_len(nDays[[i]]), qs[2,], col = 'red', lty = 2)
}
dev.off()


## state transitions

## by year

transEmpCatsYear <- transSimCatsYear <- mns <- qs <- list()

for(i in 1:8) {
    print(i)
    transEmpCatsYear[[i]] <- apply(rMissing[[i]], 1, calc_trans, c(dryCut/cmPerInch, 1/cmPerInch))
    transSimCatsYear[[i]] <- apply(rSimMissing[[i]], c(1,3), calc_trans, c(dryCut/cmPerInch, 1/cmPerInch))

    mns[[i]] <- apply(transSimCatsYear[[i]], c(1,2), mean, na.rm = TRUE)
    qs[[i]] <- apply(transSimCatsYear[[i]], c(1,2), quantile, c(.05, .95), na.rm = TRUE)

    use <- rowSums(!miss[[i]])/nDays[[i]] > 0.75
    transEmpCatsYear[[i]][, !use] <- NA
    mns[[i]][ , !use] <- NA
    qs[[i]][1, , !use] <- NA
    qs[[i]][2, , !use] <- NA
}

save(transEmpCatsYear, transSimCatsYear, file = file.path(resultsDir, 'trans_by_year.Rda'))

## from dry

cols <- c(1,4,7)
pdf(file.path(plotDir, "diag_trans_from_dry_year.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(8, 3))
for(i in sets) {
    if(i==1) {
        main  <- rep(c('Transition to dry (0-0.3 cm)',
                   'Transition to moist (0.3-1.0 cm)',
                   'Transition to wet (> 1.0 cm)'), each = 3)
    } else main  <- ''

    for(col in cols) {
        plot(years, transEmpCatsYear[[i]][col, ], type = "l",
             xlab = 'year', ylab = 'transition prob.', main = main[col], cex.main = 1,
             ylim = c(min(c(qs[[i]][1,col,],transEmpCatsYear[[i]][col, ]), na.rm = TRUE),
                      max(c(qs[[i]][2,col,],transEmpCatsYear[[i]][col, ]), na.rm = TRUE)))
        lines(years, mns[[i]][col, ], col = 'red', lty = 1)
        lines(years, qs[[i]][1,col,], col = 'red', lty = 2)
        lines(years, qs[[i]][2,col,], col = 'red', lty = 2)
        if(col == cols[1])
            mtext(labs[i], side = 2, line = 3, cex = 0.8)
    }
}
dev.off()

## from moist
cols <- c(2,5,8)

pdf(file.path(plotDir, "diag_trans_from_moist_year.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(8, 3))
for(i in sets) {
    if(i==1) {
        main  <- rep(c('Transition to dry (0-0.3 cm)',
                   'Transition to moist (0.3-1.0 cm)',
                   'Transition to wet (> 1.0 cm)'), each = 3)
    } else main  <- ''

    for(col in cols) {
        plot(years, transEmpCatsYear[[i]][col, ], type = "l",
             xlab = 'year', ylab = 'transition prob.', main = main[col], cex.main = 1,
             ylim = c(min(c(qs[[i]][1,col,],transEmpCatsYear[[i]][col, ]), na.rm = TRUE),
                      max(c(qs[[i]][2,col,],transEmpCatsYear[[i]][col, ]), na.rm = TRUE)))
        lines(years, mns[[i]][col, ], col = 'red', lty = 1)
        lines(years, qs[[i]][1,col,], col = 'red', lty = 2)
        lines(years, qs[[i]][2,col,], col = 'red', lty = 2)
        if(col == cols[1])
            mtext(labs[i], side = 2, line = 3, cex = 0.8)
    }
}
dev.off()

## from wet

cols <- c(3,6,9)

pdf(file.path(plotDir, "diag_trans_from_wet_year.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(8, 3))
for(i in sets) {
    if(i==1) {
        main  <- rep(c('Transition to dry (0-0.3 cm)',
                   'Transition to moist (0.3-1.0 cm)',
                   'Transition to wet (> 1.0 cm)'), each = 3)
    } else main  <- ''

    for(col in cols) {
        plot(years, transEmpCatsYear[[i]][col, ], type = "l",
             xlab = 'year', ylab = 'transition prob.', main = main[col], cex.main = 1,
             ylim = c(min(c(qs[[i]][1,col,],transEmpCatsYear[[i]][col, ]), na.rm = TRUE),
                      max(c(qs[[i]][2,col,],transEmpCatsYear[[i]][col, ]), na.rm = TRUE)))
        lines(years, mns[[i]][col, ], col = 'red', lty = 1)
        lines(years, qs[[i]][1,col,], col = 'red', lty = 2)
        lines(years, qs[[i]][2,col,], col = 'red', lty = 2)
        if(col == cols[1])
            mtext(labs[i], side = 2, line = 3, cex = 0.8)
    }
}
dev.off()


## now by day of season

transEmpCatsSeas <- transSimCatsSeas <- mns <- qs <- list()

for(i in 1:8) {
    print(i)
    transEmpCatsSeas[[i]] <- apply(rMissing[[i]], 2, calc_trans, c(dryCut/cmPerInch, 1/cmPerInch))
    transSimCatsSeas[[i]] <- apply(rSimMissing[[i]], c(2,3), calc_trans, c(dryCut/cmPerInch, 1/cmPerInch))

    mns[[i]] <- apply(transSimCatsSeas[[i]], c(1,2), mean, na.rm = TRUE)
    qs[[i]] <- apply(transSimCatsSeas[[i]], c(1,2), quantile, c(.05, .95), na.rm = TRUE)

    use <- colSums(!miss[[i]])/length(years) > 0.75
    transEmpCatsSeas[[i]][, !use] <- NA
    mns[[i]][ , !use] <- NA
    qs[[i]][1, , !use] <- NA
    qs[[i]][2, , !use] <- NA
}

save(transEmpCatsSeas, transSimCatsSeas, file = file.path(resultsDir, 'trans_by_seas.Rda'))


## from dry

cols <- c(1,4,7)

pdf(file.path(plotDir, "diag_trans_from_dry_seas.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(8, 3))
for(i in sets) {
    if(i==1) {
        main  <- rep(c('Transition to dry (0-0.3 cm)',
                   'Transition to moist (0.3-1.0 cm)',
                   'Transition to wet (> 1.0 cm)'), each = 3)
    } else main  <- ''

    for(col in cols) {
        plot(seq_len(nDays[[i]]), transEmpCatsSeas[[i]][col, ], type = "l",
             xlab = 'day of season', ylab = 'transition prob.', main = main[col], cex.main = 1,
             ylim = c(min(c(qs[[i]][1,col,],transEmpCatsSeas[[i]][col, ]), na.rm = TRUE),
                      max(c(qs[[i]][2,col,],transEmpCatsSeas[[i]][col, ]), na.rm = TRUE)))
        lines(seq_len(nDays[[i]]), mns[[i]][col, ], col = 'red', lty = 1)
        lines(seq_len(nDays[[i]]), qs[[i]][1,col,], col = 'red', lty = 2)
        lines(seq_len(nDays[[i]]), qs[[i]][2,col,], col = 'red', lty = 2)
        if(col == cols[1])
            mtext(labs[i], side = 2, line = 3, cex = 0.8)
    }
}
dev.off()

## from moist

cols <- c(2,5,8)

pdf(file.path(plotDir, "diag_trans_from_moist_seas.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(8, 3))
for(i in sets) {
    if(i==1) {
        main  <- rep(c('Transition to dry (0-0.3 cm)',
                   'Transition to moist (0.3-1.0 cm)',
                   'Transition to wet (> 1.0 cm)'), each = 3)
    } else main  <- ''

    for(col in cols) {
        plot(seq_len(nDays[[i]]), transEmpCatsSeas[[i]][col, ], type = "l",
             xlab = 'day of season', ylab = 'transition prob.', main = main[col], cex.main = 1,
             ylim = c(min(c(qs[[i]][1,col,],transEmpCatsSeas[[i]][col, ]), na.rm = TRUE),
                      max(c(qs[[i]][2,col,],transEmpCatsSeas[[i]][col, ]), na.rm = TRUE)))
        lines(seq_len(nDays[[i]]), mns[[i]][col, ], col = 'red', lty = 1)
        lines(seq_len(nDays[[i]]), qs[[i]][1,col,], col = 'red', lty = 2)
        lines(seq_len(nDays[[i]]), qs[[i]][2,col,], col = 'red', lty = 2)
        if(col == cols[1])
            mtext(labs[i], side = 2, line = 3, cex = 0.8)
    }
}
dev.off()

## from wet

cols <- c(3,6,9)

pdf(file.path(plotDir, "diag_trans_from_wet_seas.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(8, 3))
for(i in sets) {
    if(i==1) {
        main  <- rep(c('Transition to dry (0-0.3 cm)',
                   'Transition to moist (0.3-1.0 cm)',
                   'Transition to wet (> 1.0 cm)'), each = 3)
    } else main  <- ''

    for(col in cols) {
        plot(seq_len(nDays[[i]]), transEmpCatsSeas[[i]][col, ], type = "l",
             xlab = 'day of season', ylab = 'transition prob.', main = main[col], cex.main = 1,
             ylim = c(min(c(qs[[i]][1,col,],transEmpCatsSeas[[i]][col, ]), na.rm = TRUE),
                      max(c(qs[[i]][2,col,],transEmpCatsSeas[[i]][col, ]), na.rm = TRUE)))
        lines(seq_len(nDays[[i]]), mns[[i]][col, ], col = 'red', lty = 1)
        lines(seq_len(nDays[[i]]), qs[[i]][1,col,], col = 'red', lty = 2)
        lines(seq_len(nDays[[i]]), qs[[i]][2,col,], col = 'red', lty = 2)
        if(col == cols[1])
            mtext(labs[i], side = 2, line = 3, cex = 0.8)
    }
}
dev.off()

## Q-Q plots

### set up averaging matrices for 3-day, 10-day precip
qq1 <- qq3 <- qq10 <- r3 <- r10 <- list()

for(i in sets) {

    d <- ncol(rMissing[[i]])
    mat3 <- matrix(0, d, d-2)
    for(j in 1:(d-2)) 
        mat3[j:(j+2) ,j] <- 1
    
    mat10 <- matrix(0, d, d - 9)
    for(j in 1:(d-9)) 
        mat10[j:(j+9) ,j] <- 1

    ## single day
    
    qq1[[i]] <- sapply(seq_len(nIts), function(j) {
        out <- qqplot(rMissing[[i]], rSimMissing[[i]][,,j], plot.it = FALSE)
        return(out$y)
    })
    
    ## three day
    
    r3[[i]] <- rMissing[[i]]%*%mat3
    qq3[[i]] <- sapply(seq_len(nIts), function(j) {
        tmp <- rSimMissing[[i]][,,j]%*%mat3
        out <- qqplot(r3[[i]], tmp, plot.it = FALSE)
        return(out$y)
    })
    
    ## ten day
    
    r10[[i]] <- rMissing[[i]]%*%mat10
    qq10[[i]] <- sapply(seq_len(nIts), function(j) {
        tmp <- rSimMissing[[i]][,,j]%*%mat10
        out <- qqplot(r10[[i]], tmp, plot.it = FALSE)
        return(out$y)
    })
}

    
## if real data not consistent with simulation dist, we'd expect the sim-based interval
## to lie away from the 1:1 line, which is what we expect for data with the same dist'n as the real data

pdf(file.path(plotDir, "diag_qq.pdf"), height = 11, width = 8.5)
## par(mai = c(.4,.6,.2,.1), mgp = c(1.8, 0.7, 0))
par(mai = c(.4,.2,.1,.1), mgp = c(1.8, 0.7, 0), omi = c(0,.35,.3,0))
par(mfrow = c(8, 3))

for(i in sets) {
    d <- ncol(rMissing[[i]])
    mat3 <- matrix(0, d, d-2)
    for(j in 1:(d-2)) 
        mat3[j:(j+2) ,j] <- 1
    
    mat10 <- matrix(0, d, d - 9)
    for(j in 1:(d-9)) 
        mat10[j:(j+9) ,j] <- 1

    if(i==1) {
        main  <- c('Daily precipitation',
                   'Three-day precipitation',
                   'Ten-day precipitation')
    } else main  <- rep('', 3)

    base1 <- qqplot(rMissing[[i]], rSimMissing[[i]][,,1], plot.it = FALSE)
    qus1 <- apply(qq1[[i]], 1, quantile, c(.025,.5, .975), na.rm = TRUE)

    plot(base1$x*cmPerInch, qus1[2,]*cmPerInch, pch=16, cex = .5,
         xlab = 'daily precipitation (cm)', ylab = '', 
         ylim = range(base1$x*cmPerInch))
    mtext(main[1], side = 3, line = 1)
    lines(base1$x*cmPerInch, qus1[1,]*cmPerInch, lty = 2)
    lines(base1$x*cmPerInch, qus1[3,]*cmPerInch, lty = 2)
    abline(0,1)
    mtext(labs[i], side = 2, line = 2, cex = 0.8)
    
    tmp <- rSimMissing[[i]][,,1]%*%mat3
    base3 <- qqplot(r3[[i]], tmp, plot.it = FALSE)
    qus3 <- apply(qq3[[i]], 1, quantile, c(.025,.5, .975), na.rm = TRUE)
    
    plot(base3$x*cmPerInch, qus3[2,]*cmPerInch, pch=16,
         xlab = 'three-day precipitation (cm)', ylab = '', 
         cex = .5, ylim = range(base3$x*cmPerInch))
    mtext(main[2], side = 3, line = 1)
    lines(base3$x*cmPerInch, qus3[1,]*cmPerInch, lty = 2)
    lines(base3$x*cmPerInch, qus3[3,]*cmPerInch, lty = 2)
    abline(0,1)
    
    tmp <- rSimMissing[[i]][,,1]%*%mat10
    base10 <- qqplot(r10[[i]], tmp, plot.it = FALSE)
    qus10 <- apply(qq10[[i]], 1, quantile, c(.025,.5, .975), na.rm = TRUE)
    
    plot(base10$x*cmPerInch, qus10[2,]*cmPerInch, pch=16,
         xlab = 'ten-day precipitation (cm)', ylab = '', 
         cex = .5, ylim = range(base10$x*cmPerInch))
    mtext(main[3], side = 3, line = 1)
    lines(base10$x*cmPerInch, qus10[1,]*cmPerInch, lty = 2)
    lines(base10$x*cmPerInch, qus10[3,]*cmPerInch, lty = 2)
    abline(0,1)

}
dev.off()

## Dry spell length


pdf(file.path(plotDir, "diag_dsl.pdf"), height = 11, width = 8.5)
#postscript(file.path(plotDir, "diag_dsl.ps"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {

    tmp <- rMissing[[i]]
    tmp[miss[[i]]] <- 0 # otherwise get_spells fails
    spellsEmp <- apply(tmp, 1, get_spells, cutoff = dryCut/cmPerInch)
    
    dryLenEmp <- sapply(spellsEmp, function(x) {
        mean(x$dryLen[ , 'dryLen'])
    })

    mn <- colMeans(dryLenSimFull[[i]], na.rm = TRUE)
    qs <- apply(dryLenSimFull[[i]], 2, quantile, c(.05,.95), na.rm = TRUE)

    use <- rowSums(!miss[[i]])/nDays[[i]] == 1
    dryLenEmp[!use] <- NA
    ## We do plot simulated values where there are missing empirical values to illustrate trends.
    ## mn[!use] <- NA
    ## qs[1, !use] <- NA
    ## qs[2, !use] <- NA

    plot(years, dryLenEmp, type = 'l', main = labs[i],
         xlab = 'year', ylab = 'mean dry spell length (days)')
    lines(years, mn, col = 'red')
    lines(years, qs[1, ], col = 'red', lty = 2)
    lines(years, qs[2, ], col = 'red', lty = 2)
}
dev.off()
   
## Dry spell length - QQ

pdf(file.path(plotDir, "diag_dsl_qq.pdf"), height = 11, width = 8.5)
par(mai = c(.4,.4,.3,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(4, 2))
for(i in sets) {
    
    spellsSimTmp <- spellsSimFull[[i]]
    spellsSimTmp <- spellsSimTmp[rowSums(miss[[i]]) == 0]
    lens <- sapply(seq_len(nIts), function(i)
        length(unlist(sapply(seq_along(spellsSimTmp), function(idx) spellsSimTmp[[idx]][[i]]$dryLen[,'dryLen']))))

    tmp <- rMissing[[i]]
    tmp[miss[[i]]] <- 0 # otherwise get_spells fails
    spellsEmp <- apply(tmp, 1, get_spells, cutoff = dryCut/cmPerInch)
    spellsEmp <- spellsEmp[rowSums(miss[[i]]) == 0]
    spellsEmp <- unlist(sapply(seq_along(spellsEmp), function(i) spellsEmp[[i]]$dryLen[,'dryLen']))
    
    ## need no more empirical dry spells than min number over simulated so that have constant x-axis
    spellsEmpSub <- sample(spellsEmp, min(lens), replace = FALSE)

    qqLen <- sapply(seq_len(nIts), function(i) {
        vals <- unlist(sapply(seq_along(spellsSimTmp), function(idx)
            spellsSimTmp[[idx]][[i]]$dryLen[,'dryLen']))
        return(qqplot(spellsEmpSub, vals, plot.it = FALSE)$y)
    })

    base <- qqplot(spellsEmpSub, unlist(sapply(seq_along(spellsSimTmp), function(idx) spellsSimTmp[[idx]][[1]]$dryLen[,'dryLen'])),
               plot.it = FALSE)

    qs <- apply(qqLen, 1, quantile, c(.025,.5,.975))
    
    plot(base$x, qs[2,], pch=16, cex = .5, main = labs[i], xlab = 'observed dry spell length (days)',
         ylab = 'simulated dry spell length (days)',
         ylim = range(base$x))
    lines(base$x, qs[1,], lty = 2)
    lines(base$x, qs[3,], lty = 2)
    abline(0,1)
}
dev.off()

    
  
