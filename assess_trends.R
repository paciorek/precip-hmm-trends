library(ggplot2)

sets <- 1:8

trend <- 'trend'
config_model <- 5
config_holdout <- 0
niter <- 20000

stns <- c(rep('AZ_Winslow', 4), rep('CA_Berkeley', 3), 'CA_Quincy')
seas <- c(1:4, 1,2,4,5)

labs <- c("Winslow AZ, DJF", "Winslow AZ, MAM", "Winslow AZ, JJA", "Winslow AZ, SON",
          "Berkeley CA, DJF", "Berkeley CA, MAM", "Berkeley CA, SON", "Quincy CA, wet")
seasLabs <- c("DJF","MAM","JJA","SON","DJF","MAM","SON","wet")

plotDir <- file.path('..', 'paper1', 'plots')

source(paste0('config-model-', config_model, '.R'))
source(paste0('config-holdout-', config_holdout, '.R'))
source('config-dirs.R')
source('functions.R')

cmPerInch <- 2.54
dryCut <- 0.3
## maxSimValRatio <- 2  # set sims greater than twice the maximum observed to that value


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
nDays <- list()
miss <- list()
spellsImpFull  <- dryLenImpFull <- wetAvgImpFull <- list()
spellsSimFull <- dryLenSimFull <- wetAvgSimFull <- list()

trendsDryImpFull <- pvalsDryImpFull <- trendsDrySimFull <- trendsWetImpFull <- pvalsWetImpFull <- trendsWetSimFull <- list()

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

    fn <- paste('spells', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
    load(file.path(output_dir, paste0(fn, '.Rda')))

    dryLenSimFull[[i]] <- dryLenSim
    wetAvgSimFull[[i]] <- wetAvgSim
    spellsSimFull[[i]] <- spellsSim
    
    dryLenImpFull[[i]] <- dryLenImp
    wetAvgImpFull[[i]] <- wetAvgImp 
    spellsImpFull[[i]] <- spellsImp

    trendsDryImpFull[[i]] <- trendsDryImp
    pvalsDryImpFull[[i]] <- pvalsDryImp
    trendsDrySimFull[[i]] <- trendsDrySim

    trendsWetImpFull[[i]] <- trendsWetImp
    pvalsWetImpFull[[i]] <- pvalsWetImp
    trendsWetSimFull[[i]] <- trendsWetSim
}


pdf(file.path(plotDir, "trends.pdf"), height = 11, width = 8.5)
    
par(mai = c(.4,.5,.05,.1), mgp = c(1.8, 0.7, 0), omi = c(0,.3,.3,0))
par(mfrow = c(8, 2))

main  <- c('Mean dry spell length', 'Precipitation intensity')

for(i in sets) {
    mn <- colMeans(dryLenImpFull[[i]][ , rg])
    qs <- apply(dryLenImpFull[[i]][ , rg], 2, quantile, c(.05, .95))
    plot(years, mn, pch = 16, cex = 0.75, xlab = 'year', ylab = 'dry spell len. (days)')
    lines(years, qs[1, ], lty = 2, lwd = 0.7)
    lines(years, qs[2, ], lty = 2)
    lines(years, mn, lwd = 1.3)
    if(i == 1)
        mtext(main[1], side = 3, line = 1)

    mtext(labs[i], side = 2, line = 3, cex = 0.8)

    mn <- colMeans(wetAvgImpFull[[i]][ , rg]) * cmPerInch
    qs <- apply(wetAvgImpFull[[i]][ , rg], 2, quantile, c(.05, .95)) * cmPerInch
    plot(years, mn, pch = 16, cex = 0.75, xlab = 'year', ylab = 'precip. intensity (cm)')
    lines(years, qs[1, ], lty = 2, lwd = 0.7)
    lines(years, qs[2, ], lty = 2)
    lines(years, mn, lwd = 1.3)
    if(i == 1)
        mtext(main[2], side = 3, line = 1)  
}
dev.off()

## Winslow, Berkeley inference

yrBlocks <- sapply(timeInts, function(x) paste(min(x), max(x), sep = '-'))
sz <- 2.5

szPt <- .5

pWet <- pDry <- list()

ylimsDry <- list(c(-0.7,0.7), c(-1,1), c(-0.5,0.5), c(-0.5,0.5), c(-.15,.15), c(-.3,.3), c(-.5,.5))
ylimsWet <- list(c(-.006,.006), c(-.008,.008), c(-.008,.008), c(-.01, .01), c(-.008,.008), c(-.01,.01), c(-.015,.015))

nIts <- length(trendsDrySim[[1]])

for(j in 1:7) {
    dfImp <- data.frame(x = yrBlocks, y = 10*sapply(trendsDryImpFull[[j]], mean), significance = pvalsDryImpFull[[j]] < 0.05)
    dfSim <- data.frame(x = rep(yrBlocks, each = nIts), y = 10*unlist(trendsDrySimFull[[j]]))
    dfText <- data.frame(x = yrBlocks, y = 10 * rep(ylimsDry[[j]][1], 3), pval = paste0("p = ", sapply(pvalsDryImpFull[[j]], signif, 2)))
    pDry[[j]] <- ggplot(dfSim, aes(factor(x), y)) + geom_violin() + ylab("Sen's slope (days/decade)") + xlab('') + ylim(10 * ylimsDry[[j]]) +
        geom_point(data = dfImp, aes(factor(x), y, shape = significance, size = szPt, color = "rd")) + scale_shape_manual(values = c(21,16)) +
        theme_minimal() + guides(shape = "none", size = "none", color = "none") +
        geom_text(data = dfText, aes(factor(x), y, label = pval), size = sz, color = "red") + labs(subtitle = seasLabs[j])
    if(j == 1 || j == 5)
        pDry[[j]] <- pDry[[j]] + ggtitle("Mean dry spell length") 
    
    dfImp <- data.frame(x = yrBlocks, y = cmPerInch * 10 * sapply(trendsWetImpFull[[j]], mean), significance = pvalsWetImpFull[[j]] < 0.05)
    dfSim <- data.frame(x = rep(yrBlocks, each = nIts), y = cmPerInch * 10 * unlist(trendsWetSimFull[[j]]))
    dfText <- data.frame(x = yrBlocks, y = cmPerInch * 10 * rep(ylimsWet[[j]][1], 3), pval = paste0("p = ", sapply(pvalsWetImpFull[[j]], signif, 2)))
    pWet[[j]] <- ggplot(dfSim, aes(factor(x), y)) + geom_violin() + ylab("Sen's slope (cm/decade)") + xlab('') + ylim(cmPerInch * 10 * ylimsWet[[j]]) + 
        geom_point(data = dfImp, aes(factor(x), y, shape = significance, size = szPt, color = 'red')) + scale_shape_manual(values = c(21,16)) +
        theme_minimal() + guides(shape = "none", size = "none", color = "none") + 
        geom_text(data = dfText, aes(factor(x), y, label = pval), size = sz, color = 'red') + labs(subtitle = seasLabs[j])
    if(j == 1 || j == 5)
        pWet[[j]] <- pWet[[j]] + ggtitle("Precipitation intensity")
    
}

### warnings are from values being cutoff by the ylim I am setting.

    
## Winslow

pdf(file.path(plotDir, 'inference_winslow.pdf'), width = 8.5, height = 9)
gridExtra::grid.arrange(pDry[[1]], pWet[[1]],
                        pDry[[2]], pWet[[2]],
                        pDry[[3]], pWet[[3]],
                        pDry[[4]], pWet[[4]], nrow = 4, ncol = 2)
dev.off()

## Berkeley

pdf(file.path(plotDir, 'inference_berkeley.pdf'), width = 8.5, height = 7)
gridExtra::grid.arrange(pDry[[5]], pWet[[5]],
                        pDry[[6]], pWet[[6]],
                        pDry[[7]], pWet[[7]], nrow = 3, ncol = 2)
dev.off()

## Quincy metrics

i <- 8
data_fn <- stns[i]
season <- seas[i]

source('load_data.R')

## get imputations
imputations <- combine_chains(data_fn, season, trend, config_model, config_holdout, niter, output_dir = output_dir, nburnin = 0)$imputations
imputations <- imputations[keep, ]

## imputations[imputations > maxSimValRatio * max(r)] <- maxSimValRatio * max(r)


nIts <- nrow(imputations)
rImpute <- array(r, c(dim(r), nrow(imputations)))
for(j in 1:nIts) {
    tmp <- t(rImpute[ , , j])
    tmp[t(missing)] <- imputations[j, ]
    rImpute[ , , j] <- t(tmp)
}

nT <- dim(rImpute)[1]
nS <- dim(rImpute)[2]

wetLongImp <- apply(rImpute, c(1,3), function(data) {
    max(sapply(1:(nS-39), function(i) sum(data[i:(i+39)]))) })
wetLongImp <- t(wetLongImp) * cmPerInch

wetAmtImp <- sapply(spellsImpFull[[i]], function(x) {
    sapply(x, function(y) mean(y$wetAmt[ , 'wetAmt']))
})
wetAmtImp <- wetAmtImp * cmPerInch

wetNumImp <- sapply(spellsImpFull[[i]], function(x) {
    sapply(x, function(y) nrow(y$wetLen))
})

rg <- match(timeInts[[1]], yrs)

pdf(file.path(plotDir, "trends_quincy.pdf"), height = 2.5, width = 10)
    
par(mai = c(.4,.5,.05,.1), mgp = c(1.8, 0.7, 0))
par(mfrow = c(1,3))

mn <- colMeans(wetNumImp[,rg])
qs <- apply(wetNumImp[,rg], 2, quantile, c(.05, .95))
plot(years, mn, pch = 16, cex = 0.75, xlab = 'year', ylab = 'number of events', ylim = c(7,40))
lines(years, qs[1, ], lty = 2, lwd = 0.7)
lines(years, qs[2, ], lty = 2)
lines(years, mn, lwd = 1.3)

mn <- colMeans(wetAmtImp[,rg])
qs <- apply(wetAmtImp[,rg], 2, quantile, c(.05, .95))
plot(years, mn, pch = 16, cex = 0.75, xlab = 'year', ylab = 'mean event precipitation (cm)', ylim = c(1,10))
lines(years, qs[1, ], lty = 2, lwd = 0.7)
lines(years, qs[2, ], lty = 2)
lines(years, mn, lwd = 1.3)

mn <- colMeans(wetLongImp[,rg])
qs <- apply(wetLongImp[,rg], 2, quantile, c(.05, .95))
plot(years, mn, pch = 16, cex = 0.75, xlab = 'year', ylab = 'maximum 40-day precipitation (cm)', ylim = c(5,100))
lines(years, qs[1, ], lty = 2, lwd = 0.7)
lines(years, qs[2, ], lty = 2)
lines(years, mn, lwd = 1.3)

dev.off()


## Oroville analysis

### Inference based on imputations

trendsLongImp <- list(); length(trendsLongImp) <- length(timeInts)
pvalsLongImp <- list(); length(pvalsLongImp) <- length(timeInts)
statsLongImp <- list(); length(statsLongImp) <- length(timeInts)
trendsAmtImp <- list(); length(trendsAmtImp) <- length(timeInts)
pvalsAmtImp <- list(); length(pvalsAmtImp) <- length(timeInts)
statsAmtImp <- list(); length(statsAmtImp) <- length(timeInts)
trendsNumImp <- list(); length(trendsNumImp) <- length(timeInts)
pvalsNumImp <- list(); length(pvalsNumImp) <- length(timeInts)
statsNumImp <- list(); length(statsNumImp) <- length(timeInts)

for(j in seq_along(timeInts)) {
    rg <- match(timeInts[[j]], yrs)

    statsLongImp[[j]] <- apply(wetLongImp[ , rg], 1, analyze_trend, full = TRUE)
    trendsLongImp[[j]] <- statsLongImp[[j]][2,]
    S <- statsLongImp[[j]][3,]
    S[S>0] <- S[S>0]-1
    S[S<0] <- S[S<0]+1
    v <- mean(statsLongImp[[j]][4,]) + (1+1/nIts)*var(S)
    z <- mean(S)/sqrt(v)
    pvalsLongImp[[j]] <- 2*pnorm(-abs(z)) 

    statsAmtImp[[j]] <- apply(wetAmtImp[ , rg], 1, analyze_trend, full = TRUE)
    trendsAmtImp[[j]] <- statsAmtImp[[j]][2,]
    S <- statsAmtImp[[j]][3,]
    S[S>0] <- S[S>0]-1
    S[S<0] <- S[S<0]+1
    v <- mean(statsAmtImp[[j]][4,]) + (1+1/nIts)*var(S)
    z <- mean(S)/sqrt(v)
    pvalsAmtImp[[j]] <- 2*pnorm(-abs(z)) 

    statsNumImp[[j]] <- apply(wetNumImp[ , rg], 1, analyze_trend, full = TRUE)
    trendsNumImp[[j]] <- statsNumImp[[j]][2,]
    S <- statsNumImp[[j]][3,]
    S[S>0] <- S[S>0]-1
    S[S<0] <- S[S<0]+1
    v <- mean(statsNumImp[[j]][4,]) + (1+1/nIts)*var(S)
    z <- mean(S)/sqrt(v)
    pvalsNumImp[[j]] <- 2*pnorm(-abs(z)) 
}

## Inference based on simulated

fn <- paste('sims', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
load(file.path(output_dir, paste0(fn, '.Rda')))

rSim <- rSim[ , , keep]
## rSim[rSim > maxSimValRatio * max(r)] <- maxSimValRatio * max(r)

nT <- dim(rSim)[1]

wetLongSim <- apply(rSim, c(1,3), function(data) {
    max(sapply(1:51, function(i) sum(data[i:(i+39)]))) })
wetLongSim <- t(wetLongSim) * cmPerInch

wetAmtSim <- sapply(spellsSimFull[[i]], function(x) {
    sapply(x, function(y) mean(y$wetAmt[ , 'wetAmt']))
})
wetAmtSim <- wetAmtSim * cmPerInch

wetNumSim <- sapply(spellsSimFull[[i]], function(x) {
    sapply(x, function(y) nrow(y$wetLen))
})

trendsLongSim <- list(); length(trendsLongSim) <- length(timeInts)
trendsAmtSim <- list(); length(trendsAmtSim) <- length(timeInts)
trendsNumSim <- list(); length(trendsNumSim) <- length(timeInts)

for(j in seq_along(timeInts)) {
    rg <- match(timeInts[[j]], yrs)
    trendsLongSim[[j]] <- apply(wetLongSim[ , rg], 1, analyze_trend)
    trendsAmtSim[[j]] <- apply(wetAmtSim[ , rg], 1, analyze_trend)
    trendsNumSim[[j]] <- apply(wetNumSim[ , rg], 1, analyze_trend)
}

output_fn <- paste('results', data_fn, season,
                   ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, sep = '-')
output_fn <- file.path(output_dir, paste0(output_fn, '.Rda'))

save(trendsLongImp, pvalsLongImp, statsLongImp,
     trendsAmtImp, pvalsAmtImp, statsAmtImp,
     trendsNumImp, pvalsNumImp, statsNumImp,
     trendsLongSim, trendsAmtSim, trendsNumSim, file = output_fn)



yrBlocks <- sapply(timeInts, function(x) paste(min(x),max(x), sep = '-'))
sz <- 2.5
szPt <- 0.5

dfImp <- data.frame(x = yrBlocks, y = 10*sapply(trendsNumImp, mean), significance = pvalsNumImp < 0.05,
                    pval = paste0("p = ", sapply(pvalsNumImp, signif, 2)))
dfSim <- data.frame(x = rep(yrBlocks, each = nIts), y = 10*unlist(trendsNumSim))
pNum <- ggplot(dfSim, aes(factor(x), y)) + geom_violin() + ylab("Sen's slope (num/decade)") + xlab('') + ylim(c(-3.5,3.5)) +
    geom_point(data = dfImp, aes(factor(x), y, shape = significance, size = szPt, color = "red")) + scale_shape_manual(values = c(21,16)) +
    theme_minimal() + ggtitle("Number of events") + guides(shape = "none", size = "none", color = "none") +
    geom_text(data = dfImp, aes(factor(x), y =  -3, label = pval), size = sz, color = "red")

dfImp <- data.frame(x = yrBlocks, y = 10*sapply(trendsAmtImp, mean), significance = pvalsAmtImp < 0.05,
                    pval = paste0("p = ", sapply(pvalsAmtImp, signif, 2)))
dfSim <- data.frame(x = rep(yrBlocks, each = nIts), y = 10*unlist(trendsAmtSim))
pAmt <- ggplot(dfSim, aes(factor(x), y)) + geom_violin() + ylab("Sen's slope (cm/decade)") + xlab('') + ylim(c(-1,1)) +
    geom_point(data = dfImp, aes(factor(x), y, shape = significance, size = szPt, color = "red")) + scale_shape_manual(values = c(21,16)) +
    theme_minimal() + ggtitle("Mean wet spell precipitation") + guides(shape = "none", size = "none", color = "none") + 
    geom_text(data = dfImp, aes(factor(x), y =  -1, label = pval), size = sz, color = "red")

dfImp <- data.frame(x = yrBlocks, y = 10*sapply(trendsLongImp, mean), significance = pvalsLongImp < 0.05,
                    pval = paste0("p = ", sapply(pvalsLongImp, signif, 2)))
dfSim <- data.frame(x = rep(yrBlocks, each = nIts), y = 10*unlist(trendsLongSim))
pLong <- ggplot(dfSim, aes(factor(x), y)) + geom_violin() + ylab("Sen's slope (cm/decade)") + xlab('') + ylim(c(-10,10)) +
    geom_point(data = dfImp, aes(factor(x), y, shape = significance, size = szPt, color = "red")) + scale_shape_manual(values = c(21,16)) +
    theme_minimal() + ggtitle("Max. 40-day precipitation") + guides(shape = "none", size = "none", color = "none") + 
    geom_text(data = dfImp, aes(factor(x), y =  -10, label = pval), size = sz, color = "red")
             
pdf(file.path(plotDir, 'inference_quincy.pdf'), width = 8.5, height = 3)
gridExtra::grid.arrange(pNum, pAmt, pLong, nrow = 1, ncol = 3)
dev.off()


