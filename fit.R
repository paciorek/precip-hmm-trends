## version of full s-t code that uses jagam to directly create a single s-t main effect + interaction 

if((exists('build_only') && build_only) || (exists('model_only') && model_only)) compile <- FALSE else compile <- TRUE

output_fn <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, replicate, sep = '-')
output_fn <- file.path(output_dir, paste0(output_fn, '.Rda'))

## test save
if(!file.exists(output_fn)) {
    cat("Saving to file: ", output_fn, "\n")
    save(data_fn, season, trend, config_model, config_holdout,
       niter, replicate, file = output_fn)
}

if(!exists('replicate'))
    seed <- 1 else seed <- replicate
set.seed(seed)

D=3 # Number of clone dry states.
W=2 # Number of wet states.

if(dens_type == 3) INCLUDE_KAPPA <- TRUE else INCLUDE_KAPPA <- FALSE
if(dens_type == 2) GAMMA_MODEL <- TRUE else GAMMA_MODEL <- FALSE

library(nimble) # , lib.loc='/tmp/nim-test')  # allows return of MCMC history with more than 10 params, I think
print(packageVersion('nimble'))

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)


source('nimble_functions.R')

source('load_data.R')
source('set_holdout.R')

## Spline basis matrices created on different machines can have sign
## differences in eigenvectors, so load in pre-calculated values.
load(file.path(data_dir, "basis_matrices.Rda"))
X1 <- X1[[season]]
X2 <- X2[[season]]
K1 <- ncol(X1)
K2 <- ncol(X2)

IMPUTE <- impute_missing

logitScale <- 1.612

source("model_code.R")

pD_SPLINE  <- pD_SEAS || pD_TIME
pDW_SPLINE  <- pDW_SEAS || pDW_TIME

pW_SPLINE  <- pW_SEAS | pW_TIME
pWW_SPLINE  <- pWW_SEAS | pWW_TIME
pWD_SPLINE  <- pWD_SEAS | pWD_TIME

if(!exists('max_sigma'))
    max_sigma <- 100

if(!exists('beta_sd'))
    beta_sd <- 1e6

if(!exists('block_ints'))
    block_ints <- FALSE

if(!exists('constrain_ints'))
    constrain_ints <- FALSE

if(!exists('use_slice'))
    use_slice <- FALSE

if(!exists('constrain_full_means'))
    constrain_full_means <- FALSE

if(!exists('constrain_means_max'))
    constrain_means_max <- FALSE

rain_constants <- list(r = r, n = n, nT = nT, D = D, W = W, prior = rep(1,D+W),
                       K1 = K1, K2 = K2, X1 = X1, X2 = X2, 
                       zeros = rep(0, max(K1, K2)), logitScale = logitScale,
                       pD_SEAS = pD_SEAS, pDW_SEAS = pDW_SEAS, 
                       pW_SEAS = pW_SEAS, pWW_SEAS = pWW_SEAS, pWD_SEAS = pWD_SEAS, 
                       any_pW_SEAS = any_pW_SEAS, any_pWW_SEAS = any_pWW_SEAS, any_pWD_SEAS = any_pWD_SEAS,
                       pi_SEAS = pi_SEAS, sigma_SEAS = sigma_SEAS, xi_SEAS = xi_SEAS,
                       any_pi_SEAS = any_pi_SEAS, any_sigma_SEAS = any_sigma_SEAS, any_xi_SEAS = any_xi_SEAS,
                       pD_TIME = pD_TIME, pDW_TIME = pDW_TIME, 
                       pW_TIME = pW_TIME, pWW_TIME = pWW_TIME, pWD_TIME = pWD_TIME, 
                       any_pW_TIME = any_pW_TIME, any_pWW_TIME = any_pWW_TIME, any_pWD_TIME = any_pWD_TIME,
                       pi_TIME = pi_TIME, sigma_TIME = sigma_TIME, xi_TIME = xi_TIME,
                       any_pi_TIME = any_pi_TIME, any_sigma_TIME = any_sigma_TIME, any_xi_TIME = any_xi_TIME,
                       pD_SPLINE = pD_SPLINE, pDW_SPLINE = pDW_SPLINE,
                       missingness = missing, rounded = rounded, dens_type = dens_type,
                       INCLUDE_KAPPA = INCLUDE_KAPPA, holdout = holdout, doy = doy, max_sigma = max_sigma,
                       CONSTRAIN_FULL_MEANS = constrain_full_means, CONSTRAIN_MEANS_MAX = constrain_means_max)

rain_data <- list(proxy_r = rep(1, nT), constraint=1, constraintWD = rep(1,2))
if(constrain_means_max) {
    rain_data$constraint_max = 1
    rain_data$max_r = max(r)
}

if(impute_missing) 
    rain_constants <- c(rain_constants, list(n_missing = n_missing))

    
## 1-row matrices are converted to vectors and cause dimension mismatches in function calls.
## The code ignores the 2nd row, so this is fine.
rain_dims <- list(pD = c(nT, 1 ,D), pDW = c(nT, 1), pW = c(nT, 1, W), pWW = c(nT, 1, W), pWD = c(nT, 1, D-1),
                   pi = c(nT, n, W+1), sigma = c(nT, n, W+1), xi = c(nT, n, W+1))

if(pD_SPLINE) rain_dims$pD <- c(nT, n, D)
if(pDW_SPLINE) rain_dims$pDW <- c(nT, n)
if(any_pW_SEAS || any_pW_TIME) rain_dims$pW <- c(nT, n, W)
if(any_pWW_SEAS || any_pWW_TIME) rain_dims$pWW <- c(nT, n, W)
if(any_pWD_SEAS || any_pWD_TIME) rain_dims$pWD <- c(nT, n, D-1)

source("initialize.R")
rain_model <- nimbleModel(rain_code, rain_constants, rain_data, rain_inits, dimensions = rain_dims)

while(!is.finite(rain_model$calculate())) {
    seed <- seed + 20
    cat("Trying seed equal to ", seed, "\n")
    source("initialize.R")
    rain_model <- nimbleModel(rain_code, rain_constants, rain_data, rain_inits, dimensions = rain_dims)
    if(seed > 400) stop("Cannot initialize despite multiple tries.")
}

#if(compile)
#    rain_model_compiled <- compileNimble(rain_model)

if(exists('model_only') && model_only)
    stop()

if(exists('do_load')) {
    load(output_fn)
    stop()
}


## Configure the MCMC.
monitors <- c('logProb_proxy_r', 'P_zero','lpD', 'lpDW', 'lpW', 'lpWW', 'lpWD', 'eta', 'alpha', 'gamma')
if(INCLUDE_KAPPA) monitors <- c(monitors, 'phi')

if(pD_SEAS) monitors <- c(monitors, 'beta1_pD', 'beta1_pD_sd')
if(pDW_SEAS) monitors <- c(monitors, 'beta1_pDW', 'beta1_pDW_sd')
if(any_pW_SEAS) monitors <- c(monitors, 'beta1_pW', 'beta1_pW_sd')
if(any_pWW_SEAS) monitors <- c(monitors, 'beta1_pWW', 'beta1_pWW_sd')
if(any_pWD_SEAS) monitors <- c(monitors, 'beta1_pWD', 'beta1_pWD_sd')
if(any_pi_SEAS) monitors <- c(monitors, 'beta1_eta', 'beta1_eta_sd')
if(any_sigma_SEAS) monitors <- c(monitors, 'beta1_alpha', 'beta1_alpha_sd')
if(any_xi_SEAS) monitors <- c(monitors, 'beta1_gamma', 'beta1_gamma_sd')

if(pD_TIME) monitors <- c(monitors, 'beta2_pD', 'beta2_pD_sd')
if(pDW_TIME) monitors <- c(monitors, 'beta2_pDW', 'beta2_pDW_sd')
if(any_pW_TIME) monitors <- c(monitors, 'beta2_pW', 'beta2_pW_sd')
if(any_pWW_TIME) monitors <- c(monitors, 'beta2_pWW', 'beta2_pWW_sd')
if(any_pWD_TIME) monitors <- c(monitors, 'beta2_pWD', 'beta2_pWD_sd')
if(any_pi_TIME) monitors <- c(monitors, 'beta2_eta', 'beta2_eta_sd')
if(any_sigma_TIME) monitors <- c(monitors, 'beta2_alpha', 'beta2_alpha_sd')
if(any_xi_TIME) monitors <- c(monitors, 'beta2_gamma', 'beta2_gamma_sd')

monitors <- monitors[monitors != '']

if(is.function(scale)) scale <- 0.1
if(!exists('adaptFactorExponent')) adaptFactorExponent <- 0.25
if(!exists('adaptInterval')) adaptInterval <- 50

if(!exists('use_new_blockRW') || !use_new_blockRW) {
    control <- list(scale = 0.01, adaptInterval = 50)
} else control = list(scale = scale, adaptFactorExponent = adaptFactorExponent, adaptInterval = adaptInterval)

if(!exists('thin'))
    thin <- 10
if(!exists('thin2'))
    thin2 <- 100
rain_mcmc_conf <- configureMCMC(rain_model, monitors=monitors, thin=thin, thin2 = thin2,
                                control = control, enableWAIC = TRUE)

if(impute_missing) {
    rain_mcmc_conf$removeSamplers('imputed')
    rain_mcmc_conf$addSampler(paste0('imputed[1:', n_missing, ']'), 'impute_ffbs',
                              control = list(n_missing_by_year = rowSums(missing), missing = missing,
                                             D = D, W = W, rounded = rounded,
                                             pD_TIME = pD_SPLINE, pDW_TIME = pDW_SPLINE,
                                             pW_TIME = pW_SPLINE, pWW_TIME = pWW_SPLINE,
                                             pWD_TIME = pWD_SPLINE, thin = thin2, dens_type = dens_type))
    rain_mcmc_conf$addMonitors2('imputed')
}

## For now fit seasonal splines unconstrained
if(FALSE) {
    if(pD_SEAS) rain_mcmc_conf$removeSamplers('beta1_pD_sd')
    if(pDW_SEAS) rain_mcmc_conf$removeSamplers('beta1_pDW_sd')
    if(any_pW_SEAS) rain_mcmc_conf$removeSamplers('beta1_pW_sd')
    if(any_pWW_SEAS) rain_mcmc_conf$removeSamplers('beta1_pWW_sd')
    if(any_pWD_SEAS) rain_mcmc_conf$removeSamplers('beta1_pWD_sd')
    if(any_pi_SEAS) rain_mcmc_conf$removeSamplers('beta1_eta_sd')
    if(any_sigma_SEAS) rain_mcmc_conf$removeSamplers('beta1_alpha_sd')
    if(any_xi_SEAS) rain_mcmc_conf$removeSamplers('beta1_gamma_sd')
}

samplers <- rain_mcmc_conf$getSamplers()

for(sampler in samplers) {
    target <- sampler$target
    if(length(grep("_sd", target))) {
        rain_mcmc_conf$removeSamplers(target)
        rain_mcmc_conf$addSampler(target, type = 'RW', control = c(control, log = TRUE))
    }}


if(any(pW_SEAS) && !all(pW_SEAS))
  for(i in seq_along(pW_SEAS))
    if(!pW_SEAS[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta1_pW_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta1_pW[1:', K1, ', ', i, ']'))
    }
if(any(pWW_SEAS) && !all(pWW_SEAS))
  for(i in seq_along(pWW_SEAS))
    if(!pWW_SEAS[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta1_pWW_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta1_pWW[1:', K1, ', ', i, ']'))
    }
if(any(pWD_SEAS) && !all(pWD_SEAS))
  for(i in seq_along(pWD_SEAS))
    if(!pWD_SEAS[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta1_pWD_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta1_pWD[1:', K1, ', ', i, ']'))
    }
if(any(pi_SEAS) && !all(pi_SEAS))
  for(i in seq_along(pi_SEAS))
    if(!pi_SEAS[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta1_eta_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta1_eta[1:', K1, ', ', i, ']'))
    }
if(any(sigma_SEAS) && !all(sigma_SEAS))
  for(i in seq_along(sigma_SEAS))
    if(!sigma_SEAS[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta1_alpha_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta1_alpha[1:', K1, ', ', i, ']'))
    }
if(any(xi_SEAS) && !all(xi_SEAS))
  for(i in seq_along(xi_SEAS))
    if(!xi_SEAS[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta1_gamma_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta1_gamma[1:', K1, ', ', i, ']'))
    }

  
if(any(pW_TIME) && !all(pW_TIME))
  for(i in seq_along(pW_TIME))
    if(!pW_TIME[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta2_pW_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta2_pW[1:', K2, ', ', i, ']'))
    }
if(any(pWW_TIME) && !all(pWW_TIME))
  for(i in seq_along(pWW_TIME))
    if(!pWW_TIME[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta2_pWW_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta2_pWW[1:', K2, ', ', i, ']'))
    }
if(any(pWD_TIME) && !all(pWD_TIME))
  for(i in seq_along(pWD_TIME))
    if(!pWD_TIME[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta2_pWD_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta2_pWD[1:', K2, ', ', i, ']'))
    }
if(any(pi_TIME) && !all(pi_TIME))
  for(i in seq_along(pi_TIME))
    if(!pi_TIME[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta2_eta_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta2_eta[1:', K2, ', ', i, ']'))
    }
if(any(sigma_TIME) && !all(sigma_TIME))
  for(i in seq_along(sigma_TIME))
    if(!sigma_TIME[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta2_alpha_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta2_alpha[1:', K2, ', ', i, ']'))
    }
if(any(xi_TIME) && !all(xi_TIME))
  for(i in seq_along(xi_TIME))
    if(!xi_TIME[i]) {
      rain_mcmc_conf$removeSamplers(paste0('beta2_gamma_sd[', i, ']'))
      rain_mcmc_conf$removeSamplers(paste0('beta2_gamma[1:', K2, ', ', i, ']'))
    }


if(exists('block_ints') && block_ints) {
    ints <- c(paste0('lpD[', 1:D, ']'), 'lpDW', paste0('lpW[', 1:W, ']'),
                 paste0('lpWW[', 1:W, ']'), paste0('lpWD[', 1:(D-1), ']'),
                 paste0('eta[', 1:(W+1), ']'), paste0('alpha[', 1:(W+1), ']'),
                 paste0('gamma[', 1:(W+1), ']'))
    for(nm in ints)
        rain_mcmc_conf$removeSamplers(nm)
    rain_mcmc_conf$addSampler(ints, 'RW_block', control = control)
}


if(use_slice) {
    samplers <- rain_mcmc_conf$getSamplers()
    for(i in seq_along(samplers)) 
        if(samplers[[i]]$name == 'RW') {
            rain_mcmc_conf$removeSampler(samplers[[i]]$target)
            rain_mcmc_conf$addSampler(samplers[[i]]$target, 'slice')
        }
}
        
stop()
rain_mcmc <- buildMCMC(rain_mcmc_conf)
if(compile)
    rain_mcmc_compiled <- compileNimble(rain_mcmc, project = rain_model)

print(rain_mcmc_conf)

if(exists('build_only') && build_only)
    stop()

print(system.time(out <- runMCMC(rain_mcmc_compiled, nburnin = 0, niter = niter, WAIC = TRUE))) # nburnin=1000, niter=11000)))

## Note that NA warnings when start MCMC are probably because when some terms are not time-varying, we have vectors of length 2 where the 2nd element is NA by design.

samples <- out$samples
waic <- out$WAIC
if(impute_missing) {
    imputations <- as.matrix(rain_mcmc_compiled$mvSamples2)
} else imputations <- NULL

if(!exists('nburnin'))
    nburnin <- niter*burnin_waic_frac
waic_burned <- calculateWAIC(rain_mcmc_compiled, nburnin = nburnin)

cat("WAIC: \n")
print(waic)
cat("WAIC (burned): \n")
print(waic_burned)

calc_negLogLik <- function(compiled_model, samples) {
    m  <- nrow(samples)
    negLogLik <- rep(0, m)
    nm <- dimnames(samples)[[2]]
    varnm  <- gsub("\\[.*", "", nm)
    vars <- unique(varnm)
    for(i in seq_len(m)) {
        if(i %% 100 == 0) print(i)
        for(v in vars)
            compiled_model[[v]] <- samples[i, which(v == varnm)]
        compiled_model$calculate()
        ## for the moment can only get entire log-likelihood
        negLogLik[i] <- -compiled_model$getLogProb('proxy_r')
    }
    return(negLogLik)
}


n <- nrow(samples)
secondhalf <- (n/2 + 1):n

ll_train <- calc_negLogLik(rain_model_compiled, samples)
cat("Training ll: ", mean(ll_train[secondhalf]), "\n")

if(holdout_style != 'none') {
    ## for holdout cases
    rain_model$test  <- rain_model_compiled$test <- 1
    rain_model$one_step  <- rain_model_compiled$one_step <- 1
    ll_test1 <- calc_negLogLik(rain_model_compiled, samples)
    cat("One step test ll: ", mean(ll_test1[secondhalf], na.rm = TRUE), "\n")
    cat("Missing values: ", sum(is.na(ll_test1[secondhalf])), "\n")
    
    rain_model$one_step  <- rain_model_compiled$one_step <- 0
    ll_test2 <- calc_negLogLik(rain_model_compiled, samples)
    cat("Multi-step test ll: ", mean(ll_test2[secondhalf], na.rm = TRUE), "\n")
    cat("Missing values: ", sum(is.na(ll_test2[secondhalf])), "\n")
} else {
    ll_test1 <- NA
    ll_test2 <- NA
}

cat("Output file is ", output_fn, "\n")

save(r, missing, doy, dayCovar,
     samples, imputations, waic, waic_burned, ll_train, ll_test1, ll_test2,
     dens_type, rounded, allow_missing, use_inches, use_splines, use_halfflat, use_slice,
     IMPUTE, GAMMA_MODEL, constrain_full_means, constrain_ints, 
     pD_SEAS, pDW_SEAS, pW_SEAS, pWW_SEAS, pWD_SEAS, pi_SEAS, sigma_SEAS, xi_SEAS,
     pD_TIME, pDW_TIME, pW_TIME, pWW_TIME, pWD_TIME, pi_TIME, sigma_TIME, xi_TIME,
     K1, K2, X1, X2,
     rain_code, rain_constants, rain_data, rain_inits, rain_dims,
     beta_sd, use_new_blockRW, block_ints, 
     data_fn, season, trend, config_model, config_holdout, holdout,
     niter, replicate,
     file = output_fn)

