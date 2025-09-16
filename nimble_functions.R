forAlg2D <- nimbleFunction(run=function(p_init=double(1), pD = double(2), pDW = double(1), pW = double(2), pWW = double(2), pWD = double(2), 
                                     pD_time = double(0), pDW_time = double(0), pW_time = double(1), pWW_time = double(1), pWD_time = double(1),
                                     n = double(0), D = double(0), W = double(0),
                                     dens=double(2),missing=double(1), holdout = double(1), one_step = integer(0), test = integer(0)){
  P <- nimMatrix(0, D+W, D+W)
  Pnew <- rep(0, D+W)
  c <- numeric(n)
  if(!test) {
      if(missing[1] | holdout[1]) {
          alpha <- p_init
          c[1] <- 1
      } else {
          c[1] <- sum(dens[1,]*p_init)
          alpha <- (dens[1,]*p_init)/c[1]
      }
  } else {
      if(missing[1]) {
          alpha <- p_init
          c[1] <- 1
      } else {
          if(holdout[1]) {
              c[1] <- sum(dens[1,]*p_init)  # calc logLik from test points
              if(one_step) {
                  alpha <- (dens[1,]*p_init)/c[1] # learn from data
              } else alpha <- p_init # don't learn from data
          } else {
              c[1] <- sum(dens[1,]*p_init)
              alpha <- (dens[1,]*p_init)/c[1]
              c[1] <- 1  # don't include logLik from training points
          }
      }
  }
  idx_pD <- 1; idx_pDW <- 1; idx_pW <- rep(1, W); idx_pWW <- rep(1, W); idx_pWD <- rep(1, D-1)
  for(t in 2:n){
    for(j in 1:D)
      P[j, j] <- pD[idx_pD, j]
    P[1:D, 4] <- (1-pD[idx_pD,])*pDW[idx_pDW]
    P[1:D, 5]<- (1-pD[idx_pD,])*(1-pDW[idx_pDW])
    P[4, 4] <- pW[idx_pW[1], 1]
    P[5, 5] <- pW[idx_pW[2], 2]
    P[4, 5] <- pWW[idx_pWW[1], 1] * (1-P[4, 4])
    P[5, 4] <- pWW[idx_pWW[2], 2] * (1-P[5, 5])
    tmp <- 1-P[4:5,4]-P[4:5,5]
    tmp[tmp < 0] <- 0  # avoid negative probs from numerical imprecision
    P[4:5, 1] <- pWD[idx_pWD[1], 1] * tmp
    tmp <- tmp - P[4:5, 1]
    P[4:5, 2] <- pWD[idx_pWD[2], 2] * tmp
    P[4:5, 3] <- tmp - P[4:5, 2]
    
    for(j in 1:(D+W)){
      Pnew[j] <- sum(alpha*P[,j])
    }
    ## CHECK THIS
    if(!test) {
        if(missing[t] | holdout[t]) {
            delta <- Pnew
            c[t] <- 1
        } else {
            delta <- dens[t,]*Pnew
            c[t] <- sum(delta)
        }
        alpha <- delta/c[t]
    } else {
         if(missing[t]) {
            delta <- Pnew
            c[t] <- 1  # no log-lik
            alpha <- Pnew # can't learn from data
         } else {
             if(holdout[t]) {
                 delta <- dens[t,]*Pnew
                 c[t] <- sum(delta) # calc logLik from test points
                 if(one_step) {
                     alpha <- delta/c[t] # learn from data
                 } else alpha <- Pnew # don't learn from data
             } else {
                 delta <- dens[t,]*Pnew
                 c[t] <- 1  # don't include logLik from training points
                 alpha <- delta/sum(delta) # learn from data
             }
         }
    }        
    idx_pD <- idx_pD + pD_time
    idx_pDW <- idx_pDW + pDW_time
    idx_pW <- idx_pW + pW_time
    idx_pWW <- idx_pWW + pWW_time
    idx_pWD <- idx_pWD + pWD_time
  }
  returnType(double(0))
  return(sum(log(c)))
})

pow_vec <- nimbleFunction(run = function(x=double(1),y=double(1)){
  returnType(double(1))
  return(exp(y*log(x)))
})

dzigpd <- nimbleFunction(run = function(x=double(1),p=double(1),scale=double(1),shape=double(1),
                                      rounded=double(0),log=integer(0)){
  returnType(double(1))
  n <- length(x)
  imis <- x<0
  ix <- x>0
  y <- p
  y[imis] <- 1
  z <- 1 + shape[ix]*x[ix]/scale[ix]
  if(any(z < 0)) {
    y[ix] <- 0
  } else {
    if(rounded) {
      y[ix] <- (1-p[ix]) * ( pow_vec(1 + shape[ix]*(x[ix]-.005)/scale[ix], -1/shape[ix]) - 
                               pow_vec(1 + shape[ix]*(x[ix]+.005)/scale[ix], -1/shape[ix]) ) / 
                              pow_vec(1 + shape[ix]*.005/scale[ix], -1/shape[ix])
    } else y[ix] <- (1-p[ix])*(1/scale[ix])*pow_vec(1 + x[ix]*shape[ix]/scale[ix], -1/shape[ix]-1)
  }
  if(log) y <- log(y)
  return(y)
})

dzigamma <- nimbleFunction(run=function(x=double(1),p=double(1),scale=double(1),shape=double(1),
                                        rounded=double(0),log=integer(0)){
  returnType(double(1))
  n <- length(x)
  imis <- x<0
  ix <- x>0
  y <- p
  y[imis]=1
  if(sum(ix) > 0) {  # otherwise get mis-matched dimensions run-time warning; why not needed for dzigpd?
      if(rounded) {
          y[ix] <- (1-p[ix])*(pgamma(x[ix]+0.005, shape = shape[ix], scale = scale[ix]) -
                              pgamma(x[ix]-0.005, shape = shape[ix], scale = scale[ix])) /
              pgamma(0.005, shape = shape[ix], scale = scale[ix], lower.tail = FALSE)
      } else y[ix] <- (1-p[ix])*dgamma(x[ix], shape = shape[ix], scale = scale[ix])
  }
  if(log) y <- log(y)
  return(y)
})

dzinav <- nimbleFunction(run=function(x=double(1),p=double(1),scale=double(1),shape=double(1),kappa=double(0),
                                      rounded=double(0),log=integer(0)){
  returnType(double(1))
  n <- length(x)
  imis <- x<0
  ix <- x>0
  y <- p
  y[imis]=1
  z <- 1 + shape[ix]*x[ix]/scale[ix]
  if(any(z < 0)) {
    y[ix] <- 0
  } else {
    if(rounded) {
      y[ix] <- (1-p[ix]) * (pow(1-pow_vec(1+shape[ix]*(x[ix]+.005)/scale[ix], -1/shape[ix]), kappa) -
                          pow(1-pow_vec(1+shape[ix]*(x[ix]-.005)/scale[ix], -1/shape[ix]), kappa)) /
        (1 - pow(1-pow_vec(1+shape[ix]*.005/scale[ix], -1/shape[ix]),kappa))
    } else y[ix] <- (1-p[ix])*kappa*pow(1-pow_vec(z,-1/shape[ix]), kappa - 1)*pow_vec(z, -1/shape[ix]-1)/scale[ix]
  }
  if(log) y <- log(y)
  return(y)
})


dproxy <- nimbleFunction(
    run = function(x = double(0), dens = double(2), N = double(0), D = double(0), W = double(0),P_zero=double(1),
                 pD = double(2), pD_time = double(0), 
                 pDW = double(1), pDW_time = double(0), 
                 pW = double(2), pW_time = double(1), 
                 pWW = double(2), pWW_time = double(1), 
                 pWD = double(2), pWD_time = double(1), missing = double(1), 
                 holdout = double(1), one_step = integer(0), test = integer(0), log = integer(0)) {
                     
        returnType(double(0))
        return(forAlg2D(P_zero, pD, pDW, pW, pWW, pWD, 
                        pD_time, pDW_time, pW_time, pWW_time, pWD_time,
                        N, D, W, dens, missing, holdout, one_step, test))
    }
)

calc_dens <- nimbleFunction(
  run = function(x=double(1), N = double(0), D = double(0), W = double(0), 
                 pi = double(2), 
                 sigma = double(2), 
                 xi = double(2), 
                 kappa = double(1), include_kappa = double(0), 
                 rounded = double(0), dens_type = double(0)) {
    
    dens = matrix(nrow=N, ncol= D+W)

    idx <- c(1, D+1, D+2)  # this assumes W=2
    for(j in 1:(W+1)){
      if(dens_type == 1) {
        dens[, idx[j]] = dzigpd(x, pi[ , j], sigma[ , j], xi[ , j], rounded, log=FALSE)
      } else if(dens_type == 2) {
        dens[, idx[j]] = dzigamma(x, pi[ , j], sigma[ , j], xi[ , j], rounded, log=FALSE)
      } else dens[, idx[j]] = dzinav(x, pi[ , j], sigma[ , j], xi[ , j], kappa[j], rounded, log=FALSE)
    }
    for(j in 2:D)
      dens[ , j] <- dens[, 1]
    
      returnType(double(2))
      return(dens)
  }
)


drain <- nimbleFunction(
  run = function(x=double(1), N = double(0), D = double(0), W = double(0), P_zero=double(1),
                 pD = double(2), pD_time = double(0), 
                 pDW = double(1), pDW_time = double(0), 
                 pW = double(2), pW_time = double(1), 
                 pWW = double(2), pWW_time = double(1), 
                 pWD = double(2), pWD_time = double(1), 
                 pi = double(2), 
                 sigma = double(2), 
                 xi = double(2), 
                 kappa = double(1), include_kappa = double(0), 
                 missing=double(1), rounded = double(0), dens_type = double(0), holdout = double(1), one_step = integer(0), test = integer(0), log = integer(0)) {
    
    dens = matrix(nrow=N, ncol= D+W)

    idx <- c(1, D+1, D+2)  # this assumes W=2
    for(j in 1:(W+1)){
      if(dens_type == 1) {
        dens[, idx[j]] = dzigpd(x, pi[ , j], sigma[ , j], xi[ , j], rounded, log=FALSE)
      } else if(dens_type == 2) {
        dens[, idx[j]] = dzigamma(x, pi[ , j], sigma[ , j], xi[ , j], rounded, log=FALSE)
      } else dens[, idx[j]] = dzinav(x, pi[ , j], sigma[ , j], xi[ , j], kappa[j], rounded, log=FALSE)
    }
    for(j in 2:D)
      dens[ , j] <- dens[, 1]
    
    returnType(double(0))
     return(forAlg2D(P_zero, pD, pDW, pW, pWW, pWD, 
                    pD_time, pDW_time, pW_time, pWW_time, pWD_time,
                    N, D, W, dens, missing, holdout, one_step, test))
  }
)


drain_imputed <- nimbleFunction(
  run = function(x=double(1), log = integer(0)) {    
      returnType(double(0))
      if(log) return(0) else return(1)
  }
)


sampler_impute_ffbs <- nimbleFunction(
    name = 'sampler_impute_ffbs',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        n_missing_by_year <- control$n_missing_by_year
        n_missing <- sum(n_missing_by_year)
        missingness <- control$missing  # weird compilation failure if name variable as 'missing'
        D <- control$D
        W <- control$W
        rounded <- control$rounded
        pD_TIME <- control$pD_TIME
        pDW_TIME <- control$pDW_TIME
        pW_TIME <- control$pW_TIME
        pWW_TIME <- control$pWW_TIME
        pWD_TIME <- control$pWD_TIME
        thin <- control$thin
        dens_type <- control$dens_type
        nT <- nrow(model$r)
        n <- ncol(model$r)
        idx_thin <- 1
        idx <- c(1, D+1, D+2)  # this assumes W=2
        alpha <- matrix(0, D+W, n)
        P <- array(0, c(D+W, D+W, n))  # P[,,t] is probability for transitioning from time t to time t+1
        timesRan <- 0
    },
    run = function() {
        ## This combines density calculation per drain and forward algorithm per forAlg2D
        timesRan <<- timesRan + 1
        if(idx_thin %% thin == 0) {
            idx_miss <- 1
            for(tt in 1:nT) {
                if(n_missing_by_year[tt] > 0) {
                    ## Forward filtering
                    Pnew <- rep(0, D+W)
                    if(missingness[tt, 1]) {
                        alpha[ , 1] <<- model[['P_zero']]
                        if(idx_miss > n_missing) stop("idx_miss out of bounds")
                        idx_miss <- idx_miss + 1
                    } else {
                        tmp1 <- model[['dens']][tt, 1,]*model[['P_zero']]
                        alpha[ , 1] <<- tmp1/sum(tmp1)
                    }
                    idx_pD <- 1; idx_pDW <- 1; idx_pW <- rep(1, W); idx_pWW <- rep(1, W); idx_pWD <- rep(1, D-1)
                    for(t in 2:n){
                        for(j in 1:D)
                            P[j, j, t-1] <<- model[['pD']][tt, idx_pD, j]
                        P[1:D, 4, t-1] <<- (1-model[['pD']][tt, idx_pD,])*model[['pDW']][tt, idx_pDW]
                        P[1:D, 5, t-1] <<- (1-model[['pD']][tt, idx_pD,])*(1-model[['pDW']][tt, idx_pDW])
                        P[4, 4, t-1] <<- model[['pW']][tt, idx_pW[1], 1]
                        P[5, 5, t-1] <<- model[['pW']][tt, idx_pW[2], 2]
                        P[4, 5, t-1] <<- model[['pWW']][tt, idx_pWW[1], 1] * (1-P[4, 4, t-1])
                        P[5, 4, t-1] <<- model[['pWW']][tt, idx_pWW[2], 2] * (1-P[5, 5, t-1])
                        tmp2 <- 1-P[4:5,4, t-1]-P[4:5,5, t-1]
                        ## if(tmp2[1] < 0 | tmp2[2] < 0) 
                        ##    nimPrint(P[4,4, t-1], " ", P[5,4, t-1], " ",
                        ##             P[4,5, t-1], " ", P[5,5, t-1], " ",
                        # #            tmp2[1], " ", tmp2[2])
                        tmp2[tmp2 < 0] <- 0  # numerical imprecision adjustment to avoid rcat -> NaN
                        P[4:5, 1, t-1] <<- model[['pWD']][tt, idx_pWD[1], 1] * tmp2
                        tmp2 <- tmp2 - P[4:5, 1, t-1]
                        P[4:5, 2, t-1] <<- model[['pWD']][tt, idx_pWD[2], 2] * tmp2
                        P[4:5, 3, t-1] <<- tmp2 - P[4:5, 2, t-1]
                        for(j in 1:(D+W)){
                            Pnew[j] <- sum(alpha[ , t-1]*P[, j, t-1])
                        }
                        if(missingness[tt, t]) {
                            alpha[ , t] <<- Pnew
                            if(idx_miss > n_missing) stop("idx_miss out of bounds")
                            idx_miss <- idx_miss + 1
                       } else {
                           tmp1 <- model[['dens']][tt, t, ]*Pnew
                           alpha[ , t] <<- tmp1/sum(tmp1)
                       }
                       idx_pD <- idx_pD + pD_TIME
                       idx_pDW <- idx_pDW + pDW_TIME
                       idx_pW <- idx_pW + pW_TIME
                       idx_pWW <- idx_pWW + pWW_TIME
                       idx_pWD <- idx_pWD + pWD_TIME
                    }
                    idx_miss_current <- idx_miss
                    idx_miss <- idx_miss - 1
                    ## Backward sampling
                    z <- 1  # otherwise get complaint about z not created yet
                    for(trev in 1:n) {  # no reverse indexing
                        t <- n-trev+1
                        if(t == n) { # don't condition on the future
                            probs <- alpha[ , t]
                        } else {
                            probs <- alpha[ , t] * P[ , z, t]
                        }
                        z <- rcat(1, probs)
                        if(missingness[tt, t]) {
                            zReal <- z - D + 1
                            if(zReal < 1)
                                zReal <- 1
                            dry <- rbinom(1, 1, model[['pi']][tt, t, zReal])
                            if(idx_miss < 1 | idx_miss > n_missing) stop("Problem with missingness indexing.")
                            if(dry != 1) {
                                scale_tmp <- model[['sigma']][tt, t, zReal]
                                shape_tmp <- model[['xi']][tt, t, zReal]
                                if(dens_type != 2) {
                                    rainfall <- (scale_tmp/shape_tmp)*
                                        ((1-runif(1)^(1/model[['kappa']][zReal]))^(-shape_tmp)-1)
                                } else rainfall  <- rgamma(1, shape_tmp, scale = scale_tmp)
                                model[[target]][idx_miss] <<- rainfall
                            } else {
                                model[[target]][idx_miss] <<- 0                                
                            }
                            idx_miss <- idx_miss - 1
                        }
                    }
                    idx_miss <- idx_miss_current
                }
            }
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = FALSE)
        idx_thin <<- idx_thin + 1
    },
    methods = list(
        reset = function() {}
    )
)    
    


## This uses the filtering distribution so (inappropriately) does not account for information from future time points.
sampler_impute_filter <- nimbleFunction(
    name = 'sampler_impute_filter',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        n_missing_by_year <- control$n_missing_by_year
        n_missing <- sum(n_missing_by_year)
        missingness <- control$missing  # weird compilation failure if name variable as 'missing'
        D <- control$D
        W <- control$W
        rounded <- control$rounded
        pD_TIME <- control$pD_TIME
        pDW_TIME <- control$pDW_TIME
        pW_TIME <- control$pW_TIME
        pWW_TIME <- control$pWW_TIME
        pWD_TIME <- control$pWD_TIME
        thin <- control$thin
        dens_type <- control$dens_type
        nT <- nrow(model$r)
        n <- ncol(model$r)
        idx_thin <- 1
        idx <- c(1, D+1, D+2)  # this assumes W=2
        alpha <- rep(0, D+W)
    },
    run = function() {
        ## This combines density calculation per drain and forward algorithm per forAlg2D
        if(idx_thin %% thin == 0) {
            idx_miss <- 1
            for(tt in 1:nT) {
                if(n_missing_by_year[tt] > 0) {
                    x <- model[['r']][tt, ]
                    dens <- matrix(nrow = n, ncol = D+W)
                    for(j in 1:(W+1)){
                        if(dens_type == 1) {
                            dens[, idx[j]] <- dzigpd(x, model[['pi']][tt, , j], model[['sigma']][tt, , j], model[['xi']][tt, , j], rounded, log=FALSE)
                        } else if(dens_type == 2) {
                            dens[, idx[j]] <- dzigamma(x, model[['pi']][tt, , j], model[['sigma']][tt, , j], model[['xi']][tt, , j], rounded, log=FALSE)
                        } else dens[, idx[j]] <- dzinav(x, model[['pi']][tt, , j], model[['sigma']][tt, , j], model[['xi']][tt, , j], model[['kappa']][j], rounded, log=FALSE)
                    }
                    for(j in 2:D)
                        dens[ , j] <- dens[, 1]
                    P <-  nimMatrix(0, D+W, D+W)
                    Pnew <- rep(0, D+W)
                    if(missingness[tt, 1]) {
                        if(idx_miss > n_missing) stop("idx_miss out of bounds")
                        z <- rcat(1, model[['P_zero']])
                        alpha <<- rep(0, D+W)
                        alpha[z] <<- 1
                        zReal <- z - D + 1
                        if(zReal < 1)
                            zReal <- 1
                        dry <- rbinom(1, 1, model[['pi']][tt, 1, zReal])
                        if(dry != 1) {
                            scale_tmp <- model[['sigma']][tt, 1, zReal]
                            shape_tmp <- model[['xi']][tt, 1, zReal]
                            if(dens_type != 2) {
                                rainfall <- (scale_tmp/shape_tmp)*
                                    ((1-runif(1)^(1/model[['kappa']][zReal]))^(-shape_tmp)-1)
                            } else rainfall  <- rgamma(1, shape = shape_tmp, scale = scale_tmp)
                            model[[target]][idx_miss] <<- rainfall
                        } else model[[target]][idx_miss] <<- 0
                        idx_miss <- idx_miss + 1
                    } else {
                        tmp1 <- dens[1,]*model[['P_zero']]
                        alpha <<- tmp1/sum(tmp1)
                    }
                    idx_pD <- 1; idx_pDW <- 1; idx_pW <- rep(1, W); idx_pWW <- rep(1, W); idx_pWD <- rep(1, D-1)
                    for(t in 2:n){
                        for(j in 1:D)
                            P[j, j] <- model[['pD']][tt, idx_pD, j]
                        P[1:D, 4] <- (1-model[['pD']][tt, idx_pD,])*model[['pDW']][tt, idx_pDW]
                        P[1:D, 5]<- (1-model[['pD']][tt, idx_pD,])*(1-model[['pDW']][tt, idx_pDW])
                        P[4, 4] <- model[['pW']][tt, idx_pW[1], 1]
                        P[5, 5] <- model[['pW']][tt, idx_pW[2], 2]
                        P[4, 5] <- model[['pWW']][tt, idx_pWW[1], 1] * (1-P[4, 4])
                        P[5, 4] <- model[['pWW']][tt, idx_pWW[2], 2] * (1-P[5, 5])
                        tmp2 <- 1-P[4:5,4]-P[4:5,5]
                        P[4:5, 1] <- model[['pWD']][tt, idx_pWD[1], 1] * tmp2
                        tmp2 <- tmp2 - P[4:5, 1]
                        P[4:5, 2] <- model[['pWD']][tt, idx_pWD[2], 2] * tmp2
                        P[4:5, 3] <- tmp2 - P[4:5, 2]
                        for(j in 1:(D+W)){
                            Pnew[j] <- sum(alpha*P[,j])
                        }
                        if(missingness[tt, t]) {
                            if(idx_miss > n_missing) stop("idx_miss out of bounds")
                            z <- rcat(1, Pnew)
                            alpha <<- rep(0, D+W)
                            alpha[z] <<- 1
                            zReal <- z - D + 1
                            if(zReal < 1)
                                zReal <- 1
                            dry <- rbinom(1, 1, model[['pi']][tt, t, zReal])
                            if(dry != 1) {
                                scale_tmp <- model[['sigma']][tt, t, zReal]
                                shape_tmp <- model[['xi']][tt, t, zReal]
                                if(dens_type != 2) {
                                    rainfall <- (scale_tmp/shape_tmp)*
                                        ((1-runif(1)^(1/model[['kappa']][zReal]))^(-shape_tmp)-1)
                                } else rainfall  <- rgamma(1, shape_tmp, scale = scale_tmp)
                                model[[target]][idx_miss] <<- rainfall
                            } else model[[target]][idx_miss] <<- 0
                            idx_miss <- idx_miss + 1
                       } else {
                           tmp1 <- dens[t,]*Pnew
                           alpha <<- tmp1/sum(tmp1)
                       }
                        idx_pD <- idx_pD + pD_TIME
                        idx_pDW <- idx_pDW + pDW_TIME
                        idx_pW <- idx_pW + pW_TIME
                        idx_pWW <- idx_pWW + pWW_TIME
                        idx_pWD <- idx_pWD + pWD_TIME
                    }
                }
             }
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = FALSE)
        idx_thin <<- idx_thin + 1
    },
    methods = list(
        reset = function() {}
    )
)    
    

linexp <- nimbleFunction( run = function(x = double(1)) {
    returnType(double(1))
    neg <- x < 0
    x[neg] <- exp(x[neg])
    x[!neg] <- 1 + x[!neg]
    return(x)
})
