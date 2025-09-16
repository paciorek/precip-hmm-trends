combine_chains <- function(data_fn, season, trend = 'trend', config_model = 1, config_holdout = 0, niter = 20000, output_dir = 'output', nburnin = 0, thin = 1) {
    output_pattern <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, sep = '-')
    samples_full <- NULL
    imputations_full <- NULL
    files <- list.files(output_dir, output_pattern)
    for(i in seq_along(files)) {
        cat("Loading: ", files[i], ".\n")
        load(file.path(output_dir, files[i]))
        postburn <- (nburnin+1):nrow(samples)
        nr <- length(postburn)
        if(i == 1) {
            samples_full <- matrix(0, nrow = length(files)*nr, ncol = ncol(samples))
            colnames(samples_full) <- colnames(samples)
        }
        samples_full[(1+(i-1)*nr):(nr*i), ] <- samples[postburn, ]
        if(exists('imputations') && !is.null(imputations)) {
            if(i == 1) {
                imputations_full <- matrix(0, nrow = length(files)*nr, ncol = ncol(imputations))
                colnames(imputations_full) <- colnames(imputations)
            }
            imputations_full[(1+(i-1)*nr):(nr*i), ] <- imputations[postburn, ]
        }
    }
    samples <- samples_full[seq(1, nrow(samples_full), by = thin), ]
    if(exists('imputations') && !is.null(imputations))
        imputations <- imputations_full[seq(1, nrow(imputations_full), by = thin), ] else imputations <- NULL
    return(list(samples = samples, imputations = imputations))
}

combine_ll <- function(data_fn, season, trend = 'trend', config_model = 1, config_holdout = 0, niter = 20000, output_dir = 'output', nburnin = 0, thin = 1) {
    output_pattern <- paste('fit', data_fn, season, ifelse(trend == 'trend', 'trend', 'notrend'), config_model, config_holdout, niter, sep = '-')
    ll_train_full <- ll_test1_full <- ll_test2_full <-  NULL
    files <- list.files(output_dir, output_pattern)
    for(i in seq_along(files)) {
        cat("Loading: ", files[i], ".\n")
        load(file.path(output_dir, files[i]))
        postburn <- (nburnin+1):length(ll_train)
        ll_train_full <- c(ll_train_full, ll_train[postburn])
        ll_test1_full <- c(ll_test1_full, ll_test1[postburn])
        ll_test2_full <- c(ll_test2_full, ll_test2[postburn])
    }
    ll_train_full <- ll_train_full[seq(1, length(ll_train_full), by = thin)]
    ll_test1_full <- ll_test1_full[seq(1, length(ll_test1_full), by = thin)]
    ll_test2_full <- ll_test2_full[seq(1, length(ll_test2_full), by = thin)]
    return(list(ll_train = ll_train_full, ll_test1 = ll_test1_full, ll_test2 = ll_test2_full))
}


create_Pmat <- function(pD, pDW, pW, pWW, pWD) {
  D <- length(pD)
  W <- length(pW)
  ## samples by time by matrix
  Pmat <- array(0, c(dim(pDW), D+W, D+W))
  for(j in 1:D) {
      Pmat[ , , , j, j] <- pD[[j]]
      Pmat[ , , , j, 4] <- (1-pD[[j]])*pDW
      Pmat[ , , , j, 5] <- (1-pD[[j]])*(1-pDW)
  }
  Pmat[ , , , 4, 4] <- pW[[1]]
  Pmat[ , , , 5, 5] <- pW[[2]]
  Pmat[ , , , 4, 5] <- pWW[[1]] * (1-Pmat[ , , , 4, 4])
  Pmat[ , , , 5, 4] <- pWW[[2]] * (1-Pmat[ , , , 5, 5])
  tmp <- 1-Pmat[ , , , 4:5,4]-Pmat[ , , , 4:5,5]
  for(j in 1:W) {
      Pmat[ , , , j+3, 1] <- pWD[[1]] * tmp[,,,j]
      tmp[,,,j] <- tmp[,,,j] - Pmat[ , , , j+3, 1]
      Pmat[ , , , j+3, 2] <- pWD[[2]] * tmp[,,,j]
      Pmat[  ,, , j+3, 3] <- tmp[,,,j] - Pmat[ , , , j+3, 2]
  }
  return(Pmat)
}



simulate_vec_rep <- function(P_zero, Pmat, pi_input, scale_input, shape_input, kappa = rep(1, 3), nT, n, 
                         nreps = dim(Pmat)[1], round = TRUE, seed = 1){
  pi <- array(0, c(length(pi_input), dim(pi_input[[1]])))
  scale <- array(0, c(length(scale_input), dim(scale_input[[1]])))
  shape <- array(0, c(length(shape_input), dim(shape_input[[1]])))
  for(j in 1:length(pi_input)) {
      pi[j,,,] <- pi_input[[j]]
      scale[j,,,] <- scale_input[[j]]
      shape[j,,,] <- shape_input[[j]]
  }
  nreps_vec <- 1:nreps
  set.seed(seed)
  data <- array(0, c(nT, n, nreps))

  cdfMat <- Pmat
  cdfMat[,,,,2] <- Pmat[,,,,1]+Pmat[,,,,2]
  cdfMat[,,,,3] <- cdfMat[,,,,2]+Pmat[,,,,3]
  cdfMat[,,,,4] <- cdfMat[,,,,3]+Pmat[,,,,4]
  ## cdfMat[,,,,5] is unused and incorrect

  P_zero_cdf <- t(apply(P_zero, 1, cumsum))
      
  for(i in 1:nT) {
    if(i%%10 == 0) print(i)
    u <- runif(nreps)  
    z <- rep(1, nreps)
    z[u > P_zero_cdf[ , 1]] <- 2
    z[u > P_zero_cdf[ , 2]] <- 3
    z[u > P_zero_cdf[ , 3]] <- 4
    z[u > P_zero_cdf[ , 4]] <- 5

    for(t in 1:n) {
        cb <- cbind(nreps_vec, z)
        u <- runif(nreps)
        cdfMatSub <- cdfMat[, i, t, , ]
        cdf <- cbind(cdfMatSub[,,1][cb],cdfMatSub[,,2][cb],
                     cdfMatSub[,,3][cb],cdfMatSub[,,4][cb])
        z <- rowSums(u > cdf)+1
        zReal <- z - D + 1
        zReal[zReal < 1] <- 1
        idx <- cbind(zReal, nreps_vec)
        pi_tmp <- pi[ , , i, t][idx]
        dry <- rbinom(nreps, 1, pi_tmp)
        scale_tmp <- scale[ , , i, t][idx]
        shape_tmp <- shape[ , , i, t][idx]
        if(dens_type != 2) {
            rainfall <- (scale_tmp/shape_tmp)*
                ((1-runif(nreps)^(1/kappa[zReal]))^(-shape_tmp)-1)
        } else rainfall  <- rgamma(nreps, shape_tmp, scale = scale_tmp)
        if(max(rainfall) > 1000) {
            mu <- shape_tmp*scale_tmp
            sigma <- sqrt(shape_tmp)*scale_tmp
            wh <- which.max(rainfall)
            cat('found high precip: ', i, ' ', t, ' ', rainfall[wh], ' ', mu[wh], ' ', sigma[wh], ' ', zReal[wh], "\n")
        }
        data[i , t, ] <- rainfall*(dry == 0)
        # if(any(data[ , , it] > 1000)) {}
        ## tt <- tt + nonstationary
    }
  }
  ifelse(round, return(round(data, 2)), return(data))
}


get_values <- function(samples, intname, betaname1, betaname2, X1, X2, int_element = 1, beta_element = int_element,
                       doy, link = 'ilogit') {
  nm <- colnames(samples)
  if(is.null(int_element)) {
    int <- samples[ , grep(intname, nm), drop = FALSE]
  } else int <- samples[ , grep(paste0(intname, '\\[', int_element, '\\]'), nm), drop = FALSE]
  if(is.null(beta_element)) {
    betas1 <- samples[ , grep(paste0(betaname1, '\\['), nm)]
    betas2 <- samples[ , grep(paste0(betaname2, '\\['), nm)]
  } else {
      betas1 <- samples[ , grep(paste0(betaname1, '\\[.*', beta_element, '\\]$'), nm)]
      betas2 <- samples[ , grep(paste0(betaname2, '\\[.*', beta_element, '\\]$'), nm)]
  }
  if(ncol(int)) {
      tmp <- int[ , 1]
  } else tmp <- rep(0, nrow(samples))

  tmp <- array(tmp, c(length(tmp), nrow(X2), nrow(X1)))
  if(ncol(betas1) && ncol(betas2)) {
      tmp1 <- X1%*%t(betas1)
      tmp2 <- X2%*%t(betas2)
      ## assume both betas1 and betas2
      for(i in 1:nrow(tmp)) 
          tmp[i, , ] = tmp[i, , ] + c(outer(tmp2[,i], tmp1[,i], '+'))
  }
  fun <- match.fun(link)(tmp)
  return(fun)
}


get_spells <- function(x, cutoff = 0) {
  nT <- length(x)
  wetAmt <- rep(0, nT)
  wetAmtStart <- rep(0, nT)
  wetAmtEnd <- rep(0, nT)
  lengths <- matrix(0, nr = 2, nT)
  lengthsStart <- matrix(0, nr = 2, nT)
  lengthsEnd <- matrix(0, nr = 2, nT)
  idx <- c(0, 0)
  t <- 1
  state = 1 + (x[1] > cutoff)
  idx[state] <- 1
  lengths[state, 1] <- 1
  lengthsStart[state, idx[state]] <- 1
  if(state == 2) {
      wetAmt[1] <- x[1]
      wetAmtStart[1] <- 1
      wetAmtEnd[1] <- 1
  }
  for (t in 2:nT) {
    newState <- 1 + (x[t] > cutoff)
    if(newState == state) {
      lengths[state, idx[state]] <- lengths[state, idx[state]] + 1
      if(state == 2) {
          wetAmt[idx[2]] <- wetAmt[idx[2]] + x[t]
          wetAmtEnd[idx[2]] <- t  # increment
      }
    } else {
      idx[newState] = idx[newState] + 1
      lengths[newState, idx[newState]] <- 1
      lengthsStart[newState, idx[newState]] <- t
      lengthsEnd[state, idx[state]] <- t-1
      if(newState == 2) {
          wetAmt[idx[2]] <- x[t]
          wetAmtStart[idx[2]] <- t
          wetAmtEnd[idx[2]] <- t
      } 
    }
    state <- newState
  }
  ## assign end date as last date for last event
  if(state == 2) 
      wetAmtEnd[idx[2]] <- t
  lengthsEnd[state, idx[state]] <- t
      
  # cut off empty elements
  wetAmt <- wetAmt[wetAmt > cutoff]  
  dryLen <- lengths[1, ] 
  dryLen <- dryLen[dryLen > 0]
  wetLen <- lengths[2, ]
  wetLen <- wetLen[wetLen > 0]
  wetAmtStart <- wetAmtStart[wetAmtStart != 0]
  wetAmtEnd <- wetAmtEnd[wetAmtEnd != 0]
  dryLenStart <- lengthsStart[1, lengthsStart[1, ] != 0]
  dryLenEnd  <- lengthsEnd[1, lengthsEnd[1, ] != 0]
  wetLenStart <- lengthsStart[2, lengthsStart[2, ] != 0]
  wetLenEnd <- lengthsEnd[2, lengthsEnd[2, ] != 0]
  return(list(wetAmt = cbind(wetAmt, wetAmtStart, wetAmtEnd),
              dryLen = cbind(dryLen, dryLenStart, dryLenEnd),
              wetLen = cbind(wetLen, wetLenStart, wetLenEnd),
              wetAvg = mean(x[x > cutoff]),
              wetMed = median(x[x > cutoff])))
}


get_spells_zhang <- function(x, cutoff = 0) {
  if(all(is.na(x))) {
        return(list(wetAmt = NULL,
              dryLen = NULL,
              wetLen = NULL))
  }
  nT <- length(x)
  wetAmt <- rep(0, nT)
  wetAmtStart <- rep(0, nT)
  wetAmtEnd <- rep(0, nT)
  lengths <- matrix(0, nr = 2, nT)
  lengthsStart <- matrix(0, nr = 2, nT)
  lengthsEnd <- matrix(0, nr = 2, nT)
  idx <- c(0, 0)
  t <- 1
  state = 1 + (x[1] > cutoff)
  idx[state] <- 1
  lengths[state, 1] <- 1
  lengthsStart[state, idx[state]] <- 1
  if(state == 2) {
      wetAmt[1] <- x[1]
      wetAmtStart[1] <- 1
      wetAmtEnd[1] <- 1
  }
  for (t in 2:nT) {
    newState <- 1 + (x[t] > cutoff)
    if(newState == state) {
      lengths[state, idx[state]] <- lengths[state, idx[state]] + 1
      if(state == 2) {
          wetAmt[idx[2]] <- wetAmt[idx[2]] + x[t]
          wetAmtEnd[idx[2]] <- t  # increment
      }
    } else {
      idx[newState] = idx[newState] + 1
      lengths[newState, idx[newState]] <- 1
      lengthsStart[newState, idx[newState]] <- t
      lengthsEnd[state, idx[state]] <- t-1
      if(newState == 2) {
          wetAmt[idx[2]] <- x[t]
          wetAmtStart[idx[2]] <- t
          wetAmtEnd[idx[2]] <- t
      } 
    }
    state <- newState
  }
  ## assign end date as last date for last event
  if(state == 2) 
      wetAmtEnd[idx[2]] <- t
  lengthsEnd[state, idx[state]] <- t
      
  # cut off empty elements
  wetAmt <- wetAmt[wetAmt > cutoff]  
  dryLen <- lengths[1, ] 
  dryLen <- dryLen[dryLen > 0]
  wetLen <- lengths[2, ]
  wetLen <- wetLen[wetLen > 0]
  wetAmtStart <- wetAmtStart[wetAmtStart != 0]
  wetAmtEnd <- wetAmtEnd[wetAmtEnd != 0]
  dryLenStart <- lengthsStart[1, lengthsStart[1, ] != 0]
  dryLenEnd  <- lengthsEnd[1, lengthsEnd[1, ] != 0]
  wetLenStart <- lengthsStart[2, lengthsStart[2, ] != 0]
  wetLenEnd <- lengthsEnd[2, lengthsEnd[2, ] != 0]

  
  ## cut off dry spells extending into next year
  ## if(next_year == 0) {
  ##     n  <- length(dryLenEnd)
  ##     if(dryLenEnd[n] == 365) {
  ##         dryLen <- dryLen[-n]
  ##         dryLenStart <- dryLenStart[-n]
  ##         dryLenEnd <- dryLenEnd[-n]
  ##     }
  return(list(wetAmt = cbind(wetAmt, wetAmtStart, wetAmtEnd),
              dryLen = cbind(dryLen, dryLenStart, dryLenEnd),
              wetLen = cbind(wetLen, wetLenStart, wetLenEnd)))
}

analyze_trend <- function(data, full = FALSE) {
    tmp <- modifiedmk::mkttest(data)
    if(full) return(tmp) else return(tmp[2])
}


calc_trans <- function(data, cutoffs = dryCut/cmPerInch, vectorOut = TRUE) {
    n <- length(data)
    cats <- data
    cats[data < cutoffs[1]] <- 0
    for(j in seq_along(cutoffs)) 
        cats[data >= cutoffs[j]] <- j
    cats <- factor(cats, levels = 0:length(cutoffs))
    result <- table(cats[1:(n-1)], cats[2:n]) # , margin = 1)
    result <- result/rowSums(result)  ## faster than prop.table
    ## result <- prop.table(table(cats[1:(n-1)], cats[2:n]), margin = 1)
    if(vectorOut) return(c(result)) else return(result)
}
