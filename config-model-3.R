## have time-varying splines have fitted variance components

## 1: GPD
## 2: Gamma
## 3: Naveau (generalized GPD)

dens_type <- 2
rounded <- TRUE
allow_missing <- TRUE
use_inches <- TRUE
use_splines <- TRUE
use_halfflat <- FALSE

beta_sd <- 0.05
max_sigma <- 3

use_new_blockRW <- TRUE
block_ints <- FALSE

adaptFactorExponent <- 0.25

constrain_ints <- TRUE

constrain_full_means <- TRUE

D <- 3
W <- 2

## spline specifications

pD_SEAS <- TRUE
pDW_SEAS <- TRUE
pW_SEAS <- rep(TRUE, W); any_pW_SEAS <- any(pW_SEAS)
pWW_SEAS <- rep(TRUE, W); any_pWW_SEAS <- any(pWW_SEAS)
pWD_SEAS <- rep(TRUE, D-1); any_pWD_SEAS <- any(pWD_SEAS)
pi_SEAS <- rep(FALSE, W+1); any_pi_SEAS <- any(pi_SEAS)
sigma_SEAS <- c(FALSE, TRUE, TRUE); any_sigma_SEAS <- any(sigma_SEAS)
xi_SEAS <- c(FALSE, TRUE, TRUE); any_xi_SEAS <- any(xi_SEAS)

if(trend == 'trend') {
    pD_TIME <- TRUE
    pDW_TIME <- TRUE
    pW_TIME <- rep(TRUE, W); any_pW_TIME <- any(pW_TIME)
    pWW_TIME <- rep(TRUE, W); any_pWW_TIME <- any(pWW_TIME)
    pWD_TIME <- rep(TRUE, D-1); any_pWD_TIME <- any(pWD_TIME)
    pi_TIME <- rep(FALSE, W+1); any_pi_TIME <- any(pi_TIME)
    sigma_TIME <- c(FALSE, TRUE, TRUE); any_sigma_TIME <- any(sigma_TIME)
    xi_TIME <- c(FALSE, TRUE, TRUE); any_xi_TIME <- any(xi_TIME)
} else {
    pD_TIME <- FALSE
    pDW_TIME <- FALSE
    pW_TIME <- rep(FALSE, W); any_pW_TIME <- any(pW_TIME)
    pWW_TIME <- rep(FALSE, W); any_pWW_TIME <- any(pWW_TIME)
    pWD_TIME <- rep(FALSE, D-1); any_pWD_TIME <- any(pWD_TIME)
    pi_TIME <- rep(FALSE, W+1); any_pi_TIME <- any(pi_TIME)
    sigma_TIME <- c(FALSE, FALSE, FALSE); any_sigma_TIME <- any(sigma_TIME)
    xi_TIME <- c(FALSE, FALSE, FALSE); any_xi_TIME <- any(xi_TIME)
}

K1 <- 20
K2 <- 20

nburnin = 0
thin <- 10
thin2 <- 10

burnin_waic_frac <- 0.25
