## Code run on SCF machines, under R 4.1.2 in spring 2022

## run ../data/find_long_records_rnoaa.R
## clean up the file and add code for Winslow

## on SCF cluster, run `submit_strata.sh`
## standard runs use config-model-5.R, config-holdout-0.R (gamma, no holdout, with time)
## comparison runs use:
## config-model-4.R, config-holdout-1.R, 'trend' (GPD, holdout, with time)
## config-model-5.R, config-holdout-1.R, 'trend' (gamma, holdout, with time)
## config-model-5.R, config-holdout-1.R, 'notrend' (gamma, holdout, without time)
## this invokes run_strata.sh and then fit.R

## generate posterior predictive simulations for all 10k saved samples and all 122 years
## this requires a machine with > 100 GB memory, and possibly 256 GB for Quincy
## still in units of inches, which is the scaled used for fitting, given rounding to .01 inches

data_fn='AZ_Litchfield'
season=2
trend='trend'
config_model=5
config_holdout=0
niter=20000

## Loop over all 8 sets
./simulate.R data_fn season trend config_model config_holdout niter

## calculate spells; loop over all 8 sets
./calculate_spells.R data_fn season trend config_model config_holdout niter

## check mixing for all 8 dataset; generate Table 1, Table 3 (Appendix)
./assess_mixing.R

## model selection for 4 datasets: calculate test LL and WAIC (and generate imputations via FFBS):
## Loop over 4 sets
./calculate_waic_ll.R data_fn season trend config_model config_holdout niter

## compare models; generate Table 2
./compare_models.R

## model diagnostics: generate various figures (Figs 1, 2, 8-19)
./assess_models.R

## case study trend analysis: generate Figs. 3-7
./assess_trends.R
