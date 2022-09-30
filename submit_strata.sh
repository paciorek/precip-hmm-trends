#!/bin/bash
Rscript run_strata.R ${1} ${2} ${3} ${4} ${5} ${6} ${7} &> logs/run-${1}-season${2}-${3}-model${4}-holdout${5}-its${6}-rep${7}.Rout

## for (( i=1; i<=5; i++)); do sbatch -p high submit_strata.sh AZ_Prescott 4 trend 2 0 20000 $i; done
