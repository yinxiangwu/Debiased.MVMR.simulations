#!/bin/bash -l
#SBATCH --output=./io/%j.out
Rscript sim_3exp_overlap_10102023.R

# run it like:
# sbatch --array=1-1000 Condor.sh