# Debiased.MVMR.simulations
This repository contains the R code to conduct the simulation studies covered in the paper Debiased Multivariable Mendelian Randomization [arXiv](https://arxiv.org/abs/2402.00307).

The functions to perform MVMR-IVW, MVMR-dIVW, MVMR-adIVW are available in the R package mr.divw,  which is posted at https://github.com/tye27/mr.divw.

There are two main folders:
1. R code for simulations
2. R code to summarize simulation results

## Steps to generate results for Table 1, and Table S1-7: 

run simulation_study1.R first (repeat 100 times if parallel computation is used), and then run summarize_results_simulation_study1.R

## Steps to generate results for Table S8-11: 
run simulation_study1_bhp.R first (repeat 100 times if parallel computation is used), and then run summarize_results_simulation_study1_bhp.R

## Steps to generate results for Table 2: 
run simulation_study2.R first (repeat 100 times if parallel computation is used), and then run summarize_results_simulation_study2.R

## Steps to generate results for Table 3: 
run simulation_study3.R first (repeat 100 times if parallel computation is used), and then run summarize_results_simulation_study3.R

## Steps to generate results for Table S13-16: 
run simulation_study1_overlap.R first (repeat 100 times if parallel computation is used), and then run summarize_results_simulation_study1_overlap.R

Other programs in the R code for simulations folder:

1. bayes1cluster.sh: bash file for parallel computing
2. data_gen_individual.R R program for generating individual-level data and corresponding summary statistics for simulation study 3
3. F_stats_calculator.R: calculate conditional F-statistics, adapted from MVMR package.

Users are encouraged to contact the authors with any questions.
