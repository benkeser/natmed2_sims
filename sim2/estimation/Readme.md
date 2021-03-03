# Code for simulation 2: estimation

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Authors:** [David
Benkeser](https://davidbphd.com), Iván Díaz, Jialu Ran

-----

## Description

This folder contains all of the code needed to execute the simulation for producing estimation and confidence interval, plots for risk estimators, indirect effects, prop. mediated on a HPC cluster.
- `code` contains all scripts needed to execute the simulation and format the results. 
- `produce` contains all the scripts needed to produce the tables and plots in the paper.

The simulation uses the `future` package to use multicore parallelization when executing the simulation. If parallelization is not available, the simulation code will execute sequentially. On Emory University Biostatistics and Bioinformatics HPC cluster, the simulation executed in about two minutes for glm main, five minutes for glm interaction, and eight hours for Superlearner.

-----

## Dependencies

To execute this code, the following `R` packages are needed:
- [`natmed2`](https://github.com/benkeser/natmed2)
- `future`
- `future_apply`
- `SuperLearner`
- `ggplot2`
- `ggthemes`
- `grid`
- `xtable`
- `extrafont`

In `code` folder, a path of the location on the cluster is needed to run the run_sim.R and submit.sh file. The helper_fn.R should be put in the same folder as run_sim.R, submit.sh on the cluster. A output folder needed to be created to store the data from the simulation.

