# Code for simulation 2 for getting true values

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Authors:** [David
Benkeser](https://davidbphd.com), Iván Díaz, Jialu Ran

-----

## Description

This folder contains all of the code needed to execute the simulation for producing table 2 in the paper.
- `code` contains all scripts needed to execute the simulation and format the results
- `output` is an empty directory where simulation output will be saved

The simulation uses the `future` package to use multicore parallelization when executing the simulation. If parallelization is not available, the simulation code will execute sequentially. On Emory University Biostatistics and Bioinformatics HPC cluster, the simulation executed in about three minutes.

-----

## Dependencies

To execute this code, the following `R` packages are needed:
- [`natmed2`](https://github.com/benkeser/natmed2)
- `here`
- `future`
- `future_apply`
- `xtable`

