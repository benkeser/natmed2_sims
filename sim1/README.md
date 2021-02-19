# Code for simulation 1

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Authors:** [David
Benkeser](https://davidbphd.com), Iván Díaz, Jialu Ran

-----

## Description

This folder contains all of the code needed to execute the first simulation.
- `code` contains all scripts needed to execute the simulation and format the results
- `output` is an empty directory where simulation output will be saved

The simulation uses the `future` package to use multicore parallelization when executing the simulation. If parallelization is not available, the simulation code will execute sequentially. On a 14 core machine with hyperthreading, the simulation executed in about ten minutes.

-----

## Makefile

To reproduce the simulation you can simply type `make`. This will produce three files in the `output` directory:
- `latex_table.tex` = raw latex output used to create the table
- `table.pdf` = a compiled pdf displaying the table of results
- `sim_result.rds` = a `data.frame` containing raw simulation results

Type `make help` to see the available recipes for `make`.

-----

## Dependencies

To execute this code, the following `R` packages are needed:
- [`natmed2`](https://github.com/benkeser/natmed2)
- `here`
- `future`
- `future_apply`
- `xtable`

You will also need `pdflatex` installed in order to compile the standalone `pdf` table.