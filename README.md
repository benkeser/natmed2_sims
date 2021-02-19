# Supplementary code for Inference for natural mediation effects under case-cohort sampling with applications in identifying COVID-19 vaccine correlates of protection

[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Authors:** [David
Benkeser](https://davidbphd.com), Iván Díaz, Jialu Ran

-----

## Description

The `natmed2_sims` repository contains the code needed to reproduce the simulation results of the manuscript, "Inference for natural mediation effects under case-cohort sampling with applications in identifying COVID-19 vaccine correlates of protection."

-----

## The `natmed2` package

The simulations rely on the [`natmed2`](https://github.com/benkeser/natmed2) R package, which can be downloaded and installed from GitHub as follows.

```r
remotes::install_github("benkeser/natmed2")
```

-----

## Simulation 1

Code for the first simulation is included in the `sim1` directory, which includes its own `README` file detailing the contents. This simulation study verified the theoretical properties of the proposed estimators in large samples.

-----

## Simulation 2

Code for the second simulation is included in the `sim2` directory, which includes its own `README` file detailing the contents. This simulation study studied the properties of the estimators in the context of a large randomized trial of a COVID-19 vaccine.

-----

## Issues

If you encounter any bugs, please [file an issue](https://github.com/benkeser/natmed2_sim2/issues).

-----

## License

© 2021- [David Benkeser](https://davidbphd.com)

The contents of this repository are distributed under the MIT license.
See below for details:

    The MIT License (MIT)
    
    Copyright (c) 2021- David C. Benkeser
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
