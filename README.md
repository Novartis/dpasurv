dpasurv
====================

R-package for Dynamic Path Analysis of survival data.

## Overview

The *dpasurv* package is designed to perform dynamic path analysis via calls to the following principal functions:

- `dpa()` fits the dynamic path models corresponding to a given dynamic path diagram; returns an object of class dpa
- `effect()` estimates the cumulative direct and indirect effects from the fitted dpa-object; returns an object of class effect

The effect-objects can be summed together to obtain the total effect (i.e. direct effect + indirect effect) as well as plotted with bootstrap confidence intervals.

## Installation

``` r
# The easiest way to install dpasurv is to run (assuming devtools package is already installed):
devtools::install_github("Novartis/dpasurv")

# If devtools package is not installed, then prior to above you may run:
install.packages("devtools")
```
## Usage
``` r
library(dplyr)

# Perform dynamic path analysis
analysis <- dpa(survival::Surv(start, stop, event) ~ M + x, list(M ~ x), id = "subject", data = simdata, boot.n = 100)

# Extract direct, indirect and total effect
direct <- effect(x ~ outcome, analysis, alpha=0.05)
indirect <- effect(x ~ M ~ outcome, analysis, alpha=0.05)
total <- sum(direct, indirect)

# Plot the results
par(mfrow=c(1,3))
plot(direct); abline(h=0, lty=2, col=2)
plot(indirect); abline(h=0, lty=2, col=2)
plot(total); abline(h=0, lty=2, col=2)

# Plot the results with ggplot2 graphics:
ggplot.effect(list(direct, indirect, total))
```
