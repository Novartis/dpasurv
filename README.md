dpasurv
====================

R-package for Dynamic Path Analysis of survival data.

## Overview

The *dpasurv* package is designed to perform dynamic path analysis via calls to the following principal functions:

- `dpa()` fits the dynamic path models corresponding to a given dynamic path diagram; returns an object of class dpa
- `effect()` estimates the cumulative direct and indirect effects from the fitted dpa-object; returns an object of class effect

The effect-objects can be summed together to obtain the total effect (i.e. direct effect + indirect effect) as well as plotted with bootstrap confidence intervals.
