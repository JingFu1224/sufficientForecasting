
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sufficientForecasting

<!-- badges: start -->
<!-- badges: end -->

## Overview

The goal of sufficientForecasting is to forecast a single time series
when there is a large number of predictors and a possible nonlinear
effect.

## Installation

You can install the development version of sufficientForecasting like
so:

``` r
# The easiest way to install sufficientForecasting
install.packages("sufficientForecasting")
```

## Usage

The following example uses SF.CI to solve a problem: forecast a single
time series, and its upper bound and lower bound

``` r
library(sufficientForecasting)
## basic example code
SF.CI(y=dataExample$y,X=dataExample$X,newX=dataExample$newX,type="LLM",alpha = 0.05)
#>     yhat ci_lower ci_upper 
#>  -0.3568  -2.4740   1.6076
```
