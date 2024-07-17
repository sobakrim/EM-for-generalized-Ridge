# EM Algorithm for Generalized Ridge Regression with Spatial Covariates

This repository contains code and documentation for implementing the EM algorithm for Generalized Ridge Regression with spatial covariates. The method is demonstrated using a simulation example, and the estimated beta values are compared across different covariance models (Matern, CAR, and Ridge).

## Introduction

This project presents the EM algorithm for Generalized Ridge Regression with spatial covariates. The method is demonstrated using a simulation example where the estimated beta values are compared across different covariance models (Matern, CAR, and Ridge).

## Simulation Example

A simulation example is provided to demonstrate the method. The example includes:
- Simulating spatial data using a Matern covariance function.
- Fitting the EM algorithm for different covariance models.
- Comparing the estimated beta values against the true beta values.

### Required Libraries

The following R packages are required:
- Matrix
- geoR
- MASS
- corpcor
- ggplot2
- reshape2
