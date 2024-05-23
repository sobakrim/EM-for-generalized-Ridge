# EM Algorithm for Generalized Ridge Regression with Spatial Covariates

This repository contains code and documentation for implementing the EM algorithm for Generalized Ridge Regression with spatial covariates. The method is demonstrated using a simulation example, and the estimated beta values are compared across different covariance models (Matern, CAR, and Ridge).

## Table of Contents
- [Introduction](#introduction)
- [Method](#method)
- [Simulation Example](#simulation-example)
- [Results](#results)
- [Usage](#usage)
- [License](#license)

## Introduction

This project presents the EM algorithm for Generalized Ridge Regression with spatial covariates. The method is demonstrated using a simulation example where the estimated beta values are compared across different covariance models (Matern, CAR, and Ridge).

## Method

The method involves:
- Loading the required libraries.
- Defining necessary functions for covariance computation and simulation.
- Implementing the EM algorithm to estimate the regression parameters.

## Simulation Example

A simulation example is provided to demonstrate the method. The example includes:
- Simulating spatial data using a Matern covariance function.
- Fitting the EM algorithm for different covariance models.
- Comparing the estimated beta values against the true beta values.

## Results

The results of the simulation example include visualizations of the true and estimated beta spatial fields for each covariance model.

## Usage

To use the code in this repository:
1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/em_spatial_covariates.git
    ```
2. Open the `EM_Algorithm_Spatial_Covariates.Rmd` file in RStudio.
3. Run the R Markdown file to generate the HTML report with the simulation example and results.

### Required Libraries

The following R packages are required:
- Matrix
- geoR
- MASS
- corpcor
- ggplot2
- reshape2

You can install the required packages using:
```r
install.packages(c("Matrix", "geoR", "MASS", "corpcor", "ggplot2", "reshape2"))
