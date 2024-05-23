# EM Algorithm for Generalized Ridge Regression with Spatial Covariates

This repository contains code and documentation for implementing the EM algorithm for Generalized Ridge Regression with spatial covariates. The method is demonstrated using a simulation example, and the estimated beta values are compared across different covariance models (Matern, CAR, and Ridge).

## Table of Contents
- [Introduction](#introduction)
- [Method](#method)
- [Simulation Example](#simulation-example)
- [Results](#results)

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

### Required Libraries

The following R packages are required:
- Matrix
- geoR
- MASS
- corpcor
- ggplot2
- reshape2
