<h1 align="center">EM Algorithm for Generalized Ridge Regression with Spatial Covariates</h1>

<p align="center">
  <em>R implementation of a spatially-extended ridge regression using the EM algorithm</em>
</p>

---

<h2>📘 Overview</h2>

This repository provides an R implementation of the <strong>Expectation-Maximization (EM) algorithm</strong> for <strong>Generalized Ridge Regression</strong> with spatial covariates.  
The method is evaluated via a simulation study comparing regression coefficient estimates under different spatial priors:

<ul>
  <li><strong>Matérn covariance</strong></li>
  <li><strong>Conditional Autoregressive (CAR) model</strong></li>
  <li><strong>Standard (non-spatial) Ridge</strong></li>
</ul>

---

<h2>📊 Simulation Study</h2>

The provided simulation illustrates:
<ul>
  <li>Generating spatially correlated covariates and responses using a Matérn process</li>
  <li>Running the EM algorithm for each covariance structure</li>
  <li>Comparing the estimated β coefficients to the true underlying values</li>
</ul>

---

<h2>🧰 Required R Packages</h2>

Make sure to install the following R packages before running the code:

```r
install.packages(c(
  "Matrix", "geoR", "MASS", "corpcor", "ggplot2", "reshape2"
))

<h2>🚀 Getting Started</h2>

Clone this repository and run the main script:

source("simulation_EM_spatial_ridge.R")

The script outputs both plots and numerical results comparing estimation performance across models.
<h2>📂 Repository Structure</h2> <ul> <li><code>simulation_EM_spatial_ridge.R</code> — Main simulation and estimation script</li> <li><code>README.md</code> — Project documentation</li> </ul> 
