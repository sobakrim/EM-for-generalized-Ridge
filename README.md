<h1 align="center">EM Algorithm for Generalized Ridge Regression with Spatial Covariates</h1>

<p align="center">
  <em>R implementation of a spatially-extended ridge regression using the EM algorithm</em>
</p>

---

<h2> Overview</h2>

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
The full approach is described in:
</p>

<blockquote>
  <p>
    Obakrim, S. <em>et&nbsp;al.</em> “EM algorithm for generalized Ridge regression with spatial covariates,”
    <em>Environmetrics</em>, 2024.
    <a href="https://doi.org/10.1002/env.2871">doi:10.1017/eds.2022.35</a>
  </p>
</blockquote>
