# A Hierarchical Three-Population Mortality Model: An Extension of Two-Population Frameworks for Assessing Longevity Basis Risk

This repository implements a hierarchical three-population mortality modelling framework using Human Mortality Database (HMD) data. It provides a fully reproducible pipeline for estimating mortality dynamics, generating forecasts, and evaluating index-based longevity hedging strategies under model uncertainty.

The workflow is organised into five modular R Markdown files, enabling users to adapt and extend the code to evaluate alternative populations, model specifications, and hedging strategies.

---

## Overview

The project constructs and fits a hierarchical mortality model across three populations:

- **R1**: First reference population  
- **R2**: Second reference population  
- **B**: Book population (representing insurer or pension fund exposure)

The hierarchical structure enables information sharing across populations, improving coherence and stability in parameter estimation.

---

## Key Features

- Hierarchical three-population mortality modelling framework  
- Cairns–Blake–Dowd (CBD) models (M5, M6, M7)  
- Time-series forecasting of mortality parameters  
- Residual bootstrap for uncertainty quantification  
- Index-based longevity hedging framework  
- Simulation-based longevity risk evaluation (5,000 replications)

---

## Repository Structure

The project is implemented as a five-stage modular pipeline using R Markdown:

### Core workflow

- `1 - HMD Dataset.Rmd`  
  Data extraction and preprocessing from the Human Mortality Database (HMD), including population selection and dataset construction.

- `2 - Model Fitting and Plots.Rmd`  
  Estimation of hierarchical three-population mortality models (CBD family: M5, M6, M7) and visualisation of fitted results.

- `3 - Time Series Model and Forecast.Rmd`  
  Time-series modelling of mortality parameters and generation of forecasted mortality rates.

- `4 - Bootstrapping and Simulation.Rmd`  
  Residual bootstrap procedure for uncertainty quantification and simulation of mortality trajectories.

- `5 - Hedging.Rmd`  
  Construction and evaluation of index-based longevity hedging strategies, including risk reduction analysis.

---

### Data

- `HMD_Countries.csv`  
  Mapping file for HMD country codes and labels used for dataset construction.

---

### Outputs

Each R Markdown file generates a corresponding `.md` file for GitHub rendering (e.g. `1---HMD-Dataset.md`), documenting methodology and results.

---

## How to Run

1. Clone the repository  
2. Open the R project in RStudio  
3. Run the R Markdown files sequentially:

   1 → HMD Dataset  
   2 → Model Fitting and Plots  
   3 → Time Series Model and Forecast  
   4 → Bootstrapping and Simulation  
   5 → Hedging  

Each stage builds on outputs from the previous step.

---

## Reproducibility

All results are reproducible given the same HMD data and R package versions. Simulation-based procedures use fixed seeds where applicable.

---

## Required Packages

```r
install.packages(c(
  "StMoMo",
  "demography",
  "stats",
  "forecast",
  "vars",
  "MASS"
))
