# A Hierarchical Three-Population Mortality Model: An Extension of Two-Population Frameworks for Assessing Longevity Basis Risk

This repository implements a hierarchical three-population mortality modelling framework using Human Mortality Database (HMD) data. The framework is designed for reproducible analysis of mortality dynamics, forecasting, and index-based longevity hedging. The five R Markdown files provide a modular workflow that allows users to flexibly modify populations, model specifications, and hedging assumptions to suit their own applications.

---

## Overview

The project constructs and fits a **hierarchical mortality model** across three populations:

- **R1**: First reference population  
- **R2**: Second reference population  
- **B**: Book population (insurer/pension fund exposure)

The hierarchical structure allows information sharing across populations, improving stability and coherence in parameter estimation.

---

## Mortality Models

We consider the Cairns–Blake–Dowd (CBD) family of mortality models.

### M5
```math
logit(q_{x,t}) = \kappa_{t,1} + (x - \bar{x})\kappa_{t,2} = \eta_{x,t}
```

### M6
```math
logit(q_{x,t}) = \kappa_{t,1} + (x - \bar{x})\kappa_{t,2} + \gamma_{t-x} = \eta_{x,t}
```

### M7
```math
logit(q_{x,t}) = \kappa_{t,1} + (x - \bar{x})\kappa_{t,2} + \left((x - \bar{x})^2 - \sigma^2\right)\kappa_{t,3} + \gamma_{t-x} = \eta_{x,t}
```

---

## Hierarchical Structure

The model is estimated sequentially across populations:

### R1
```math
logit(q^{R1}_{x,t}) = \eta^{R1}_{x,t}
```

### R2
```math
logit(q^{R2}_{x,t}) - logit(q^{R1}_{x,t}) = \eta^{R2}_{x,t}
```

### B
```math
logit(q^{B}_{x,t}) - logit(q^{R2}_{x,t}) = \eta^{B}_{x,t}
```

This structure ensures later populations are informed by earlier reference populations, improving estimation stability.

---

## Time-Series Modelling

Estimated mortality parameters are forecast using the following time-series specifications.

### Period effects

**R1 (MRWD model):**
```math
\kappa^{R1}_{t} = d + \kappa^{R1}_{t-1} + \epsilon_{t}
```

**R2 and Book Population (VAR model):**
```math
\boldsymbol{\kappa}^{i}_{t} = \mathbf{c} + \mathbf{A}\boldsymbol{\kappa}^{i}_{t-1} + \boldsymbol{\epsilon}^{i}_{t}, 
\quad i \in \{R2, B\}
```

---

### Cohort effects

**R1 (ARIMA(1,1,0)):**
```math
\gamma^{R1}_{t-x} - \gamma^{R1}_{t-x-1}
=
\phi_{0} + \phi_{1}\left(\gamma^{R1}_{t-x-1} - \gamma^{R1}_{t-x-2}\right) + \xi_{t}
```

**R2 and Book Population (AR(1)):**
```math
\gamma^{i}_{t-x} = \phi_{0} + \phi_{1}\gamma^{i}_{t-x-1} + \xi_{t}, 
\quad i \in \{R2, B\}
```
---

## Bootstrap Procedure

Uncertainty is quantified using a **residual bootstrap**:

- Deviance residuals are resampled with replacement
- Samples are converted into death counts (Binomial framework)
- The model is refitted for each replication (B = 5000)
- Parameters are projected forward using time-series models

This procedure takes approximately **55 minutes** to run.

---

## Longevity Hedging

We evaluate **index-based longevity swaps** calibrated to the two reference populations (R1 and R2) to hedge the book population.

### Longevity Risk Reduction
```math
\text{Longevity Risk Reduction} = \left(1 - \frac{\text{Hedged Risk}}{\text{Unhedged Risk}}\right)\times 100\%
```

### Portfolio Value (Simulation i)
```math
PV^{(i)} =
\sum_{t=1}^{25} \frac{l^{B(i)}_{65+t,\,2020+t}}{(1+r)^t}
- w_1 \sum_{t=1}^{25} \frac{\left({}_tp^{R1(i)}_{65,2020} - {}_tp^{R1,\text{forward}}_{65,2020}\right)}{(1+r)^t}
- w_2 \sum_{t=1}^{25} \frac{\left({}_tp^{R2(i)}_{65,2020} - {}_tp^{R2,\text{forward}}_{65,2020}\right)}{(1+r)^t},
\quad i = 1,\ldots,5000.
```

The weights w1 and w2 represent the hedge weights of the two index-based longevity swaps, which are calibrated to minimise the 99.5% Value-at-Risk (VaR), minus the mean, of the present value of the hedged portfolio under the distribution generated from the 5,000 simulations.
---

## Data

- Human Mortality Database (HMD)
- Population-specific mortality rates and exposures
- Custom CSV mapping for country codes and labels

---

## Key Features

- Hierarchical multi-population mortality modelling
- CBD model implementation (M5–M7)
- Time-series forecasting of mortality parameters
- Residual bootstrap uncertainty quantification
- Longevity hedging via index-based swaps
- Simulation-based risk assessment

---

## Required Packages

The following R packages are required to run this project:

```r
install.packages(c(
  "StMoMo",
  "demography",
  "stats",
  "forecast",
  "vars",
  "MASS"
))
```

### Overview
- **StMoMo**: Stochastic mortality modelling (CBD framework implementation)  
- **demography**: Access and manipulation of mortality data (HMD integration)  
- **stats**: Base statistical functions for modelling and estimation  
- **forecast**: Time-series forecasting of model parameters  
- **vars**: Vector autoregressive (VAR) modelling for multi-population dynamics  
- **MASS**: Additional statistical tools used in estimation and simulation  

---

## Output

- Fitted mortality models
- Forecasted mortality rates
- Bootstrap distributions of parameters
- Longevity Risk reduction estimates

---

## Notes

- Model selection (M5/M6/M7) is based on BIC across populations.
- Hierarchical modelling may not always yield globally optimal BIC for all populations, depending on demographic structure.
