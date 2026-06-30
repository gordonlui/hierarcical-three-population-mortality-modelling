2 - Model Fitting and Plots
================
Yiu Chung Lui (Gordon)

## Prerequisite - Load Previous Environment

``` r
load("environment.RData")
```

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

## 2.1 - Model Fitting and BIC Evaluation

This section allows users to fit a hierarchical three-population
mortality model using population datasets derived from the HMD data,
which involves the first reference population (R1), second reference
population (R2), and the book population (B). Currently, this framework
supports only the Cairns–Blake–Dowd (CBD) family of models, specifically
the M5, M6, and M7 models. The Cairns–Blake–Dowd (CBD) family of models
considered in this framework are defined as follows:

M5:

``` math
logit(q_{x,t}) = \kappa_{t,1} + (x - \bar{x})\kappa_{t,2} = \eta_{x,t}
```

M6:

``` math
logit(q_{x,t}) = \kappa_{t,1} + (x - \bar{x})\kappa_{t,2} + \gamma_{t-x} = \eta_{x,t}
```

M7:

``` math
logit(q_{x,t}) = \kappa_{t,1} + (x - \bar{x})\kappa_{t,2} + \left((x - \bar{x})^2 - \sigma^2\right)\kappa_{t,3} + \gamma_{t-x} = \eta_{x,t}
```

The hierarchical approach proceeds by sequentially fitting models across
the three populations, where parameter estimates from earlier
(reference) populations are used to inform and constrain the estimation
of subsequent populations, thereby improving coherence and stability
across the joint modelling framework.

First Reference Population (R1):

``` math
logit(q^{R1}_{x,t}) = \eta^{R1}_{x,t}
```

Second Reference Population (R2):

``` math
logit(q^{R2}_{x,t}) - logit(q^{R1}_{x,t}) = \eta^{R2}_{x,t}
```

Book Population (B):

``` math
logit(q^{B}_{x,t}) - logit(q^{R2}_{x,t}) = \eta^{B}_{x,t}
```

### Load packages

``` r
library(StMoMo)
```

    ## Loading required package: gnm

    ## Loading required package: forecast

``` r
library(demography)
```

    ## Registered S3 methods overwritten by 'demography':
    ##   method      from 
    ##   print.lca   e1071
    ##   summary.lca e1071

    ## This is demography 2.0

``` r
library(stats)
```

### Define country labels

``` r
labels <- names(stmomo_datasets)
```

### Specify Model Equations

``` r
# Model equations
M7 <- m7(link = "logit")
M6 <- m6(link = "logit")
M5 <- cbd(link = "logit")
```

### Inverse Logit Function (“inv_logit”)

``` r
# Inverse logit function
inv_logit <- function(x) {
  1/(1 + exp(-x))
}
```

### Specify R1, R2, and B Population (UK, ENW, and sim_ENW example)

``` r
# Specify model populations
model_populations <- c("UK", "ENW", "sim_ENW")

# Specify dataset into fit 
R1_data <- stmomo_datasets[[model_populations[1]]]
R2_data <- stmomo_datasets[[model_populations[2]]]
B_data  <- stmomo_datasets[[model_populations[3]]]

# R1, R2 and B qxt
R1_qxt <- R1_data$Dxt / R1_data$Ext
R2_qxt <- R2_data$Dxt / R2_data$Ext
B_qxt  <- B_data$Dxt / B_data$Ext
```

### BIC table

``` r
# Create model list
model_list <- list(
  M7 = M7,
  M6 = M6,
  M5 = M5
)

# Create matrix to store BIC
models_BIC <- matrix(NA, nrow = 27, ncol = 7)
colnames(models_BIC) <- c("R1 Model", "R2 Model", "B Model", "R1 BIC", "R2 BIC", "B BIC", "Joint BIC")

# Loop for all model specifications
row <- 1

for (i in 1:length(model_list)) {
  # Specify R1 model
  R1 <- names(model_list[i])
  R1_model <- model_list[[i]]

  for (j in 1:length(model_list)) {
    # Specify R2 model
    R2 <- names(model_list[j])
    R2_model <- model_list[[j]]
    
    for (k in 1:length(model_list)) {
      # Specify B_model
      B_3 <- names(model_list[k])
      B_3_model <- model_list[[k]]
      
      # R1 Fit
      hide_output_dummy <- capture.output(R1_fit <- fit(R1_model, data = R1_data))
      R1_fitted <- fitted(R1_fit)
      # R2 Fit
      hide_output_dummy <- capture.output(R2_fit <- fit(R2_model, data = R2_data, oxt = R1_fitted))
      R2_fitted <- fitted(R2_fit)
      # B_3 Fit
      hide_output_dummy <- capture.output(B_3_fit <- fit(B_3_model, data = B_data, oxt = R2_fitted))
      B_3_fitted <- fitted(B_3_fit)
      
      # Store results
      models_BIC[row, "R1 Model"] <- R1
      models_BIC[row, "R2 Model"] <- R2
      models_BIC[row, "B Model"]<- B_3
      models_BIC[row, "R1 BIC"] <- round(BIC(R1_fit), 2)
      models_BIC[row, "R2 BIC"] <- round(BIC(R2_fit), 2)
      models_BIC[row, "B BIC"]<- round(BIC(B_3_fit), 2)
      models_BIC[row, "Joint BIC"] <- round(sum(BIC(R1_fit), BIC(R2_fit), BIC(B_3_fit)), 2)
      
      # Move on to next row
      row <- row + 1
    }
  }
}

print(model_populations)
```

    ## [1] "UK"      "ENW"     "sim_ENW"

``` r
print(models_BIC)
```

    ##       R1 Model R2 Model B Model R1 BIC     R2 BIC     B BIC      Joint BIC  
    ##  [1,] "M7"     "M7"     "M7"    "21511.29" "21275.48" "18479.88" "61266.65" 
    ##  [2,] "M7"     "M7"     "M6"    "21511.29" "21275.48" "18171.13" "60957.9"  
    ##  [3,] "M7"     "M7"     "M5"    "21511.29" "21275.48" "17670.68" "60457.45" 
    ##  [4,] "M7"     "M6"     "M7"    "21511.29" "20926.33" "18479.88" "60917.51" 
    ##  [5,] "M7"     "M6"     "M6"    "21511.29" "20926.33" "18172.9"  "60610.52" 
    ##  [6,] "M7"     "M6"     "M5"    "21511.29" "20926.33" "17672.45" "60110.08" 
    ##  [7,] "M7"     "M5"     "M7"    "21511.29" "20397.71" "18479.88" "60388.89" 
    ##  [8,] "M7"     "M5"     "M6"    "21511.29" "20397.71" "18172.9"  "60081.9"  
    ##  [9,] "M7"     "M5"     "M5"    "21511.29" "20397.71" "17671.26" "59580.26" 
    ## [10,] "M6"     "M7"     "M7"    "21871.77" "21275.48" "18479.88" "61627.13" 
    ## [11,] "M6"     "M7"     "M6"    "21871.77" "21275.48" "18171.13" "61318.38" 
    ## [12,] "M6"     "M7"     "M5"    "21871.77" "21275.48" "17670.68" "60817.93" 
    ## [13,] "M6"     "M6"     "M7"    "21871.77" "21534.01" "18479.88" "61885.66" 
    ## [14,] "M6"     "M6"     "M6"    "21871.77" "21534.01" "18265.59" "61671.37" 
    ## [15,] "M6"     "M6"     "M5"    "21871.77" "21534.01" "17765.1"  "61170.87" 
    ## [16,] "M6"     "M5"     "M7"    "21871.77" "21004.6"  "18479.88" "61356.25" 
    ## [17,] "M6"     "M5"     "M6"    "21871.77" "21004.6"  "18265.59" "61141.96" 
    ## [18,] "M6"     "M5"     "M5"    "21871.77" "21004.6"  "17764.28" "60640.65" 
    ## [19,] "M5"     "M7"     "M7"    "43610.44" "21275.48" "18479.88" "83365.8"  
    ## [20,] "M5"     "M7"     "M6"    "43610.44" "21275.48" "18171.13" "83057.05" 
    ## [21,] "M5"     "M7"     "M5"    "43610.44" "21275.48" "17670.68" "82556.6"  
    ## [22,] "M5"     "M6"     "M7"    "43610.44" "21534.01" "18479.88" "83624.33" 
    ## [23,] "M5"     "M6"     "M6"    "43610.44" "21534.01" "18265.59" "83410.04" 
    ## [24,] "M5"     "M6"     "M5"    "43610.44" "21534.01" "17765.1"  "82909.55" 
    ## [25,] "M5"     "M5"     "M7"    "43610.44" "41703.72" "18479.88" "103794.04"
    ## [26,] "M5"     "M5"     "M6"    "43610.44" "41703.72" "18265.59" "103579.75"
    ## [27,] "M5"     "M5"     "M5"    "43610.44" "41703.72" "21944.07" "107258.23"

### Select model specifications from BIC table

Users can apply any of the criteria above to determine the most
appropriate model specification for each population. In this case, the
M7–M5–M5 model is selected as it achieves the lowest BIC across all
populations. However, it should be noted that, depending on the
population, the hierarchical approach does not always yield the lowest
BIC for subsequent populations, and this warrants further investigation.

``` r
R1_model  <- M7
R2_model  <- M5
B_3_model <- M5
```

## 2.2 - Fitted Model Plots

The models are refitted based on the selected populations (R1_model,
R2_model, and B_3_model). These fitted models generate time series of
parameter estimates, which vary according to the chosen model
specification.

### R1 model fit and residuals plot

``` r
if (identical(R1_model, M7)) {
  # M7 Fit
  hide_output_dummy <- capture.output(
    R1_fit <- fit(R1_model, data = R1_data)
  )

  # Storing R1 fitted value
  R1_fitted <- fitted(R1_fit)

  # Cohort effect
  par(mfrow = c(1, 2))
  plot(seq(min(year_range - age_upper), max(year_range - age_lower)), R1_fit$gc, 
       type = "l", xlab = "t-x", ylab = "gamma", main = paste0("R1: ", R1_data$label))
  plot(seq(min(year_range - age_upper) + 1, max(year_range - age_lower)), diff(R1_fit$gc), 
       type = "l", xlab = "t-x", ylab = "1st differenced gamma", main = paste0("R1: ", R1_data$label))

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, R1_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("R1: ", R1_data$label)) # level
  plot(year_range, R1_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("R1: ", R1_data$label)) # slope
  plot(year_range, R1_fit$kt[3, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 3", main = paste0("R1: ", R1_data$label)) # curvature

  # residuals plot
  plot(residuals(R1_fit), main = paste0("R1: ", R1_data$label))
  
} else if (identical(R1_model, M6)) {
  # M6 fit
  hide_output_dummy <- capture.output(
    R1_fit <- fit(R1_model, data = R1_data)
  )

  # Storing R1 fitted value
  R1_fitted <- fitted(R1_fit)

  # Cohort effect
  par(mfrow = c(1, 2))
  plot(seq(min(year_range - age_upper), max(year_range - age_lower)), R1_fit$gc, 
       type = "l", xlab = "t-x", ylab = "gamma", main = paste0("R1: ", R1_data$label))
  plot(seq(min(year_range - age_upper) + 1, max(year_range - age_lower)), diff(R1_fit$gc), 
       type = "l", xlab = "t-x", ylab = "1st differenced gamma", main = paste0("R1: ", R1_data$label))

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, R1_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("R1: ", R1_data$label)) # level
  plot(year_range, R1_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("R1: ", R1_data$label)) # slope

  # residuals plot
  plot(residuals(R1_fit), main = paste0("R1: ", R1_data$label))
} else if (identical(R1_model, M5)) {
  # M5 fit
  hide_output_dummy <- capture.output(
    R1_fit <- fit(R1_model, data = R1_data)
  )

  # Storing R1 fitted value
  R1_fitted <- fitted(R1_fit)

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, R1_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("R1: ", R1_data$label)) # level
  plot(year_range, R1_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("R1: ", R1_data$label)) # slope

  # residuals plot
  plot(residuals(R1_fit), main = paste0("R1: ", R1_data$label))
} else {
  print("R1 Model specification not included/correct")
}
```

![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

### R2 model fit and residuals plot

``` r
if (identical(R2_model, M7)) {
  # M7 Fit
  hide_output_dummy <- capture.output(
    R2_fit <- fit(R2_model, data = R2_data, oxt = R1_fitted) # oxt is the offset)
  )

  # Storing R2 fitted value
  R2_fitted <- fitted(R2_fit)

  # Cohort effect
  par(mfrow = c(1, 2))
  plot(seq(min(year_range - age_upper), max(year_range - age_lower)), R2_fit$gc, 
       type = "l", xlab = "t-x", ylab = "gamma", main = paste0("R2: ", R2_data$label))
  plot(seq(min(year_range - age_upper) + 1, max(year_range - age_lower)), diff(R2_fit$gc), 
       type = "l", xlab = "t-x", ylab = "1st differenced gamma", main = paste0("R2: ", R2_data$label))

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, R2_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("R2: ", R2_data$label)) # level
  plot(year_range, R2_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("R2: ", R2_data$label)) # slope
  plot(year_range, R2_fit$kt[3, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 3", main = paste0("R2: ", R2_data$label)) # curvature

  # residuals plot
  plot(residuals(R2_fit), main = paste0("R2: ", R2_data$label))
  
} else if (identical(R2_model, M6)) {
  # M6 fit
  hide_output_dummy <- capture.output(
    R2_fit <- fit(R2_model, data = R2_data, oxt = R1_fitted) # oxt is the offset
  )

  # Storing R2 fitted value
  R2_fitted <- fitted(R2_fit)

  # Cohort effect
  par(mfrow = c(1, 2))
  plot(seq(min(year_range - age_upper), max(year_range - age_lower)), R2_fit$gc, 
       type = "l", xlab = "t-x", ylab = "gamma", main = paste0("R2: ", R2_data$label))
  plot(seq(min(year_range - age_upper) + 1, max(year_range - age_lower)), diff(R2_fit$gc), 
       type = "l", xlab = "t-x", ylab = "1st differenced gamma", main = paste0("R2: ", R2_data$label))

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, R2_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("R2: ", R2_data$label)) # level
  plot(year_range, R2_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("R2: ", R2_data$label)) # slope

  # residuals plot
  plot(residuals(R2_fit), main = paste0("R2: ", R2_data$label))
} else if (identical(R2_model, M5)) {
  # M5 fit
  hide_output_dummy <- capture.output(
    R2_fit <- fit(R2_model, data = R2_data, oxt = R1_fitted) # oxt is the offset
  )

  # Storing R2 fitted value
  R2_fitted <- fitted(R2_fit)

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, R2_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("R2: ", R2_data$label)) # level
  plot(year_range, R2_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("R2: ", R2_data$label)) # slope

  # residuals plot
  plot(residuals(R2_fit), main = paste0("R2: ", R2_data$label))
} else {
  print("R2 Model specification not included/correct")
}
```

![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

### B model fit and residuals plot

``` r
if (identical(B_3_model, M7)) {
  # M7 Fit
  hide_output_dummy <- capture.output(
    B_3_fit <- fit(B_3_model, data = B_data, oxt = R2_fitted) # oxt is the offset)
  )

  # Storing B_3 fitted value
  b_3_fitted <- fitted(B_3_fit)

  # Cohort effect
  par(mfrow = c(1, 2))
  plot(seq(min(year_range - age_upper), max(year_range - age_lower)), B_3_fit$gc, 
       type = "l", xlab = "t-x", ylab = "gamma", main = paste0("B: ", B_data$label))
  plot(seq(min(year_range - age_upper) + 1, max(year_range - age_lower)), diff(B_3_fit$gc), 
       type = "l", xlab = "t-x", ylab = "1st differenced gamma", main = paste0("B: ", B_data$label))


  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, B_3_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("B: ", B_data$label)) # level
  plot(year_range, B_3_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("B: ", B_data$label)) # slope
  plot(year_range, B_3_fit$kt[3, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 3", main = paste0("B: ", B_data$label)) # curvature


  # residuals plot
  plot(residuals(B_3_fit), main = paste0("B_3: ", B_data$label))
  
} else if (identical(B_3_model, M6)) {
  # M6 fit
  hide_output_dummy <- capture.output(
    B_3_fit <- fit(B_3_model, data = B_data, oxt = R2_fitted) # oxt is the offset
  )

  # Storing B_3 fitted value
  B_3_fitted <- fitted(B_3_fit)

  # Cohort effect
  par(mfrow = c(1, 2))
  plot(seq(min(year_range - age_upper), max(year_range - age_lower)), B_3_fit$gc, 
       type = "l", xlab = "t-x", ylab = "gamma", main = paste0("B: ", B_data$label))
  plot(seq(min(year_range - age_upper) + 1, max(year_range - age_lower)), diff(B_3_fit$gc), 
       type = "l", xlab = "t-x", ylab = "1st differenced gamma", main = paste0("B: ", B_data$label))

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, B_3_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("B: ", B_data$label)) # level
  plot(year_range, B_3_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("B: ", B_data$label)) # slope

  # residuals plot
  plot(residuals(B_3_fit), main = paste0("B_3: ", B_data$label))
} else if (identical(B_3_model, M5)) {
  # M5 fit
  hide_output_dummy <- capture.output(
    B_3_fit <- fit(B_3_model, data = B_data, oxt = R2_fitted) # oxt is the offset
  )

  # Storing B_3 fitted value
  B_3_fitted <- fitted(B_3_fit)

  # age effects with time
  par(mfrow = c(1, 3))
  plot(year_range, B_3_fit$kt[1, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 1", main = paste0("B: ", B_data$label)) # level
  plot(year_range, B_3_fit$kt[2, ], type = "l", 
       xlab = "Calendar Year (t)", ylab = "kappa 2", main = paste0("B: ", B_data$label)) # slope

  # residuals plot
  plot(residuals(B_3_fit), main = paste0("B: ", B_data$label))
} else {
  print("B_3 Model specification not included/correct")
}
```

![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](2---Model-Fitting-and-Plots_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

### R1 model plots

``` r
# Open a larger window (width & height in inches)
dev.new(width = 12, height = 10)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3, 1))  # margins: bottom, left, top, right

# Plotting R1 fit
plot_ages <- c("60", "70", "80", "90") # Change plot ages here
calendar_year <- as.numeric(colnames(R1_qxt))

# par(mfrow = c(2, 2))
for (i in 1:length(plot_ages)) {
  age <- plot_ages[i]
  # Plot raw mortality rate vs fitted
  plot(calendar_year, R1_qxt[age, ], 
       main = paste("R1: ", R1_data$label, " Population Mortality \n vs Fitted (Age ", age, ")", sep = ""),
       xlab = "Calendar Year (t)", ylab = paste(R1_data$label, " q", age, ",t", sep = ""))
  lines(calendar_year, inv_logit(R1_fitted[age, ]), type = "l")
}
```

### R1 model plots

``` r
# Open a larger window (width & height in inches)
dev.new(width = 12, height = 10)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3, 1))  # margins: bottom, left, top, right

# Plotting R2 fit
plot_ages <- c("60", "70", "80", "90") # Change plot ages here
calendar_year <- as.numeric(colnames(R2_qxt))

#par(mfrow = c(2, 2))
for (i in 1:length(plot_ages)) {
  age <- plot_ages[i]
  # Plot raw mortality rate vs fitted
  plot(calendar_year, R2_qxt[age, ], 
       main = paste("R2: ", R2_data$label, " Population Mortality \n vs Fitted (Age ", age, ")", sep = ""),
       xlab = "Calendar Year (t)", ylab = paste(R2_data$label, " q", age, ",t", sep = ""))
  lines(calendar_year, inv_logit(R2_fitted[age, ]), type = "l")
}
```

### B model plots

``` r
# Open a larger window (width & height in inches)
dev.new(width = 12, height = 10)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3, 1))  # margins: bottom, left, top, right

# Plotting B_3 fit
plot_ages <- c("60", "70", "80", "90") # Change plot ages here
calendar_year <- as.numeric(colnames(B_qxt))

# par(mfrow = c(1, 1))
for (i in 1:length(plot_ages)) {
  age <- plot_ages[i]
  # Plot raw mortality rate vs fitted
  plot(calendar_year, B_qxt[age, ], 
       main = paste("B: ", B_data$label, " Population Mortality \n vs Fitted (Age ", age, ")", sep = ""),
       xlab = "Calendar Year (t)", ylab = paste(B_data$label, " q", age, ",t", sep = ""))
  lines(calendar_year, inv_logit(B_3_fitted[age, ]), type = "l")
}
```

### BIC values for selected models

``` r
CBD_BIC <- matrix(nrow = 3, ncol = 1)
rownames(CBD_BIC) <- c("R1", "R2", "B")
colnames(CBD_BIC) <- c("BIC")

CBD_BIC["R1", "BIC"] <- BIC(R1_fit)
CBD_BIC["R2", "BIC"] <- BIC(R2_fit)
CBD_BIC["B", "BIC"] <- BIC(B_3_fit)

print(CBD_BIC)
```

    ##         BIC
    ## R1 21511.29
    ## R2 20397.71
    ## B  17671.26

## Save Environment

``` r
save.image(file = "environment.RData")
```
