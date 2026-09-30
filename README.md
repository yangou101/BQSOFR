# GAL Project: Bayesian Functional Quantile Regression with Measurement Error

[![R](https://img.shields.io/badge/Language-R-blue)](https://www.r-project.org/)
[![Stan](https://img.shields.io/badge/Computation-Stan-orange)](https://mc-stan.org/)
[![GAL](https://img.shields.io/badge/Likelihood-GAL-purple)](#methods)
[![NHANES](https://img.shields.io/badge/Data-NHANES-brightgreen)](#data)
[![Status](https://img.shields.io/badge/Status-Research%20Project-success)](#overview)

This repository contains analysis code for the GAL project, which studies Bayesian scalar-on-function quantile regression with measurement error using the generalized asymmetric Laplace (GAL) distribution.

The application examines the association between weekday physical activity profiles and conditional quantiles of body mass index (BMI) using NHANES data.

---

## Overview

Repeated measurements of physical activity contain both information about an individual's usual activity pattern and measurement variability. This project studies how to account for measurement error when relating a functional exposure to different conditional quantiles of a scalar outcome.

The main analysis tasks are to:

- Prepare and align the NHANES outcome, scalar covariates, and repeated functional measurements.
- Construct functional exposure representations and measurement-error-corrected curves.
- Fit Bayesian functional quantile regression models using a GAL working likelihood.
- Compare fully Bayesian joint modelling, two-stage regression calibration, and naive analysis.
- Summarize scalar and functional regression coefficients, evaluate predictive performance, and generate figures.

---

## Repository Structure

The current documentation focuses on the NHANES application.

| Location | Contents |
| --- | --- |
| `Application/` | Application code, including function definitions and analysis scripts |
| `README.md` | Project overview, data download instructions, and analysis workflow |
| Release `data-v1` | The separately uploaded `GAL.RData` data file |

The data file is distributed through GitHub Releases rather than stored in the ordinary repository file tree.

---

## Data

The application uses NHANES data to study physical activity profiles and BMI, with adjustment for demographic, dietary, and health-related covariates.

The analysis data are distributed in `GAL.RData`. Because the file exceeds GitHub's size limit for ordinary Git-tracked files, it is provided as a **Release asset**.

### Download the Data

1. Open **Releases** on the repository homepage.
2. Select the release titled **GAL project data**, associated with the tag **`data-v1`**.
3. Expand **Assets**.
4. Click **`GAL.RData`** to download the data.

**Download `GAL.RData` itself.** The automatically generated **Source code (zip)** and **Source code (tar.gz)** archives do not include this separately uploaded data asset.

Download the analysis code from the **main** branch to obtain the current scripts.

---

## Methods

The project uses the GAL distribution as a flexible working likelihood for Bayesian scalar-on-function quantile regression.

The study considers three approaches:

### Fully Bayesian Joint Modelling (FBQ)

The joint approach models the repeated functional measurements and the outcome together, treating the underlying functional exposure as latent and accounting for uncertainty within the joint model.

### Two-Stage Regression Calibration (RC)

The two-stage approach first estimates the underlying functional exposure using a mixed-effects-model-based correction. The estimated curves are then treated as fixed inputs to Bayesian functional quantile regression. The application code uses FUI-corrected curves for this approach.

### Naive Analysis

The naive approach uses averaged observed functional measurements without an explicit measurement-error correction and serves as a comparison.

---

## Analysis Workflow

### 1. Download the Code and Data

Download or clone the code from the **main** branch and download `GAL.RData` separately from the `data-v1` release.

### 2. Set Up the R Environment

Install the packages used by the relevant scripts and configure the Stan/RStan toolchain before fitting models.

### 3. Load Function Definitions and Prepare the Data

Load the function definitions required by the application scripts. Update the data-loading path to the location of the downloaded file, for example:

```r
# Replace this example path with the actual location of GAL.RData.
load("path/to/GAL.RData")
```

Run the data preparation code and ensure that the outcome, scalar covariates, and functional measurements use the same subjects and ordering.

### 4. Fit the Models

Run the relevant application scripts in `Application/` for the desired methods and quantile levels.

### 5. Summarize and Compare Results

Use the corresponding post-processing code to summarize posterior estimates, plot functional coefficients, and calculate LOO where pointwise log-likelihood draws are available.

---

## Software Requirements

The analysis uses **R** and **Stan through RStan**.

Packages used by the supplied model-fitting and post-processing code include:

```r
install.packages(c(
  "rstan",
  "brms",
  "refund",
  "mgcv",
  "loo",
  "ggplot2"
))
```

Additional dependencies may be required by individual data preparation or first-stage correction scripts. Check their package-loading statements before running them.

RStan also requires a compatible C++ toolchain.

---

## Reproducibility Notes

- Update local file paths before running the scripts.
- Load the required function definitions before executing model-fitting calls.
- Keep subject selection, subject ordering, and time-grid ordering consistent across models.
- For LOO comparisons, use the same observations, quantile level, and response scale.
- Save fitted objects and record the random seed, sampling settings, and package versions used for each analysis.

---

## Contact

**Yang Ou**  
Indiana University Bloomington  
Email: [yangou@iu.edu](mailto:yangou@iu.edu)
