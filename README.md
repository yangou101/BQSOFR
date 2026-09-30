# GAL Project: Bayesian Scalar-on-Function Quantile Regression with Measurement Error

This repository provides R and Stan code for the GAL project, which studies Bayesian scalar-on-function quantile regression with measurement error using the generalized asymmetric Laplace (GAL) distribution.

## Methods

The study considers:

- **Fully Bayesian joint modelling (FBQ):** jointly models the outcome and repeated functional measurements to account for measurement error.
- **Two-stage regression calibration (RC):** estimates the underlying functional exposure using a mixed-effects-model-based correction, followed by Bayesian quantile regression. The application code uses FUI-corrected curves.
- **Naive analysis:** uses uncorrected averaged functional measurements as a comparison.

## NHANES Application

The application examines associations between weekday physical activity profiles and conditional quantiles of body mass index (BMI) using NHANES data, with adjustment for demographic, dietary, and health-related covariates.

Application scripts are located in the `Application/` folder.

## Data Download

The data file, `GAL.RData`, is provided separately as a GitHub Release asset because it exceeds the size limit for ordinary Git-tracked files.

To download the data:

1. Open **Releases** on the repository homepage.
2. Select **GAL project data**, associated with the tag **`data-v1`**.
3. Expand **Assets** and click **`GAL.RData`**.

Download the RData asset itself. The automatically generated **Source code (zip)** and **Source code (tar.gz)** archives do not include this separately uploaded data asset.

## Running the Analysis

1. Download or clone the analysis code from the **main** branch.
2. Download `GAL.RData` from the `data-v1` release.
3. Install the R packages required by the analysis scripts and configure the Stan/RStan toolchain.
4. Update the data-loading path in the scripts to the local location of `GAL.RData`.
5. Load the function definitions, prepare the analysis data, and run the application scripts in `Application/`.

The supplied model-fitting functions use packages including `rstan`, `brms`, `refund`, and `mgcv`. LOO evaluation and plotting additionally use `loo` and `ggplot2`; consult the relevant scripts for other dependencies.

## Contact

Yang Ou  
Indiana University Bloomington  
Email: yangou@iu.edu

