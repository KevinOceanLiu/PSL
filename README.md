[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19284370.svg)](https://doi.org/10.5281/zenodo.19284370)

# PSL: Code and Data for Spatial Model Comparison

## Overview

This repository provides the dataset and R script required to reproduce the spatial model comparisons presented in our manuscript. It evaluates six spatial prediction models: **SGWR, GOS, LM, RF, PS, and PSL**. The case study focuses on the spatial prediction of Soil Organic Carbon (SOC).

The script first constructs a 40-dimensional spatial pattern feature space from the original 10 environmental covariates, and then evaluates the six models under the same cross-validation setting.

## Repository Contents

- **`data.csv`**  
  Input dataset containing spatial coordinates, cross-validation fold IDs, the response variable SOC, and 10 original environmental covariates.

- **`6_model_calculation_code.R`**  
  Main R script for constructing 40D spatial pattern features, training the six models, generating predictions and residuals, and reporting model performance metrics. The random seed is fixed at `42` for reproducibility.

- **`6models_residuals.csv`**  
  Output file generated after running the script. It contains observed SOC values, model predictions, and residuals for all six models.

## Required Data Columns

The input file `data.csv` should include the following columns:

```text
longitude, latitude, SOC, hex_id, foldID,
Slope, Biomass, Precipitation, TPI, Temperature,
NPP, SolarRadiation, WindSpeed, SinAspect, CosAspect
