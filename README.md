[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19284370.svg)](https://doi.org/10.5281/zenodo.19284370)

# PSL: Code and Data for Spatial Model Comparison

## Overview
This repository provides the dataset and R script required to reproduce the spatial model comparisons presented in our manuscript. It evaluates 6 spatial models: **LM, RF, BCS, GOS, PS, and PSL (Two-Stream Architecture)**, specifically applied to the spatial prediction of Soil Organic Carbon (SOC).

## Repository Contents
* **`data.csv`**: The dataset containing spatial coordinates, cross-validation fold IDs, the response variable (SOC), 10 original covariates, and 40 expanded features.
* **`6 Model Calculation Code.R`**: The main R script designed to train, evaluate, and compare all 6 models. Hyperparameters and random seeds (seed=42) are strictly aligned with the manuscript for full reproducibility.

## Dependencies
To execute the code, please ensure you have **R** installed, along with the following packages:
```R
install.packages(c("readr", "dplyr", "tibble", "geosimilarity", "ranger"))
