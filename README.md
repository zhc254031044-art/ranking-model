# Statistical Models for Ranking Data

This repository contains the R code used for the MSc dissertation Statistical Models for Ranking Data.

The repository is organised into two folders:

## data_analysis

Contains the APA 2008 election data and the empirical analysis code.

- `APA.soi` — APA 2008 election dataset.
- `Data_analysis.R` — descriptive analysis and figures.
- `APA_Thurstone_BT.R` — Thurstone and Bradley--Terry analysis.
- `APA_PL_MM.R` — Plackett--Luce and Mallows analysis.
- `APA_PL_BayesianPL.R` — classical and Bayesian Plackett--Luce analysis.

## simulation

Contains the simulation code used in the dissertation.

- `Alpha for Mallows.R` — calibration of Mallows concentration parameters.
- `Pairwise Data Generation BT mechanism.R` — pairwise data generated under the Bradley--Terry mechanism.
- `Pairwise Data Generation Thurstone mechanism.R` — pairwise data generated under the Thurstone mechanism.
- `Complete Data Generation PL mechanism.R` — complete rankings generated under the Plackett--Luce mechanism.
- `Complete Data Generation MM mechanism.R` — complete rankings generated under the Mallows mechanism.
- `Bayesian PL - PL mechanism.R` — Bayesian Plackett--Luce analysis using partial rankings generated from the Plackett--Luce mechanism.
- `Bayesian PL - Mallows mechanism.R` — Bayesian Plackett--Luce analysis using partial rankings generated from the Mallows mechanism.

## Software

All analyses were implemented in R.
