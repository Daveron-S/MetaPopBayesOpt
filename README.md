# Hybrid SDM–metapopulation modelling

This repository contains simulation analyses used to investigate how sampling and parameter uncertainty affect inference in a hybrid species distribution model (SDM)–metapopulation framework.

The analyses examine how well a known metapopulation occupancy surface can be reconstructed from incomplete occurrence data, how reliably the parameters underlying that surface can be recovered, and the extent to which different parameter combinations can produce similar spatial predictions.

## Overview

The repository contains two complementary simulation analyses.

### 1. Sampling and reconstruction of occupancy patterns

The first analysis examines how:

- incomplete occurrence information; and
- minimum spatial separation among sampling locations;

affect the ability of the hybrid SDM–metapopulation model to reconstruct a known landscape-wide occupancy surface.

The analysis also examines recovery of the fitted SDM and metapopulation parameters as sampling effort changes.

### 2. Bayesian optimisation and equifinality

The second analysis uses Bayesian optimisation to search a biologically constrained parameter space for combinations of SDM and metapopulation parameters that reproduce a known reference occupancy surface.

Candidate parameterisations are evaluated against the complete reference occupancy surface, reducing the influence of sampling variation when examining parameter recovery.

The analysis is used to assess:

- how closely Bayesian optimisation can reproduce the reference occupancy pattern;
- which parameters are consistently recovered;
- variation in parameter estimates among independent optimisation runs; and
- equifinality, where different parameter combinations generate similarly accurate occupancy predictions.

## Hybrid model

The simulations combine an SDM component with a patch-based metapopulation model.

The SDM component describes habitat suitability using responses to:

- conifer woodland (`p1`);
- broadleaf woodland (`p2`);
- vegetation height (`p3`); and
- the SDM intercept (`p4`).

The metapopulation component includes parameters describing:

- baseline extinction risk (`mu`);
- scaling of extinction with local population size (`x`);
- dispersal (`alpha`); and
- scaling of colonisation with connectivity (`ygamma`).

Together, these components allow predicted occurrence to reflect both environmental suitability and metapopulation processes.

## Bayesian optimisation

Bayesian optimisation is implemented using `rBayesianOptimization`.

Each independent optimisation consists of:

- 50 randomly selected initial parameter combinations;
- 50 sequential optimisation iterations;
- an upper-confidence-bound acquisition function (`kappa = 2.576`);
- a Matérn covariance kernel (`nu = 2.5`).

For each candidate parameter set, the hybrid model is simulated repeatedly and the resulting mean occupancy surface is compared with the known reference surface using Warren's I.

The optimisation is repeated independently 20 times to examine variation among high-performing parameter combinations.

The highest-scoring parameter set is subsequently used in 50 independent simulations to generate a mean predicted occupancy surface for comparison with the reference distribution.

## Model evaluation

Agreement between reference and reconstructed occupancy surfaces is assessed using:

- **Warren's I**, describing spatial similarity between occupancy surfaces; and
- **Pearson's correlation coefficient**, describing correspondence in cell-level predicted occupancy.

Parameter recovery is assessed by comparing estimates from independent optimisation runs with the known generating parameter values.

## Main outputs

The analyses generate:

- reconstructed occupancy surfaces;
- parameter estimates from repeated model fits;
- summaries of optimisation performance;
- comparisons between known and reconstructed occupancy surfaces; and
- publication-ready figures showing occupancy reconstruction and variation in parameter estimates.

## R packages

The analyses use the following R packages:

- `PopScape`
- `terra`
- `rBayesianOptimization`
- `Rcpp`
- `ggplot2`
- `dplyr`
- `tidyr`
- `purrr`
- `ggtext`

Additional C++ functions used by the metapopulation model are compiled through `Rcpp`.

## Reproducibility

The simulations are stochastic. Re-running the complete workflow may therefore produce small differences in individual parameter estimates and occupancy surfaces unless random-number seeds are fixed.

The optimisation results reported in the associated analysis are based on 20 independent Bayesian optimisation runs, with the best-performing parameterisation subsequently evaluated across 50 metapopulation simulations.

## Repository purpose

This repository supports analyses investigating how sampling intensity and spatial separation affect reconstruction of a known metapopulation occupancy surface and parameter recovery in a hybrid SDM–metapopulation model, and how Bayesian optimisation can be used to examine parameter recovery and equifinality. And is used to produce the Figures in the manuscript Smith, Maclean, Baker 2026 Linking habitat suitability and metapopulation dynamics to predict species distributions.
