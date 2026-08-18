# Multi-Species FMSY for Gulf of Alaska Groundfish

This repository contains the code and report for a project calculating multi-species fishing mortality at maximum sustainable yield (FMSY) using the Rceattle model, applied to three interacting Gulf of Alaska groundfish species.

Advanced Population Modelling course project, December 2025.

## Overview

Fishing mortality that produces maximum sustainable yield (FMSY) can be calculated analytically in single-species models using per-recruit analysis. This analytical approach does not hold once predation mortality depends on the abundance of multiple interacting species. This project addresses that problem for three Gulf of Alaska groundfish species: walleye pollock (*Gadus chalcogrammus*), arrowtooth flounder (*Atheresthes stomias*), and Pacific cod (*Gadus macrocephalus*).

Four harvest control rules are compared by projecting the fitted model forward to equilibrium (2018-2100):

- F40: constant fishing mortality achieving 40% of unfished spawning biomass
- F40, iterative: as above, with predator species projected before prey
- CMSY: unconstrained catch maximization
- Constrained CMSY: catch maximization with a 35% spawning biomass floor

Unconstrained catch maximization produces biological collapse across all three species. Constrained catch maximization with a 35% biomass floor produces stable, biologically viable multi-species FMSY estimates.

## Data

`GOA2018SS`, distributed with the Rceattle package. Age-structured catch and survey data, biological parameters, and stomach content data for predation estimation, 1977-2018.

## Model

Rceattle (R-based Climate-Enhanced Age-based model with Temperature-specific Trophic Linkages and Energetics), implemented in Template Model Builder (TMB) (Holsman et al., 2016; Adams et al., 2022). Predation mortality is estimated as a function of predator and prey abundance, size, and suitability, and is partitioned from fishing mortality and baseline natural mortality.

Models were fit in sequence: single-species with fixed natural mortality, single-species with estimated natural mortality, then multi-species with estimated predation mortality (`msmMode = 1`, `suitMode = 0`, `niter = 3`).

## Results

Multi-species FMSY, constrained CMSY (35% biomass floor):

| Species | FMSY | MSY (kt) | SSB<sub>MSY</sub> (kt) | Depletion | F40 |
|---|---|---|---|---|---|
| Pollock | 0.372 | 3,353 | 2,166 | 35.0% | 0.317 |
| Arrowtooth flounder | 0.696 | 194 | 270 | 35.0% | 0.542 |
| Pacific cod | 0.270 | 135 | 123 | 35.0% | 0.189 |
| Total | - | 3,682 | 2,559 | - | - |

Constrained FMSY values are 17-43% higher than the corresponding F40 values, reflecting the lower equilibrium biomass permitted by the 35% constraint relative to the 40% target. Fishing mortality variance across the final ten projection years (2091-2100) is on the order of 10<sup>-27</sup> for all species, confirming equilibrium.

The full methodology, validation, and discussion are in `docs/report.pdf`.

## Repository Structure

```
R/
  01_model_fitting.R           single-species and multi-species model fitting
  02_harvest_control_rules.R   F40, F40 iterative, CMSY, and constrained CMSY
  03_analytical_fmsy_check.R   single-species per-recruit cross-check
docs/
  report.pdf                   full written report
  index.html                   project summary page
README.md
```

`03_analytical_fmsy_check.R` implements the standard single-species per-recruit method as a methodological control. It does not account for predation mortality and is included to demonstrate why this approach is not valid for a multi-species system with trophic interactions, not as an alternative estimate of FMSY.

## Requirements

- R (>= 4.0)
- [Rceattle](https://github.com/grantdadams/Rceattle)
- TMB

Install Rceattle from source:

```r
install.packages("devtools")
install.packages("TMB", type = "source")
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
devtools::install_github("grantdadams/Rceattle")
```

## Running the Analysis

Scripts are run in order. `02_harvest_control_rules.R` and `03_analytical_fmsy_check.R` each require an Rceattle session with the relevant objects from `01_model_fitting.R` (`ss_run_M`) already fit.

```r
source("R/01_model_fitting.R")
source("R/02_harvest_control_rules.R")
source("R/03_analytical_fmsy_check.R")
```

## References

Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn, M. W., Ianelli, J. N., Spies, I., ... Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. *Fisheries Research*, 251, 106303.

Gaichas, S., Gamble, R., Fogarty, M., Benoit, H., Essington, T., Fu, C., ... Link, J. (2012). Assembly rules for aggregate-species production models: simulations in support of management strategy evaluation. *Marine Ecology Progress Series*, 459, 275-292.

Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2016). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. *Deep Sea Research Part II*, 134, 360-378.
