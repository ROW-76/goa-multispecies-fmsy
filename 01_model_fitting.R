# Model fitting sequence for Gulf of Alaska groundfish.
# Fits single-species, single-species with estimated natural mortality,
# and multi-species (predation-explicit) models.
#
# Data: GOA2018SS (walleye pollock, Pacific cod, arrowtooth flounder), 1977-2018.
# Package: Rceattle (Adams et al., 2022, Fisheries Research 251).

library(Rceattle)

data("GOA2018SS")
GOA2018SS$fleet_control$proj_F_prop <- rep(1, nrow(GOA2018SS$fleet_control))

# Single-species, fixed natural mortality
ss_run <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = NULL,
  file = NULL,
  estimateMode = 0,
  random_rec = FALSE,
  msmMode = 0,
  phase = "default",
  initMode = 2,
  verbose = 1
)

# Single-species, estimated natural mortality
ss_run_M <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = NULL,
  file = NULL,
  estimateMode = 0,
  M1Fun = build_M1(M1_model = c(1, 2, 1)),
  random_rec = FALSE,
  msmMode = 0,
  phase = "default",
  initMode = 2,
  verbose = 1
)

# Multi-species, predation mortality estimated (MSVPA-based, empirical suitability)
ms_run <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = ss_run_M$estimated_params,
  file = NULL,
  estimateMode = 0,
  M1Fun = build_M1(M1_model = c(1, 2, 1)),
  niter = 3,
  random_rec = FALSE,
  msmMode = 1,
  suitMode = 0,
  initMode = 2,
  phase = "default",
  verbose = 1
)

ss_run$opt$convergence
ss_run_M$opt$convergence
ms_run$opt$convergence

mod_list <- list(ss_run, ss_run_M, ms_run)
mod_names <- c("SS", "SS-M", "MS")

plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
