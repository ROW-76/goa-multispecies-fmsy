# Harvest control rule comparison and multi-species FMSY estimation.
# Requires ss_run_M from 01_model_fitting.R as initial values.
#
# Four rules are compared:
#   F40               constant F achieving 40% of unfished SSB (target-based)
#   F40, iterative    as above, with predators projected before prey
#   CMSY              unconstrained catch maximization
#   Constrained CMSY  catch maximization with a 35% SSB floor
#
# Each model is projected from 2019 to 2100 under its resulting constant F.

library(Rceattle)

data("GOA2018SS")
GOA2018SS$projyr <- 2100
GOA2018SS$fleet_control$proj_F_prop <- rep(1, nrow(GOA2018SS$fleet_control))

# F40: SB0 derived by projecting all species simultaneously under no fishing
ms_run_f40 <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = ss_run_M$estimated_params,
  file = NULL,
  estimateMode = 0,
  niter = 3,
  random_rec = FALSE,
  HCR = build_hcr(HCR = 3, DynamicHCR = FALSE, Ftarget = 0.4),
  msmMode = 1,
  suitMode = 0,
  verbose = 1
)

# F40, iterative: arrowtooth and cod projected under no fishing first,
# then pollock projected under no fishing with predators at F40
ms_run_f40_iter <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = ss_run_M$estimated_params,
  file = NULL,
  estimateMode = 0,
  niter = 3,
  random_rec = FALSE,
  HCR = build_hcr(HCR = 3, DynamicHCR = FALSE, Ftarget = 0.4, HCRorder = c(2, 1, 1)),
  msmMode = 1,
  suitMode = 0,
  verbose = 1
)

# Unconstrained CMSY
ms_run_cmsy <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = ss_run_M$estimated_params,
  file = NULL,
  estimateMode = 0,
  niter = 3,
  random_rec = FALSE,
  HCR = build_hcr(HCR = 1),
  msmMode = 1,
  suitMode = 0,
  verbose = 1
)

# Constrained CMSY: 35% SSB floor
ms_run_concmsy <- Rceattle::fit_mod(
  data_list = GOA2018SS,
  inits = ss_run_M$estimated_params,
  file = NULL,
  estimateMode = 0,
  niter = 3,
  random_rec = FALSE,
  HCR = build_hcr(HCR = 1, Plimit = 0.35),
  msmMode = 1,
  suitMode = 0,
  verbose = 1
)

mod_list <- list(ms_run_f40, ms_run_f40_iter, ms_run_cmsy, ms_run_concmsy)
model_names <- c("F40", "F40, iterative", "CMSY", "Constrained CMSY")

plot_biomass(mod_list, model_names = model_names, incl_proj = TRUE)
plot_ssb(mod_list, model_names = model_names, incl_proj = TRUE)
plot_recruitment(mod_list, model_names = model_names, incl_proj = TRUE)
plot_catch(mod_list, model_names = model_names, incl_proj = TRUE)

# Alternative verification: uniform F tested across all species by direct
# grid search, used to cross-check the optimizer-based CMSY result.
f_values <- seq(0.1, 1.5, by = 0.05)
grid_results <- vector("list", length(f_values))
yields <- numeric(length(f_values))

for (i in seq_along(f_values)) {
  grid_results[[i]] <- try(
    Rceattle::fit_mod(
      data_list = GOA2018SS,
      inits = ss_run_M$estimated_params,
      file = NULL,
      estimateMode = 0,
      niter = 3,
      random_rec = FALSE,
      HCR = build_hcr(HCR = 5, Ftarget = rep(f_values[i], 3)),
      M1Fun = build_M1(M1_model = c(1, 2, 1)),
      msmMode = 1,
      suitMode = 0,
      verbose = 0
    ),
    silent = TRUE
  )

  if (!inherits(grid_results[[i]], "try-error")) {
    proj_years <- grid_results[[i]]$data_list$projyr
    yields[i] <- sum(grid_results[[i]]$quantities$catch[, proj_years])
  } else {
    yields[i] <- NA
  }
}

fmsy_index <- which.max(yields)
fmsy_grid <- f_values[fmsy_index]

plot(f_values, yields, type = "b",
     xlab = "Fishing mortality (F)",
     ylab = "Total yield (mt)",
     main = "Yield curve, uniform F across species")
abline(v = fmsy_grid, col = "red", lty = 2)
