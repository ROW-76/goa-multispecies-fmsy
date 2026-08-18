# Single-species per-recruit FMSY calculation, used as a check against the
# numerical equilibrium projection in 02_harvest_control_rules.R.
#
# This method assumes natural mortality is independent of predator and prey
# abundance, which does not hold for this system. It is included to confirm
# that the analytical approach breaks down once predation is accounted for.

library(Rceattle)

data("GOA2018SS")

ss_run <- fit_mod(
  data_list = GOA2018SS,
  inits = NULL,
  estimateMode = 0,
  msmMode = 0,
  verbose = 1
)

ms_run <- fit_mod(
  data_list = GOA2018SS,
  inits = ss_run$estimated_params,
  estimateMode = 0,
  msmMode = 1,
  niter = 3,
  verbose = 1
)

calc_SPR <- function(F_val, M, wt, mat, sex_ratio, sel, spawn_month) {
  nages <- length(M)
  N <- rep(0, nages)
  N[1] <- sex_ratio[1]

  for (age in 2:(nages - 1)) {
    Z <- M[age - 1] + F_val * sel[age - 1]
    N[age] <- N[age - 1] * exp(-Z)
  }

  Z_minus1 <- M[nages - 1] + F_val * sel[nages - 1]
  Z_plus <- M[nages] + F_val * sel[nages]

  if (exp(-Z_plus) < 0.9999) {
    N[nages] <- N[nages - 1] * exp(-Z_minus1) / (1 - exp(-Z_plus))
  } else {
    N[nages] <- N[nages - 1] * exp(-Z_minus1) * 100
  }

  Z_vec <- M + F_val * sel
  sum(N * wt * mat * exp(-Z_vec * spawn_month / 12))
}

calc_YPR <- function(F_val, M, wt, sel) {
  nages <- length(M)
  N <- rep(0, nages)
  N[1] <- 1

  for (age in 2:(nages - 1)) {
    Z <- M[age - 1] + F_val * sel[age - 1]
    N[age] <- N[age - 1] * exp(-Z)
  }

  Z_minus1 <- M[nages - 1] + F_val * sel[nages - 1]
  Z_plus <- M[nages] + F_val * sel[nages]

  if (exp(-Z_plus) < 0.9999) {
    N[nages] <- N[nages - 1] * exp(-Z_minus1) / (1 - exp(-Z_plus))
  } else {
    N[nages] <- N[nages - 1] * exp(-Z_minus1) * 100
  }

  YPR <- 0
  for (age in 1:nages) {
    F_age <- F_val * sel[age]
    Z <- M[age] + F_age
    if (Z > 0.0001) {
      catch_age <- F_age * N[age] * (1 - exp(-Z)) / Z
      YPR <- YPR + catch_age * wt[age]
    }
  }
  YPR
}

nspp <- ms_run$data_list$nspp
species_names <- c("Pollock", "Pcod", "ATF")
results <- data.frame(Species = species_names, FMSY = NA, BMSY = NA, MSY = NA)
F_grid <- seq(0, 2.0, by = 0.01)

for (sp in 1:nspp) {

  nages <- ms_run$data_list$nages[sp]
  minage <- ms_run$data_list$minage[sp]

  M <- ms_run$quantities$M1[sp, 1, 1:nages]

  wt_index <- ms_run$data_list$ssb_wt_index[sp]
  wt_rows <- which(ms_run$data_list$wt$Wt_index == wt_index)
  wt_data <- ms_run$data_list$wt[wt_rows[1], ]
  wt <- as.numeric(wt_data[, 6:(6 + nages - 1)])

  mat <- as.numeric(ms_run$data_list$pmature[sp, 2:(1 + nages)])
  sex_ratio <- as.numeric(ms_run$data_list$sex_ratio[sp, 2:(1 + nages)])
  spawn_month <- ms_run$data_list$spawn_month[sp]

  sel <- rep(0, nages)
  fleet_indices <- which(ms_run$data_list$fleet_control$Fleet_species == sp)

  if (length(fleet_indices) > 0) {
    nyrs <- ms_run$data_list$nyrs
    for (flt in fleet_indices) {
      sel <- sel + ms_run$quantities$sel[flt, 1, 1:nages, nyrs]
    }
    sel <- sel / length(fleet_indices)
  } else {
    sel <- rep(1, nages)
  }

  if (max(sel) > 0) {
    sel <- sel / max(sel)
  }

  alpha <- exp(ms_run$estimated_params$rec_pars[sp, 1])
  beta <- exp(ms_run$estimated_params$rec_pars[sp, 2])

  SPR0 <- calc_SPR(0, M, wt, mat, sex_ratio, sel, spawn_month)
  rceattle_SPR0 <- ms_run$quantities$SPR0[sp]

  SPR_vec <- rep(0, length(F_grid))
  YPR_vec <- rep(0, length(F_grid))
  SSB_eq <- rep(0, length(F_grid))
  R_eq <- rep(0, length(F_grid))
  Yield_eq <- rep(0, length(F_grid))

  for (i in seq_along(F_grid)) {
    F_val <- F_grid[i]

    SPR_vec[i] <- calc_SPR(F_val, M, wt, mat, sex_ratio, sel, spawn_month)
    YPR_vec[i] <- calc_YPR(F_val, M, wt, sel)

    if (alpha * SPR_vec[i] > 1.0) {
      SSB_eq[i] <- (alpha * SPR_vec[i] - 1.0) / beta
    } else {
      SSB_eq[i] <- 0
    }

    if (SSB_eq[i] > 0) {
      R_eq[i] <- (alpha * SSB_eq[i]) / (1.0 + beta * SSB_eq[i])
    } else {
      R_eq[i] <- 0
    }

    Yield_eq[i] <- R_eq[i] * YPR_vec[i]
  }

  max_idx <- which.max(Yield_eq)

  if (length(max_idx) > 0 && Yield_eq[max_idx] > 0) {
    results$FMSY[sp] <- F_grid[max_idx]
    results$BMSY[sp] <- SSB_eq[max_idx]
    results$MSY[sp] <- Yield_eq[max_idx]
  }

  png(paste0("yield_curve_", species_names[sp], ".png"), width = 800, height = 600)
  plot(F_grid, Yield_eq, type = "l", lwd = 2,
       xlab = "Fishing mortality (F)",
       ylab = "Equilibrium yield (mt)",
       main = paste(species_names[sp], "equilibrium yield curve"))
  if (!is.na(results$FMSY[sp])) {
    abline(v = results$FMSY[sp], col = "red", lty = 2, lwd = 2)
  }
  dev.off()
}

print(results, row.names = FALSE)
write.csv(results, "FMSY_BMSY_results.csv", row.names = FALSE)
