## ── ex_warfarin_dcm.R ─────────────────────────────────────────────────────
##
## Fit a Deep Compartment Model (DCM) — Janssen et al. 2022 (CPT:PSP,
## DOI 10.1002/psp4.12808) — and show *when* it actually helps over a
## classical analytical NLME fit.
##
## DCM replaces the analytical covariate model
##
##     CL = TVCL * (WT / 70)^THETA_WT * (CRCL / 100)^THETA_CRCL * exp(ETA_CL)
##
## with a small neural network whose inputs are subject covariates and
## whose outputs ARE the typical PK values. Lognormal IIV is composed on
## top in the standard mu-ref form (`TYPICAL_PK.CL * exp(ETA_CL)`).
##
## ──────────── Why this script generates its own data ───────────────────
##
## The bundled `two_cpt_oral_cov.csv` is simulated with a clean
## power-function effect — `CL ∝ (WT/70)^0.75`, `CL ∝ (CRCL/100)^0.50`.
## An analytical power-function model fits that perfectly, so a DCM
## would tie or lose on AIC/BIC. That's the wrong story to tell about
## DCM: the value-add is **non-linear covariate effects**, the kind a
## power function can't represent.
##
## So we use the bundled `warfarin_if.ferx` model — which has a sharp
## kink in clearance at WT = 70 (`if (WT > 70)` allometric, else flat)
## — to simulate a fresh dataset where:
##   * WT genuinely drives CL (with a non-linear kink)
##   * CRCL is present but doesn't affect anything (a "distractor")
##
## We then fit DCM and analytical models to that simulated data and
## compare. DCM should win on AIC because the analytical power model
## can only find a single compromise exponent across the kink.
##
## ──────────── Requirements ──────────────────────────────────────────────
## Requires ferx-r built with `--features nn` (which forwards to
## ferx-core's `nn` feature). Edit `src/rust/Cargo.toml` and add `"nn"`
## to the `default` features list, then reinstall the package.

library(ferx)
set.seed(42)

## ── 1. Simulate a dataset with a non-linear covariate effect ──────────────
##
## Strategy:
##   (i)   Take the bundled warfarin_if data file's covariate + dosing rows
##         (it already has WT spanning the kink and CRCL noise columns).
##   (ii)  Use `ferx_simulate()` with the warfarin_if MODEL file (which has
##         the `if (WT > 70)` kink in its [individual_parameters]) to
##         simulate fresh DV concentrations from that non-linear ground
##         truth.
##   (iii) Stitch the simulated DVs back into the NONMEM-format CSV,
##         write to a tempfile, and use that as our "real" dataset for
##         both fits.
gt <- ferx_example("warfarin_if")
cat("Ground-truth model: ", gt$model, "\n")
cat("Simulating fresh DVs from a non-linear kink-at-WT=70 model...\n")

sim <- ferx_simulate(gt$model, gt$data, n_sim = 1L, seed = 42L)
## sim columns: SIM, ID, TIME, IPRED, DV_SIM

## Stitch simulated DV back into the NONMEM CSV. Observation rows
## (EVID==0, MDV==0) get their DV replaced; dose rows stay as-is.
orig <- read.csv(gt$data, stringsAsFactors = FALSE)
obs_mask <- orig$EVID == 0 & orig$MDV == 0
key_orig <- paste(orig$ID[obs_mask], orig$TIME[obs_mask])
key_sim  <- paste(sim$ID,            sim$TIME)
orig$DV  <- as.character(orig$DV)   # original CSV uses "." for dose rows
orig$DV[obs_mask] <- sprintf("%.4f", sim$DV_SIM[match(key_orig, key_sim)])

dataset_path <- tempfile(fileext = ".csv")
write.csv(orig, dataset_path, row.names = FALSE, quote = FALSE)
cat("Simulated dataset:  ", dataset_path, "\n")
cat("                    ", sum(obs_mask), "observations across",
    length(unique(orig$ID)), "subjects\n")
cat("                    WT range:  ",
    sprintf("%.1f to %.1f", min(orig$WT), max(orig$WT)), "\n")
cat("                    CRCL range:",
    sprintf("%.1f to %.1f", min(orig$CRCL), max(orig$CRCL)),
    " (distractor — no PK effect in the ground truth)\n")

## ── 2. Fit the DCM ─────────────────────────────────────────────────────────
##
## The bundled `warfarin_dcm.ferx` declares method=focei + maxiter=200.
## We override outer_maxiter here for a quick demonstration; bump back
## up for a real fit.
ex_dcm <- ferx_example("warfarin_dcm")
cat("\n── Fitting DCM model:", ex_dcm$model, "──\n")
fit <- ferx_fit(
  model      = ex_dcm$model,
  data       = dataset_path,
  method     = "focei",
  covariance = FALSE,
  settings   = list(outer_maxiter = 100)
)

## ── 3. Print the fit (Option E rendering) ─────────────────────────────────
## NN-weight thetas are collapsed into the `--- NEURAL NETWORKS ---` block
## instead of cluttering the THETA Estimates table.
print(fit)

## ── 4. Programmatic access to the NN metadata ─────────────────────────────
if (length(fit$neural_networks) > 0) {
  nn <- fit$neural_networks[[1]]
  cat("\n--- Programmatic NN inspection ---\n")
  cat("NN name:           ", nn$name, "\n")
  cat("Layer shape:       [", paste(nn$shape, collapse = ", "), "]\n", sep = "")
  cat("Activations:        hidden=", nn$hidden_activation,
      "  output=", nn$output_activation, "\n", sep = "")
  cat("Weight count:      ", nn$n_weights, "\n")
  cat("Theta-vector slot: ", nn$weights_offset, "..",
      nn$weights_offset + nn$n_weights - 1L, "\n", sep = "")
  cat("Inputs:            [", paste(nn$input_names, collapse = ", "), "]\n", sep = "")
  cat("Outputs:           [", paste(nn$output_names, collapse = ", "), "]\n", sep = "")
} else {
  cat("\nNOTE: fit$neural_networks is empty. ferx-r was built without\n")
  cat("      `--features nn`. Edit src/rust/Cargo.toml's `default`\n")
  cat("      feature list to add \"nn\" and reinstall.\n")
}

## ── 5. Are the covariates useful? Interpreting the NN ─────────────────────
##
## We expect to find:
##   * `cor(WT, CL)` substantially > 0  — WT drove CL in the ground truth
##   * `cor(CRCL, CL)` close to 0       — CRCL is a distractor; NN should ignore
##
## Two interpretable lenses on the fitted NN:
##
##   (a) Spread of per-subject typical values. If the NN ignored its
##       inputs, every subject gets the same typical CL/V/etc and all
##       cross-subject variation is absorbed by etas — flat individual
##       estimates. A wide spread that *correlates with the covariates*
##       indicates the NN is doing real work.
##
##   (b) Correlation between input covariates and individual PK
##       estimates. The expected pattern for this ground truth: strong
##       `cor(WT, CL)`, near-zero `cor(CRCL, CL)`.

if (length(fit$neural_networks) > 0) {
  raw <- read.csv(dataset_path, stringsAsFactors = FALSE)
  cov_df <- aggregate(
    raw[, c("WT", "CRCL")],
    by = list(ID = raw$ID),
    FUN = function(x) suppressWarnings(as.numeric(x))[1]
  )
  est <- fit$individual_estimates
  est$ID <- as.character(est$ID)
  cov_df$ID <- as.character(cov_df$ID)
  joined <- merge(cov_df, est, by = "ID")

  cat("\n--- NN covariate-importance heuristics ---\n")

  cat("\n  (a) Per-subject individual-estimate spread\n")
  cat("      (wide spread = the NN's output drives meaningful between-subject differences)\n")
  for (pk in c("CL", "V1", "Q", "V2", "KA")) {
    if (pk %in% names(joined)) {
      v <- joined[[pk]]
      cat(sprintf("      %-4s: min=%.4f  max=%.4f  CV%% across subjects = %.1f%%\n",
                  pk, min(v), max(v), sd(v) / mean(v) * 100))
    }
  }

  cat("\n  (b) Pearson correlations: covariate vs individual PK estimate\n")
  cat("      Ground truth: WT drives CL with a kink at WT=70; CRCL has no effect.\n")
  cat("      Expected: |r(WT, CL)| substantially > 0; |r(CRCL, *)| close to 0.\n")
  cat(sprintf("      %-12s %10s %10s\n", "PK param", "r(WT)", "r(CRCL)"))
  cat(sprintf("      %s\n", strrep("-", 38)))
  for (pk in c("CL", "V1", "Q", "V2", "KA")) {
    if (pk %in% names(joined)) {
      r_wt   <- suppressWarnings(cor(joined[[pk]], joined$WT))
      r_crcl <- suppressWarnings(cor(joined[[pk]], joined$CRCL))
      cat(sprintf("      %-12s %10.3f %10.3f\n", pk, r_wt, r_crcl))
    }
  }

  ## Partial-dependence-style: typical CL vs WT, reconstructed from
  ## individual_estimates / exp(eta) (mu-ref lognormal model).
  if (!is.null(fit$ebe_etas) && "ETA_CL" %in% names(fit$ebe_etas)) {
    ebe <- fit$ebe_etas
    ebe$ID  <- as.character(ebe$ID)
    j2      <- merge(joined, ebe, by = "ID", suffixes = c("", ".eta"))
    j2$CL_tv <- j2$CL / exp(j2$ETA_CL)
    ord <- order(j2$WT)
    cat("\n  Typical CL as a function of WT (eta=0, low→high WT):\n")
    cat("  A non-linear ramp (especially a kink near WT=70) means the NN\n")
    cat("  recovered the simulated kink. A flat line means it ignored WT.\n")
    breaks <- round(seq(1, nrow(j2), length.out = 5))
    for (k in breaks) {
      cat(sprintf("      WT %5.1f  →  CL_tv = %.4f\n",
                  j2$WT[ord[k]], j2$CL_tv[ord[k]]))
    }
  }
}

## ── 6. Compare to a regular NLME fit on the same simulated data ───────────
##
## `two_cpt_oral_cov.ferx` is the textbook analytical NLME — same eta /
## omega / sigma, same `two_cpt_oral` PK, but covariates enter as
## power functions: `CL = TVCL * (WT/70)^THETA_WT * (CRCL/100)^THETA_CRCL`.
##
## On our kink-y simulated data the power function can't represent the
## sharp change in CL at WT=70 cleanly — it has to pick a single exponent
## that compromises across both regimes. The NN, with extra capacity, can
## represent the kink. So we expect DCM to win on OFV; the AIC/BIC
## verdict depends on whether the OFV improvement exceeds the NN's
## parameter-count penalty.
ex_an <- ferx_example("two_cpt_oral_cov")
cat("\n\n══════════ Comparison: DCM vs analytical NLME ══════════\n\n")
cat("Fitting analytical comparison model (", ex_an$model, ")...\n", sep = "")
fit_an <- ferx_fit(
  model      = ex_an$model,
  data       = dataset_path,
  method     = "focei",
  covariance = FALSE,
  settings   = list(outer_maxiter = 100)
)

cmp <- data.frame(
  model    = c("DCM (NN typical values)", "Analytical (power-function covariates)"),
  ofv      = c(fit$ofv,           fit_an$ofv),
  aic      = c(fit$aic,           fit_an$aic),
  bic      = c(fit$bic,           fit_an$bic),
  n_params = c(fit$n_parameters,  fit_an$n_parameters),
  iters    = c(fit$n_iterations,  fit_an$n_iterations)
)
print(cmp, row.names = FALSE)

ofv_delta <- fit_an$ofv - fit$ofv
aic_delta <- fit_an$aic - fit$aic
cat("\nΔOFV (analytical − DCM): ", sprintf("%+.2f", ofv_delta),
    "  (positive ⇒ DCM fits better)\n", sep = "")
cat("ΔAIC (analytical − DCM): ", sprintf("%+.2f", aic_delta),
    "  (positive ⇒ DCM wins overall after the parameter penalty)\n", sep = "")

cat("\nInterpretation:\n")
if (aic_delta > 0) {
  cat("  ✓ DCM wins on AIC. The NN captured the non-linear WT effect that\n")
  cat("    the analytical power model couldn't represent. The extra weight\n")
  cat("    parameters earned their keep.\n")
} else if (ofv_delta > 0) {
  cat("  ~ DCM has lower OFV but loses on AIC. The NN fit better but its\n")
  cat("    extra parameters aren't justified by the OFV gain. The simpler\n")
  cat("    analytical model is preferable here.\n")
} else {
  cat("  ✗ DCM didn't improve on the analytical model. Either the analytical\n")
  cat("    form already captured the covariate effect (so the kink wasn't\n")
  cat("    actually a problem for it) or the DCM hasn't converged — bump\n")
  cat("    outer_maxiter and retry.\n")
}

## ── 7. Per-subject IPRED sanity check ────────────────────────────────────
cat("\n--- Per-subject IPRED range (both models) ---\n")
cat(sprintf("  DCM:        [%.3f, %.3f]\n",
            min(fit$sdtab$IPRED, na.rm = TRUE),
            max(fit$sdtab$IPRED, na.rm = TRUE)))
cat(sprintf("  Analytical: [%.3f, %.3f]\n",
            min(fit_an$sdtab$IPRED, na.rm = TRUE),
            max(fit_an$sdtab$IPRED, na.rm = TRUE)))

## ── Wrap-up ───────────────────────────────────────────────────────────────
##
## For a production DCM workflow you'd now:
##   * Bump outer_maxiter to 500–1000 for actual convergence.
##   * Generate a VPC with ferx_simulate() (see ex3_two_cmt_oral_cov.R).
##   * If the AIC verdict is "analytical wins", that's a real result, not
##     a failure — DCM only helps when covariates are genuinely non-linear.
##   * If DCM wins, inspect the per-output correlations (§5b above) to
##     identify which covariates are pulling weight. Drop "distractor"
##     covariates from `inputs = [...]` to reduce parameter count.
##
## Reference: Janssen A. et al. (2022). Deep compartment models: A deep
## learning approach for the reliable prediction of time-series data in
## pharmacokinetic modeling. CPT Pharmacometrics Syst Pharmacol 11:934-945.
## DOI 10.1002/psp4.12808.
