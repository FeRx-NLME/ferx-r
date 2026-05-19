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
## with a small neural network whose inputs are subject covariates. The
## bundled `warfarin_dcm.ferx` uses the **baseline × NN-modulator**
## factoring from Janssen 2022 — baseline thetas (TVCL, TVV1, …) carry the
## absolute PK scale, and the NN learns a multiplicative covariate-driven
## modulator on top:
##
##     CL = TVCL * TYPICAL_PK.CL * exp(ETA_CL)
##           └── scale ─┘   └── NN modulator, ≈ 1.0 at init ──┘
##
## Trying to make the NN produce *absolute* PK values directly from
## Glorot-init weights fails — softplus(0) ≈ 0.69 for every output and
## the optimizer can't span the 1-100× range across CL/V/Q/V2/KA. See
## the `warfarin_dcm.ferx` block comment for the full rationale.
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
## We override maxiter for a quick demonstration; bump back up for a
## real fit.
##
## We also override `optimizer = "lbfgs"`. The default outer optimizer
## (SLSQP) currently silently terminates at iter 1 with theta unchanged
## from init when mu-referencing is active — open issue in ferx-core.
## L-BFGS, BOBYQA, and BFGS all work; we pick L-BFGS here for its
## gradient-driven convergence speed.
ex_dcm <- ferx_example("warfarin_dcm")
cat("\n── Fitting DCM model:", ex_dcm$model, "──\n")
fit <- ferx_fit(
  model      = ex_dcm$model,
  data       = dataset_path,
  method     = "focei",
  covariance = FALSE,
  settings   = list(maxiter = 100, optimizer = "lbfgs")
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
  settings   = list(maxiter = 100, optimizer = "lbfgs")
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
  cat("    maxiter and retry.\n")
}

## ── 7. Per-subject IPRED sanity check ────────────────────────────────────
cat("\n--- Per-subject IPRED range (both models) ---\n")
cat(sprintf("  DCM:        [%.3f, %.3f]\n",
            min(fit$sdtab$IPRED, na.rm = TRUE),
            max(fit$sdtab$IPRED, na.rm = TRUE)))
cat(sprintf("  Analytical: [%.3f, %.3f]\n",
            min(fit_an$sdtab$IPRED, na.rm = TRUE),
            max(fit_an$sdtab$IPRED, na.rm = TRUE)))

## ── 8. Visual diagnostics ──────────────────────────────────────────────────
##
## Standard PK goodness-of-fit (GOF) plots plus a covariate-effect plot
## that's specific to DCM: typical CL as a function of WT, for both
## models. The kink-at-WT=70 ground truth should appear as a bend in the
## DCM curve and a smooth power-law line for the analytical model — the
## visual signature of "the analytical form can't capture this".
##
## Base R only — no `ggplot2` / `vpc` dependency. Output goes to a PDF
## under tempdir; we print the path. Open with `open <path>` on macOS,
## `xdg-open` on Linux, or just attach to your R session interactively
## and re-run without `pdf()`/`dev.off()` to plot to a device.

cat("\n── Generating diagnostic plots ──\n")

## Helper: build typical-CL-by-WT data frame for one fit (lognormal mu-ref:
## CL = CL_tv * exp(ETA_CL), so CL_tv = CL_indiv / exp(ETA_CL)).
build_tv_cl <- function(this_fit, cov_df) {
  ebe <- this_fit$ebe_etas
  ebe$ID <- as.character(ebe$ID)
  est <- this_fit$individual_estimates
  est$ID <- as.character(est$ID)
  df <- merge(est, ebe, by = "ID", suffixes = c("", ".eta"))
  df <- merge(df, cov_df, by = "ID")
  df$CL_tv <- df$CL / exp(df$ETA_CL)
  df[order(df$WT), c("ID", "WT", "CRCL", "CL", "ETA_CL", "CL_tv")]
}

## Subject-level covariate frame (one row per subject) — re-pulled from
## the same simulated CSV so it stays consistent regardless of fit order.
raw <- read.csv(dataset_path, stringsAsFactors = FALSE)
cov_df <- aggregate(
  raw[, c("WT", "CRCL")],
  by = list(ID = raw$ID),
  FUN = function(x) suppressWarnings(as.numeric(x))[1]
)
cov_df$ID <- as.character(cov_df$ID)
tv_dcm <- build_tv_cl(fit,    cov_df)
tv_an  <- build_tv_cl(fit_an, cov_df)

## ── Page 1: GOF + typical-CL-vs-WT ──────────────────────────────────────
pdf_path <- tempfile(fileext = ".pdf")
pdf(pdf_path, width = 11, height = 8.5)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

## DV vs IPRED — both models, log-log so low concentrations are visible.
plot_dv_pred <- function(sdtab, title, col_pred = "IPRED") {
  xy <- range(c(sdtab$DV, sdtab[[col_pred]]), na.rm = TRUE)
  xy[1] <- max(xy[1], 0.01)  # guard against zeros for log scale
  plot(sdtab[[col_pred]], sdtab$DV, pch = 20, col = rgb(0, 0, 0, 0.35),
       log = "xy", xlim = xy, ylim = xy,
       main = title, xlab = col_pred, ylab = "DV")
  abline(0, 1, col = "red", lwd = 2)
}

plot_dv_pred(fit$sdtab,    "DCM:  DV vs IPRED",        col_pred = "IPRED")
plot_dv_pred(fit_an$sdtab, "Analytical:  DV vs IPRED", col_pred = "IPRED")

## Typical CL as a function of WT, both models overlaid.
## With the kink ground truth (CL flat below WT=70, allometric above),
## the analytical line should be smooth-monotone (it has only one
## exponent to play with), while the DCM curve should bend near WT=70
## if the NN captured the kink.
xrng <- range(c(tv_dcm$WT, tv_an$WT))
yrng <- range(c(tv_dcm$CL_tv, tv_an$CL_tv))
plot(tv_dcm$WT, tv_dcm$CL_tv, type = "b", pch = 20, col = "steelblue",
     xlim = xrng, ylim = yrng,
     main = "Typical CL vs WT", xlab = "WT", ylab = "Typical CL (eta = 0)")
lines(tv_an$WT, tv_an$CL_tv, type = "b", pch = 20, col = "firebrick")
abline(v = 70, lty = 2, col = "gray60")
text(70, yrng[1] + diff(yrng) * 0.04, " ground-truth kink (WT = 70)",
     pos = 4, col = "gray40", cex = 0.9)
legend("topleft",
       legend = c("DCM (NN)", "Analytical (power)"),
       col = c("steelblue", "firebrick"),
       lty = 1, pch = 20, bty = "n")

## DV vs PRED (population predictions — eta=0). Captures whether the
## structural model + covariate model alone are adequate, independent of
## individual eta estimates.
plot_dv_pred(fit$sdtab,    "DCM:  DV vs PRED",        col_pred = "PRED")
plot_dv_pred(fit_an$sdtab, "Analytical:  DV vs PRED", col_pred = "PRED")

## CWRES vs TIME — symmetric clouds centred on 0 ⇒ residuals are well
## behaved. Drift, fanning, or persistent sign indicate misspecification.
plot_cwres <- function(sdtab, title) {
  plot(sdtab$TIME, sdtab$CWRES, pch = 20, col = rgb(0, 0, 0, 0.35),
       main = title, xlab = "TIME", ylab = "CWRES",
       ylim = range(c(sdtab$CWRES, c(-3, 3)), na.rm = TRUE))
  abline(h = 0, col = "red", lwd = 2)
  abline(h = c(-2, 2), lty = 2, col = "gray60")
}
plot_cwres(fit$sdtab,    "DCM:  CWRES vs TIME")
plot_cwres(fit_an$sdtab, "Analytical:  CWRES vs TIME")

## ── Page 2: Individual concentration profiles ──────────────────────────
## Pick 4 representative subjects: 2 with WT < 70 (below the kink, where
## the ground truth has no covariate effect on CL) and 2 with WT > 70
## (above the kink, where allometric scaling kicks in). For each subject,
## plot observed DV + IPRED from both models over TIME.
ids_below <- head(cov_df$ID[cov_df$WT < 70 & cov_df$WT > 55], 2)
ids_above <- head(cov_df$ID[cov_df$WT > 70], 2)
selected_ids <- c(ids_below, ids_above)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (id in selected_ids) {
  wt_i <- cov_df$WT[cov_df$ID == id]
  ## sdtab.ID is numeric; cast for the comparison.
  rows_dcm <- fit$sdtab[fit$sdtab$ID == as.numeric(id), ]
  rows_an  <- fit_an$sdtab[fit_an$sdtab$ID == as.numeric(id), ]
  ord_dcm  <- order(rows_dcm$TIME)
  ord_an   <- order(rows_an$TIME)

  yrng <- range(c(rows_dcm$DV, rows_dcm$IPRED, rows_an$IPRED), na.rm = TRUE)
  plot(rows_dcm$TIME[ord_dcm], rows_dcm$DV[ord_dcm],
       pch = 19, col = "black", ylim = yrng,
       main = sprintf("Subject %s  (WT=%.1f)", id, wt_i),
       xlab = "TIME", ylab = "Concentration")
  lines(rows_dcm$TIME[ord_dcm], rows_dcm$IPRED[ord_dcm],
        col = "steelblue", lwd = 2)
  lines(rows_an$TIME[ord_an],  rows_an$IPRED[ord_an],
        col = "firebrick",  lwd = 2, lty = 2)
  if (id == selected_ids[1]) {
    legend("topright",
           legend = c("Observed", "DCM IPRED", "Analytical IPRED"),
           col   = c("black", "steelblue", "firebrick"),
           pch   = c(19, NA, NA),
           lty   = c(NA, 1, 2),
           bty = "n", cex = 0.85)
  }
}

dev.off()
cat("Diagnostic plots written to: ", pdf_path, "\n", sep = "")
cat("Open with:                   open ", shQuote(pdf_path), "\n", sep = "")

## ── 9. Visual Predictive Check (optional, slower) ──────────────────────────
##
## Set RUN_VPC <- TRUE to also produce a simple 2-panel VPC: the
## 5/50/95th percentiles of simulated concentrations from each fit,
## overlaid on observed DVs. Costs ~10-30 s extra wall time per model.
RUN_VPC <- FALSE

if (RUN_VPC) {
  cat("\n── Running visual predictive checks (200 replicates per model) ──\n")
  vpc_n  <- 200L
  sim_dcm <- ferx_simulate(ex_dcm$model, dataset_path, n_sim = vpc_n,
                           seed = 1L, fit = fit)
  sim_an  <- ferx_simulate(ex_an$model,  dataset_path, n_sim = vpc_n,
                           seed = 1L, fit = fit_an)

  ## Bin times into ~10 bins and compute observed + simulated percentiles
  ## per bin. Simple binned VPC — not as fancy as the vpc package's
  ## prediction-corrected version but enough to spot bias.
  binned_pctiles <- function(times, values, n_bins = 10) {
    breaks <- quantile(times, probs = seq(0, 1, length.out = n_bins + 1),
                       na.rm = TRUE, names = FALSE)
    breaks[1] <- breaks[1] - 1e-9
    bin <- cut(times, breaks = breaks, include.lowest = TRUE)
    centres <- tapply(times, bin, mean, na.rm = TRUE)
    q05 <- tapply(values, bin, quantile, probs = 0.05, na.rm = TRUE)
    q50 <- tapply(values, bin, quantile, probs = 0.50, na.rm = TRUE)
    q95 <- tapply(values, bin, quantile, probs = 0.95, na.rm = TRUE)
    data.frame(t = as.numeric(centres),
               q05 = as.numeric(q05),
               q50 = as.numeric(q50),
               q95 = as.numeric(q95))
  }

  obs_q <- binned_pctiles(fit$sdtab$TIME, fit$sdtab$DV)
  sim_q_dcm <- binned_pctiles(sim_dcm$TIME, sim_dcm$DV_SIM)
  sim_q_an  <- binned_pctiles(sim_an$TIME,  sim_an$DV_SIM)

  vpc_pdf <- tempfile(fileext = ".pdf")
  pdf(vpc_pdf, width = 10, height = 5)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  yrng <- range(c(obs_q$q05, obs_q$q95, sim_q_dcm$q05, sim_q_dcm$q95,
                  sim_q_an$q05, sim_q_an$q95), na.rm = TRUE)

  for (panel in list(list(q = sim_q_dcm, title = "VPC — DCM",        col = "steelblue"),
                     list(q = sim_q_an,  title = "VPC — Analytical", col = "firebrick"))) {
    sim_q <- panel$q
    plot(obs_q$t, obs_q$q50, type = "n", ylim = yrng,
         main = panel$title, xlab = "TIME",
         ylab = "DV (5/50/95th percentile)")
    polygon(c(sim_q$t, rev(sim_q$t)),
            c(sim_q$q05, rev(sim_q$q95)),
            col = adjustcolor(panel$col, alpha.f = 0.25), border = NA)
    lines(sim_q$t, sim_q$q50, col = panel$col, lwd = 2)
    lines(obs_q$t, obs_q$q05, col = "black", lty = 3)
    lines(obs_q$t, obs_q$q50, col = "black", lwd = 2)
    lines(obs_q$t, obs_q$q95, col = "black", lty = 3)
    legend("topright",
           legend = c("Obs 5/50/95%", paste("Sim 5–95% (", panel$title, ")", sep = ""),
                      "Sim median"),
           col = c("black", panel$col, panel$col),
           lty = c(1, NA, 1), lwd = c(2, NA, 2),
           fill = c(NA, adjustcolor(panel$col, alpha.f = 0.25), NA),
           border = c(NA, panel$col, NA),
           bty = "n", cex = 0.75)
  }
  dev.off()
  cat("VPC plots written to:        ", vpc_pdf, "\n", sep = "")
}

## ── Wrap-up ───────────────────────────────────────────────────────────────
##
## For a production DCM workflow you'd now:
##   * Bump maxiter to 500–1000 for actual convergence.
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
