## ── ex_warfarin_dcm.R ─────────────────────────────────────────────────────
##
## Fit a Deep Compartment Model (DCM) — Janssen et al. 2022 (CPT:PSP,
## DOI 10.1002/psp4.12808). The analytical covariate model
##
##     CL = TVCL * (WT / 70)^THETA_WT * (CRCL / 100)^THETA_CRCL * exp(ETA_CL)
##
## is replaced by a small neural network whose inputs are subject
## covariates (WT, CRCL) and whose outputs ARE the typical PK values
## (CL, V1, Q, V2, KA). Lognormal IIV is composed on top of the NN
## outputs in the standard mu-ref form (`TYPICAL_PK.CL * exp(ETA_CL)`).
##
## Requires ferx-r built with `--features nn` (which forwards to
## ferx-core's `nn` feature). To enable in a local checkout, edit
## `src/rust/Cargo.toml` and add `"nn"` to the `default` features list,
## then reinstall the package.
##
## Compare with inst/examples/ex3_two_cmt_oral_cov.R — same dataset,
## same eta/omega/sigma structure, the only difference is the covariate
## model. Run side-by-side to see the DCM output rendering vs the
## classical analytical fit.

library(ferx)

## ── 1. Load the bundled DCM example ────────────────────────────────────────
ex <- ferx_example("warfarin_dcm")
stopifnot(file.exists(ex$model), file.exists(ex$data))

cat("Model:   ", ex$model, "\n")
cat("Dataset: ", ex$data, "  (",
    length(readLines(ex$data)) - 1L, " rows)\n", sep = "")

## ── 2. Fit ─────────────────────────────────────────────────────────────────
## The bundled model declares method=focei with maxiter=200 inside
## [fit_options]. We override maxiter here for a quick demonstration;
## bump back up for a real run.
fit <- ferx_fit(
  model      = ex$model,
  data       = ex$data,
  method     = "focei"
)

## ── 3. Print the fit ───────────────────────────────────────────────────────
## The Option E rendering kicks in here: instead of 141 rows of
## `W_TYPICAL_PK_*` / `B_TYPICAL_PK_*` thetas, you'll see a compact
## --- NEURAL NETWORKS --- block summarising each NN block.
print(fit)

## ── 4. Programmatic access to the NN metadata ─────────────────────────────
## `fit$neural_networks` is a named list, one entry per [covariate_nn]
## block. Empty when ferx-r is built without `--features nn`.
if (length(fit$neural_networks) > 0) {
  nn <- fit$neural_networks[[1]]
  cat("\n--- Programmatic NN inspection ---\n")
  cat("NN name:           ", nn$name, "\n")
  cat("Layer shape:       [", paste(nn$shape, collapse = ", "), "]\n", sep = "")
  cat("Activations:        hidden=", nn$hidden_activation,
      "  output=", nn$output_activation, "\n", sep = "")
  cat("Weight count:      ", nn$n_weights, "\n")
  cat("Theta-vector slot: ", nn$weights_offset, "..",
      nn$weights_offset + nn$n_weights - 1L,
      "  (in fit$theta)\n", sep = "")
  cat("Inputs:            [", paste(nn$input_names, collapse = ", "), "]\n", sep = "")
  cat("Outputs:           [", paste(nn$output_names, collapse = ", "), "]\n", sep = "")

  ## Trained weight values (length == n_weights).
  cat("Weight summary:    min=", sprintf("%.4f", min(nn$weights)),
      "  max=",  sprintf("%.4f", max(nn$weights)),
      "  mean=", sprintf("%.4f", mean(nn$weights)),
      "  sd=",   sprintf("%.4f", stats::sd(nn$weights)),
      "\n", sep = "")

  ## Sanity check: the same weights are also reachable via fit$theta
  ## using the offset (R uses 1-based indexing; weights_offset is 0-based
  ## as exposed from Rust).
  weights_from_theta <- fit$theta[
    (nn$weights_offset + 1L):(nn$weights_offset + nn$n_weights)
  ]
  stopifnot(isTRUE(all.equal(unname(weights_from_theta), nn$weights)))
  cat("Round-trip:        fit$theta slice matches fit$neural_networks weights ✓\n")
} else {
  cat("\n",
      "NOTE: fit$neural_networks is empty. This means ferx-r was built\n",
      "      without `--features nn`. Edit src/rust/Cargo.toml's `default`\n",
      "      feature list to add \"nn\" and reinstall.\n",
      sep = "")
}

## ── 5. Are the covariates useful? Interpreting the NN ─────────────────────
##
## A classical analytical model lets us read covariate effects directly:
## a significant `THETA_WT` with %RSE under, say, 30% says "WT matters".
## A neural network is harder to read — there's no single weight to test.
## The first-order question for a DCM is therefore not "what's the weight"
## but: **is the NN actually using the inputs, or did the optimizer find a
## constant function that fits via etas alone?**
##
## Two interpretable lenses for that:
##
##   (a) Spread of per-subject typical values. If the NN ignores its
##       inputs, every subject gets the same typical CL/V/etc, and all
##       cross-subject variation is absorbed by etas. A wide range of
##       individual estimates *that correlates with the covariates*
##       indicates the NN is doing real work.
##
##   (b) Correlation between input covariates and individual PK
##       estimates. If `cor(WT, individual$CL)` is small (|r| < 0.15)
##       the NN is not propagating WT into CL. If it's substantial
##       (|r| > 0.30 ish), WT is materially shaping CL — which is what
##       you want a "deep covariate model" to do.
##
## Note: `fit$individual_estimates$CL` is at the *subject's* eta, so
## variation across subjects bundles together (i) the NN's covariate
## response and (ii) random IIV. For the correlation lens, that's fine —
## (ii) is noise around (i) and won't manufacture a covariate effect that
## isn't there.

if (length(fit$neural_networks) > 0) {
  ## Pull per-subject covariates from the source data.
  raw    <- read.csv(ex$data, stringsAsFactors = FALSE)
  cov_df <- aggregate(
    raw[, c("WT", "CRCL")],
    by = list(ID = raw$ID),
    FUN = function(x) suppressWarnings(as.numeric(x))[1]   # first non-NA per subject
  )
  ## Merge with per-subject individual estimates.
  est <- fit$individual_estimates
  est$ID <- as.character(est$ID)
  cov_df$ID <- as.character(cov_df$ID)
  joined <- merge(cov_df, est, by = "ID")

  cat("\n--- NN covariate-importance heuristics ---\n")

  ## (a) Spread of individual estimates.
  cat("\n  (a) Per-subject individual-estimate spread\n")
  cat("      (wide spread = the NN's output drives meaningful between-subject differences)\n")
  for (pk in c("CL", "V1", "Q", "V2", "KA")) {
    if (pk %in% names(joined)) {
      v <- joined[[pk]]
      cat(sprintf("      %-4s: min=%.4f  max=%.4f  CV%% across subjects = %.1f%%\n",
                  pk, min(v), max(v), sd(v) / mean(v) * 100))
    }
  }

  ## (b) Correlation lens: do the inputs explain the spread?
  cat("\n  (b) Pearson correlations: covariate vs individual PK estimate\n")
  cat("      (|r| > ~0.30 indicates the NN is propagating that covariate;\n",
      "       |r| < ~0.15 means the NN is largely ignoring it)\n", sep = "")
  cat(sprintf("      %-12s %10s %10s\n", "PK param", "r(WT)", "r(CRCL)"))
  cat(sprintf("      %s\n", strrep("-", 38)))
  for (pk in c("CL", "V1", "Q", "V2", "KA")) {
    if (pk %in% names(joined)) {
      r_wt   <- suppressWarnings(cor(joined[[pk]], joined$WT))
      r_crcl <- suppressWarnings(cor(joined[[pk]], joined$CRCL))
      cat(sprintf("      %-12s %10.3f %10.3f\n", pk, r_wt, r_crcl))
    }
  }

  ## Optional: a quick partial-dependence-style table by predicting the
  ## NN forward at each subject's actual covariates, with eta=0. We
  ## reconstruct the typical value from the individual estimate by
  ## dividing out exp(eta_i) — works because the model is mu-ref
  ## lognormal (`CL = TYPICAL_PK.CL * exp(ETA_CL)`).
  if (!is.null(fit$ebe_etas) && "ETA_CL" %in% names(fit$ebe_etas)) {
    ebe <- fit$ebe_etas
    ebe$ID <- as.character(ebe$ID)
    joined2 <- merge(joined, ebe, by = "ID", suffixes = c("", ".eta"))
    joined2$CL_tv <- joined2$CL / exp(joined2$ETA_CL)
    ord <- order(joined2$WT)
    cat("\n  Typical CL across the WT range (eta=0, low→high WT):\n")
    cat(sprintf("      WT %.0f -> CL_tv %.4f\n",
                joined2$WT[ord[1]], joined2$CL_tv[ord[1]]))
    n  <- nrow(joined2)
    cat(sprintf("      WT %.0f -> CL_tv %.4f\n",
                joined2$WT[ord[n %/% 2]], joined2$CL_tv[ord[n %/% 2]]))
    cat(sprintf("      WT %.0f -> CL_tv %.4f\n",
                joined2$WT[ord[n]], joined2$CL_tv[ord[n]]))
    cat("      (a monotonic ramp = NN learned a sensible body-weight effect;\n",
        "       flat values = NN ignored WT)\n", sep = "")
  }
}

## ── 6. Compare to a regular NLME fit on the same data ─────────────────────
##
## `examples/models/two_cpt_oral_cov.ferx` is the analytical equivalent of
## `warfarin_dcm.ferx`: same data (`data/two_cpt_oral_cov.csv`), same eta /
## omega / sigma structure, same `two_cpt_oral(...)` PK. The only
## difference is the covariate model — analytical power functions
## (`(WT/70)^THETA_WT * (CRCL/100)^THETA_CRCL`) vs the MLP.
##
## We fit it side-by-side and look at three things:
##   * OFV — same dataset, so directly comparable (smaller is better).
##   * AIC / BIC — penalise the DCM's much-higher parameter count.
##     The NN ships ~140 parameters, the analytical model ~13.
##   * Predictions — if the DCM's IPRED column tracks the analytical
##     IPRED, the NN didn't learn anything weird.
##
## For reference, `examples/models/warfarin.ferx` is the textbook
## one-compartment warfarin model with no covariates and a *different*
## dataset (warfarin.csv has no WT/CRCL). It's not data-compatible here
## but is useful as a no-covariate baseline if you load its dataset
## separately.

cat("\n\n========== Comparison: DCM vs analytical NLME ==========\n\n")
ex_an <- ferx_example("two_cpt_oral_cov")
cat("Fitting analytical comparison model (", ex_an$model, ")...\n", sep = "")

fit_an <- ferx_fit(
  model      = ex_an$model,
  data       = ex_an$data,
  method     = "focei",
  covariance = FALSE
)

cmp <- data.frame(
  model    = c("DCM (NN typical values)", "Analytical (power-function covariates)"),
  ofv      = c(fit$ofv,      fit_an$ofv),
  aic      = c(fit$aic,      fit_an$aic),
  bic      = c(fit$bic,      fit_an$bic),
  n_params = c(fit$n_parameters, fit_an$n_parameters),
  iters    = c(fit$n_iterations, fit_an$n_iterations)
)
print(cmp, row.names = FALSE)

cat("\nInterpretation:\n")
cat("  * If OFV is similar, the NN didn't extract more information from\n")
cat("    the covariates than the analytical power model — the simpler\n")
cat("    model is preferable.\n")
cat("  * If the DCM's OFV is substantially lower (>5–10 OFV units per\n")
cat("    extra parameter), the NN captured a non-linear covariate effect\n")
cat("    the analytical form missed. Check whether the AIC/BIC penalty\n")
cat("    keeps the DCM ahead.\n")
cat("  * For warfarin-style PK these two should be close — body weight\n")
cat("    and renal function on CL/V are almost linear-in-log-space, so\n")
cat("    the power model nails them. DCM tends to win on more complex\n")
cat("    drugs (e.g. monoclonal antibodies; see Janssen 2022 §3 case\n")
cat("    studies on factor VIII and edoxaban).\n")

## ── 7. Predictions look reasonable ────────────────────────────────────────
## Final sanity check: per-subject IPRED ranges from both models.
cat("\n--- Per-subject IPRED range ---\n")
cat(sprintf("  DCM:        [%.3f, %.3f]\n",
            min(fit$sdtab$IPRED, na.rm = TRUE),
            max(fit$sdtab$IPRED, na.rm = TRUE)))
cat(sprintf("  Analytical: [%.3f, %.3f]\n",
            min(fit_an$sdtab$IPRED, na.rm = TRUE),
            max(fit_an$sdtab$IPRED, na.rm = TRUE)))

## For a real DCM workflow you'd now bump `maxiter` to 500–1000 and
## generate a VPC with ferx_simulate() — see inst/examples/ex3_two_cmt_oral_cov.R
## for the analytical-side recipe; the DCM is a drop-in replacement.
