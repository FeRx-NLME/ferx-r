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
  method     = "focei",
  covariance = FALSE,        # disable SE estimation for speed
  settings   = list(outer_maxiter = 50)  # quick demo; raise for production
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

## ── 5. Predictions look reasonable ────────────────────────────────────────
## Even at 50 outer iterations with a Glorot-initialised NN, FOCEI should
## bring per-subject predictions into the right order of magnitude. For a
## real DCM workflow, fit to convergence (maxiter ~500-1000) and then
## generate a VPC with ferx_simulate() as in ex3_two_cmt_oral_cov.R.
cat("\nOFV:", sprintf("%.2f", fit$ofv), "\n")
cat("Per-subject ipred range:",
    sprintf("%.4f", min(fit$sdtab$IPRED, na.rm = TRUE)),
    "to",
    sprintf("%.4f", max(fit$sdtab$IPRED, na.rm = TRUE)),
    "\n")
