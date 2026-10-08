# ferx_sir() / ferx_covariance() check the model file against the fit's
# `model_hash` before the from-fit prologue (#492).
#
# Both re-read the model file the fit recorded. An edit that asks for
# data-derived bindings - a theta level block, a symbolic `[covariate_model]`
# centre - used to be refused by the prologue as a fit that lacks those
# bindings, which names the wrong cause: the fit is fine, the file changed. The
# glue now reads the file through ferx-core's `ModelSource::read_verified`,
# which compares the hash before it parses, so core's own `model hash mismatch`
# refusal comes first. The fixture is `two_cpt_oral_cov.csv` with
# `CL ~ WT power(center = 70)`; the fits stop after two outer iterations, as
# the assertions are about which refusal fires, not about the estimates.

hf_model_text <- "
[parameters]
  theta TVCL(4.0, 0.1, 100.0)
  theta TVV1(40.0, 1.0, 500.0)
  theta TVQ(8.0, 0.1, 100.0)
  theta TVV2(80.0, 1.0, 500.0)
  theta TVKA(1.0, 0.01, 10.0)
  omega ETA_CL ~ 0.15
  omega ETA_V1 ~ 0.15
  sigma PROP_ERR ~ 0.04 (sd)

[individual_parameters]
  CL = TVCL * exp(ETA_CL)
  V1 = TVV1 * exp(ETA_V1)
  Q  = TVQ
  V2 = TVV2
  KA = TVKA

[covariates]
  WT continuous

[covariate_model]
  CL ~ WT power(center = 70) => THETA_CL_WT(0.6, 0.01, 5.0)

[structural_model]
  pk two_cpt_oral(cl=CL, v1=V1, q=Q, v2=V2, ka=KA)

[error_model]
  DV ~ proportional(PROP_ERR)

[fit_options]
  method   = focei
  maxiter  = 2
  covariance = false
"

# The two edits of #492's table. Each asks for bindings the fit was not made
# with, so without the hash check first each is refused as a fit lacking them.
hf_edits <- list(
  level_block = c("theta TVKA(1.0, 0.01, 10.0)",
                  "theta TVKA[ID, contrast = none](1.0, 0.01, 10.0)"),
  symbolic_centre = c("center = 70", "center = median")
)

# A fit of the unedited model, in its own directory so the edit cannot reach
# another test's file. The covariance matrix is a placeholder: SIR needs one
# to run at all, and no assertion here reaches the resampling.
hf_fit <- function() {
  dir <- tempfile("ferx_hash_first_")
  dir.create(dir)
  model <- file.path(dir, "m.ferx")
  data <- file.path(dir, "d.csv")
  writeLines(hf_model_text, model)
  file.copy(ferx_example("two_cpt_oral_cov")$data, data)
  fit <- ferx_fit(model, data, verbose = FALSE)
  n <- length(fit$theta) + nrow(fit$omega) + length(fit$sigma)
  fit$cov_matrix <- diag(1e-4, n)
  list(fit = fit, model = model)
}

hf_edit <- function(model, edit) {
  text <- readLines(model)
  stopifnot(sum(grepl(edit[1L], text, fixed = TRUE)) == 1L)
  writeLines(sub(edit[1L], edit[2L], text, fixed = TRUE), model)
}

hf_msg <- function(expr) {
  tryCatch({
    expr
    NA_character_
  }, error = function(e) conditionMessage(e))
}

hf_calls <- list(
  ferx_sir = function(fit) {
    ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L, sir_seed = 1L)
  },
  ferx_covariance = function(fit) ferx_covariance(fit)
)

test_that("an edited model file is refused as edited, not as a fit lacking bindings", {
  for (edit in names(hf_edits)) {
    h <- hf_fit()
    hf_edit(h$model, hf_edits[[edit]])
    for (entry in names(hf_calls)) {
      info <- paste(edit, entry)
      msg <- hf_msg(hf_calls[[entry]](h$fit))
      # Core's text, under this entry's name, naming the file.
      expect_match(msg, paste0("^", entry, ": model hash mismatch for "),
                   info = info)
      expect_match(msg, h$fit$model_path, fixed = TRUE, info = info)
      expect_match(msg, "The .ferx file has changed since the fit was produced",
                   fixed = TRUE, info = info)
      # Not either from-fit refusal the edit used to reach.
      expect_no_match(msg, "carries no", fixed = TRUE, info = info)
      # Core gives the stale-source refusal no diagnostic code (ferx-core
      # #1746), and R does not borrow one.
      expect_no_match(msg, "\\[E_[A-Z_]+\\]$", info = info)
    }
  }
})

test_that("a fit with no model_hash is not checked: the from-fit prologue still answers", {
  # The other side of the gate. An empty hash disables the check (an older
  # fit, or one whose hash failed at fit time), so the edited level block
  # reaches the prologue's own refusal.
  h <- hf_fit()
  hf_edit(h$model, hf_edits$level_block)
  fit <- h$fit
  fit$model_hash <- NA_character_
  for (entry in names(hf_calls)) {
    msg <- hf_msg(suppressWarnings(hf_calls[[entry]](fit)))
    expect_match(msg, "this fit carries no theta level bindings", fixed = TRUE,
                 info = entry)
    expect_no_match(msg, "hash mismatch", fixed = TRUE, info = entry)
  }
})
