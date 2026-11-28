# Evaluate `expr`, muffling only the false "theta 'MIXL' ... not referenced"
# lint the engine raises on a [mixture] model's mixing-probability theta
# (ferx-core #2483). The simulate wrappers fold all of a run's diagnostics into
# one R warning, so only a warning whose sole diagnostic is that lint is
# muffled: any other diagnostic, alone or alongside it, still surfaces. Drop
# this helper and its callers' wrapping once #2483 is in the pin.
muffle_mixl_lint <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("produced 1 diagnostic:\n  theta 'MIXL' is declared in [parameters] but not referenced",
              conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# A cached FOCEI fit of fixtures/mixture_override.ferx on mixture_iv.csv: a
# two-class [mixture] model with a class-2 ETA_CL override (`omega(2)`), whose
# fitted base omega and override sit far from the model file's 0.09 / 0.30.
# Shared by the from-fit [mixture] tests (ferx-r #672).
mixture_override_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- muffle_mixl_lint(ferx_fit(
        test_path("fixtures", "mixture_override.ferx"),
        test_path("fixtures", "mixture_iv.csv"),
        method = "focei", verbose = FALSE
      ))
    }
    fit
  }
})

# The mixture_override fixture as an ODE model with an [adaptive_dosing] block.
.write_mixture_adaptive_model <- function() {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL1(1.2, 0.01, 100.0)",
    "  theta TVCL2(2.5, 0.01, 100.0)",
    "  theta TVV(10.0, 0.1, 1000.0)",
    "  theta MIXL(0.0, -10.0, 10.0)",
    "  omega ETA_CL ~ 0.09",
    "  sigma EPS ~ 0.04",
    "[mixture]",
    "  nsub = 2",
    "  logit(1) = MIXL",
    "  omega(2) ETA_CL ~ 0.30",
    "[individual_parameters]",
    "  CL = if (MIXNUM == 1) TVCL1 * exp(ETA_CL) else TVCL2 * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  ode(states=[central])",
    "[odes]",
    "  d/dt(central) = -(CL / V) * central",
    "[scaling]",
    "  y = central / V",
    "[error_model]",
    "  DV ~ proportional(EPS)",
    "[adaptive_dosing]",
    "  observe     = central / V",
    "  at          = every 12 from 0 to 48",
    "  start_dose  = 100",
    "  route       = bolus(cmt=1)",
    "  dose_bounds = [0, 200]",
    "  confirm     = 1",
    "  when signal < 1 : increase 25%"
  ), path)
  path
}

# `fit` as a `.fitrx` saved before #655 has it (#675): no override values,
# descriptors or SEs (the loader refuses descriptors without values, ferx-core
# #2833), packed estimate kept. Edits the bytes, never the list.
mixture_fit_without_override_values <- function(fit) {
  fitrx_edit_fit(fit, function(w) {
    w$mixture_overrides <- NULL
    w$mixture_override_params <- NULL
    w$se_mixture_overrides <- NULL
    w
  })
}
