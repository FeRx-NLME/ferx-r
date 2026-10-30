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
