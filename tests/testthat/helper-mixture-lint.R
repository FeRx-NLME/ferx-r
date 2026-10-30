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
