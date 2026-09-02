# ferx-core #993 / #1004: a dose attribute the engine applies at the dose event
# (`F`, `LAGTIME`/`ALAG`, and their compartment-indexed forms) may not also be
# read on the prediction path. #993 covered ODE models; #1004 extended it to
# analytical `pk ...` models, where the `pk(..., f=F)` *mapping* - not the
# parameter's name - is what binds it.
#
# These tests exist because the R package documents that rejection in NEWS.md
# but had nothing asserting it from R: a scope regression on either half would
# have landed here green.

.dad_write <- function(lines, path = tempfile(fileext = ".ferx")) {
  writeLines(lines, path)
  path
}

# Analytical one-compartment oral model. `pk_args` and `readout` are the two
# halves the check pairs up, so each test varies exactly one of them.
.dad_model <- function(pk_args = "cl=CL, v=V, ka=KA, f=F",
                       readout = c("[scaling]", "  obs_scale = V / F")) {
  c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  theta TVV(8.1, 0.1, 500.0)",
    "  theta TVKA(1.0, 0.01, 50.0)",
    "  theta TVF(0.9, 0.1, 1.0)",
    "",
    "  omega ETA_CL ~ 0.07",
    "",
    "  sigma PROP_ERR ~ 0.01 (sd)",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "  KA = TVKA",
    "  F  = TVF",
    "",
    "[structural_model]",
    paste0("  pk one_cpt_oral(", pk_args, ")"),
    "",
    readout,
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  )
}

# `ferx_model_validate()` prints its report and returns the structured result
# invisibly, so swallow the console output and hand back the value.
.dad_validate <- function(path) {
  out <- NULL
  suppressWarnings(capture.output(out <- ferx_model_validate(path)))
  out
}

.dad_codes <- function(path) .dad_validate(path)$diagnostics$code

test_that("an analytical model that maps `f=F` and reads `F` back is rejected", {
  path <- .dad_write(.dad_model())
  on.exit(unlink(path))

  out <- .dad_validate(path)
  expect_false(out$ok)
  expect_true("E_DOSE_ATTR_DOUBLE_USE" %in% out$diagnostics$code)

  # The remediation must name the *mapping*, not a rename: on the analytical
  # engine the name is inert and the mapping follows the parameter, so
  # "rename it" is the wrong repair (ferx-core #1004).
  msg <- out$diagnostics$message[out$diagnostics$code == "E_DOSE_ATTR_DOUBLE_USE"]
  expect_match(msg, "mapping", fixed = TRUE, all = FALSE)
})

test_that("ferx_fit() refuses to run the double-use model rather than fitting it", {
  ex   <- ferx_example("warfarin")
  path <- .dad_write(.dad_model())
  on.exit(unlink(path))

  # A parse error is terminal, so this returns immediately - no fit is run.
  expect_error(
    suppressWarnings(ferx_fit(path, ex$data, verbose = FALSE)),
    "dose attribute|DOSE_ATTR"
  )
})

test_that("a parameter merely named `F` that no pk() argument maps stays ordinary", {
  # The apparent-parameter convention (`CL/F`, `V/F`) must keep working: without
  # the `f=` mapping, `F` is not a dose attribute on the analytical engine.
  path <- .dad_write(.dad_model(pk_args = "cl=CL, v=V, ka=KA"))
  on.exit(unlink(path))

  expect_false("E_DOSE_ATTR_DOUBLE_USE" %in% .dad_codes(path))
})

test_that("an [initial_conditions] read of a mapped `F` is accepted", {
  # Deliberately not a rejection surface: an initial condition is not an
  # absorbed dose, so the engine seeds the amount with F = 1 and no lag and
  # `init(depot) = F * 500` applies F exactly once (ferx-core #1004).
  path <- .dad_write(.dad_model(
    readout = c("[initial_conditions]", "  init(depot) = F * 500.0")
  ))
  on.exit(unlink(path))

  expect_false("E_DOSE_ATTR_DOUBLE_USE" %in% .dad_codes(path))
})

test_that("the `alag=` spelling is bound and rejected like `lagtime=`", {
  # `lagtime` and `alag` are aliases for one slot, so both mappings arm the
  # check - and the message quotes the spelling the model actually used.
  lines <- .dad_model(
    pk_args = "cl=CL, v=V, ka=KA, alag=TLAG",
    readout = c("[scaling]", "  obs_scale = V / TLAG")
  )
  lines <- sub("^  F  = TVF$", "  TLAG = TVF", lines)
  path  <- .dad_write(lines)
  on.exit(unlink(path))

  out <- .dad_validate(path)
  expect_true("E_DOSE_ATTR_DOUBLE_USE" %in% out$diagnostics$code)
  msg <- out$diagnostics$message[out$diagnostics$code == "E_DOSE_ATTR_DOUBLE_USE"]
  expect_match(msg, "alag=TLAG", fixed = TRUE, all = FALSE)
})
