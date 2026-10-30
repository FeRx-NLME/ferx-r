
# ---- header from test-model-validate.R ----
# Helper: write a minimal valid warfarin-like .ferx file
write_ferx <- function(sections, path = tempfile(fileext = ".ferx")) {
  writeLines(sections, path)
  path
}

VALID_WARFARIN_SECTIONS <- c(
  "[parameters]",
  "  theta TVCL(0.134, 0.001, 10.0)",
  "  omega ETA_CL ~ 0.07",
  "  sigma PROP_ERR ~ 0.01",
  "",
  "[individual_parameters]",
  "  CL = TVCL * exp(ETA_CL)",
  "",
  "[structural_model]",
  "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
  "",
  "[error_model]",
  "  DV ~ proportional(PROP_ERR)"
)



test_that("valid warfarin model prints 'Validating:' header and returns invisibly", {
  ex  <- ferx_example("warfarin")
  out <- capture.output(res <- withVisible(ferx_model_validate(ex$model)))
  expect_true(any(grepl("Validating:", out)))
  expect_false(res$visible)
})
test_that("valid model returns ok = TRUE", {
  ex  <- ferx_example("warfarin")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
  expect_s3_class(res$diagnostics, "data.frame")
  expect_named(res$diagnostics,
               c("severity", "code", "message", "block", "line", "suggestion"))
})
test_that("valid model: $missing_sections is empty (no missing sections reported)", {
  # ferx_model_validate returns a scalar boolean, not a list.
  # We verify indirectly: output contains no [MISSING] markers.
  ex  <- ferx_example("warfarin")
  out <- capture.output(ferx_model_validate(ex$model))
  expect_false(any(grepl("\\[MISSING\\]", out)))
})
test_that("valid model: $unknown_sections is empty (no unknown section markers)", {
  ex  <- ferx_example("warfarin")
  out <- capture.output(ferx_model_validate(ex$model))
  expect_false(any(grepl("\\[unknown section\\]", out)))
})
test_that("model missing a required section returns FALSE", {
  path <- write_ferx(grep("^\\[error_model\\]", VALID_WARFARIN_SECTIONS,
                          value = TRUE, invert = TRUE))
  on.exit(unlink(path))
  res <- ferx_model_validate(path)
  expect_false(isTRUE(res$ok))
})
test_that("model missing a required section reports the missing section name", {
  path <- write_ferx(grep("^\\[error_model\\]|^  DV ~", VALID_WARFARIN_SECTIONS,
                          value = TRUE, invert = TRUE))
  on.exit(unlink(path))
  out <- capture.output(ferx_model_validate(path))
  expect_true(any(grepl("error_model", out)))
})
test_that("model with unknown section reports it", {
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[foo]", "  bar = 1")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  out <- capture.output(ferx_model_validate(path))
  expect_true(any(grepl("unknown section|foo", out)))
})
test_that("model with unknown section is not ok", {
  # Regression for ferx-core #1040: `unknown` was printed as `[unknown section]`
  # but left out of the returned status, so `res$ok` was TRUE for a model
  # carrying a block the engine ignores.
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[foo]", "  bar = 1")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  res <- suppressWarnings(capture.output(r <- ferx_model_validate(path)))
  expect_false(isTRUE(r$ok))
})
test_that("the engine rejects an unknown section with E_UNKNOWN_BLOCK", {
  # The status above must not rest on the R-side list alone: the engine itself
  # refuses the header, so `ok` stays FALSE even if the R check were removed.
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[fit_option]", "  method = focei")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  invisible(capture.output(res <- ferx_model_validate(path)))
  expect_false(isTRUE(res$ok))
  expect_true("E_UNKNOWN_BLOCK" %in% res$diagnostics$code)
  d <- res$diagnostics[res$diagnostics$code == "E_UNKNOWN_BLOCK", , drop = FALSE]
  expect_identical(d$block[[1L]], "fit_option")
  expect_true(grepl("fit_options", d$suggestion[[1L]], fixed = TRUE))
})
test_that("optional sections come from the engine, not a duplicated R list", {
  # The hard-coded vector this replaces had drifted: it omitted `covariates`
  # (and event_model, simulation, ...) so a valid bundled example printed
  # `covariates [unknown section]`, and it listed `initial_values`, which the
  # engine dropped long ago (ferx-core e5e934d).
  known <- ferx:::ferx_rust_known_blocks()
  expect_true(all(c("parameters", "individual_parameters", "structural_model",
                    "error_model", "covariates", "event_model", "simulation",
                    "adaptive_dosing", "mixture", "data_selection") %in% known))
  expect_false("initial_values" %in% known)
  expect_identical(known, sort(known))
})
test_that("a bundled example with a [covariates] block validates clean", {
  ex <- ferx_example("two_cpt_oral_cov")
  out <- capture.output(res <- ferx_model_validate(ex$model))
  expect_true(isTRUE(res$ok))
  expect_false(any(grepl("unknown section", out, fixed = TRUE)))
})
test_that("an instance-named header opens its own block (#306)", {
  # `[event_model cause_a]` used to be read as a body line, so it and its whole
  # body were appended to [parameters] and `event_model` was never seen.
  blocks <- ferx:::.ferx_extract_blocks(ferx_example("tte_competing_risks")$model)
  expect_true(all(c("parameters", "event_model", "fit_options") %in% names(blocks)))
  expect_false(any(grepl("^\\[|^(cmt|family|scale)\\b", blocks$parameters)))
  # Both instances are keyed on the type, in file order.
  expect_identical(grep("^cmt", blocks$event_model, value = TRUE),
                   c("cmt    = 2", "cmt    = 3"))
  expect_identical(blocks$fit_options, c("method  = focei", "maxiter = 300"))
})
test_that("a time-to-event or binary model needs no PK sections (#306)", {
  # The required set depends on the model family, which only the engine knows.
  # An R-side list reported these as missing [individual_parameters],
  # [structural_model] and [error_model].
  for (nm in c("tte_competing_risks", "tte_weibull", "binary_logistic")) {
    out <- capture.output(res <- ferx_model_validate(ferx_example(nm)$model))
    expect_true(isTRUE(res$ok), label = nm)
    expect_false(any(grepl("MISSING", out, fixed = TRUE)), label = nm)
  }
})
test_that("a missing section is the engine's E_MISSING_BLOCK, marked in the report (#306)", {
  path <- write_ferx(grep("^\\[error_model\\]|^  DV ~", VALID_WARFARIN_SECTIONS,
                          value = TRUE, invert = TRUE))
  on.exit(unlink(path))
  out <- capture.output(res <- ferx_model_validate(path))
  expect_false(isTRUE(res$ok))
  d <- res$diagnostics[res$diagnostics$code == "E_MISSING_BLOCK", , drop = FALSE]
  expect_identical(d$block, "error_model")
  expect_true(any(grepl("^  error_model +\\[MISSING\\]$", out)))
  # Present sections are not marked missing.
  expect_true(any(grepl("^  structural_model +\\[ok\\]$", out)))
})
test_that("a retired section is labelled retired, not unknown", {
  # `ferx_rust_known_blocks()` lists only what this build *accepts*, so a
  # retired name lands in the R-side `unknown` set. Printing `[unknown section]`
  # for it contradicts the engine's own E_DEPRECATED_BLOCK line right below,
  # which says the block was retired and should be deleted.
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[initial_values]", "  TVCL = 0.134")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  out <- capture.output(res <- ferx_model_validate(path))
  expect_false(isTRUE(res$ok))
  expect_true("E_DEPRECATED_BLOCK" %in% res$diagnostics$code)
  expect_true(any(grepl("initial_values +\\[retired section\\]", out)))
  expect_false(any(grepl("unknown section", out, fixed = TRUE)))
})
test_that("a feature-gated section is labelled disabled, not unknown", {
  # `markov_model` is a name the engine recognises in every build but only
  # accepts with `--features markov`, which this package does not enable
  # (src/Makevars builds ferx-core with ci,nn,survival by default). It is absent
  # from `ferx_rust_known_blocks()` while the parser still reports it precisely.
  skip_if("markov_model" %in% ferx:::ferx_rust_known_blocks(),
          "this build enables the markov feature")
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[markov_model]", "  states A, B")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  out <- capture.output(res <- ferx_model_validate(path))
  expect_false(isTRUE(res$ok))
  expect_true("E_BLOCK_FEATURE_DISABLED" %in% res$diagnostics$code)
  expect_true(any(grepl("markov_model +\\[feature not enabled\\]", out)))
  expect_false(any(grepl("unknown section", out, fixed = TRUE)))
})
test_that("syntax error in model body returns FALSE and prints errors", {
  lines <- c(
    "[parameters]",
    "  theta !!!BAD SYNTAX!!!",
    "",
    "[individual_parameters]",
    "  CL = 1.0",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "",
    "[error_model]",
    "  DV ~ proportional(0.01)"
  )
  path <- write_ferx(lines)
  on.exit(unlink(path))
  res <- ferx_model_validate(path)
  expect_false(isTRUE(res$ok))
})
test_that("data-dependent checks fire only when data is supplied", {
  # Reference a covariate that is not in the warfarin dataset.
  lines <- c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  omega ETA_CL ~ 0.07",
    "  sigma PROP_ERR ~ 0.01",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL) * (NOT_A_COL / 70)",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  )
  path <- write_ferx(lines)
  on.exit(unlink(path))
  ex <- ferx_example("warfarin")

  # Without data: the missing-covariate check cannot fire.
  res_no_data <- ferx_model_validate(path)
  expect_false(
    "E_MISSING_COVARIATE" %in% res_no_data$diagnostics$code
  )

  # With data: the check fires with a stable error code.
  res_with_data <- ferx_model_validate(path, data = ex$data)
  expect_false(isTRUE(res_with_data$ok))
  expect_true("E_MISSING_COVARIATE" %in% res_with_data$diagnostics$code)
})
test_that("diagnostics carry stable codes and severity labels", {
  lines <- c(
    "[parameters]",
    "  theta !!!BAD!!!",
    "[individual_parameters]",
    "  CL = 1.0",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "[error_model]",
    "  DV ~ proportional(0.01)"
  )
  path <- write_ferx(lines)
  on.exit(unlink(path))
  res <- ferx_model_validate(path)
  expect_false(isTRUE(res$ok))
  expect_gt(nrow(res$diagnostics), 0L)
  expect_true(all(res$diagnostics$severity %in% c("error", "warning")))
  expect_true(all(grepl("^[EW]_", res$diagnostics$code)))
})
test_that("ferx_model_validate errors on missing file", {
  expect_error(
    ferx_model_validate("no_such_file.ferx"),
    regexp = "File not found",
    ignore.case = TRUE
  )
})
test_that("ferx_model_validate errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))
  expect_error(
    ferx_model_validate(path),
    regexp = "\\.ferx",
    ignore.case = TRUE
  )
})
test_that("[derived] and [output] sections are accepted as optional, not flagged as unknown", {
  path <- write_ferx(c(
    VALID_WARFARIN_SECTIONS,
    "",
    "[derived]",
    "  KE = CL / 10",
    "",
    "[output]",
    "  CL"
  ))
  on.exit(unlink(path))
  out <- capture.output(suppressWarnings(ferx_model_validate(path)))
  expect_false(any(grepl("unknown section", out, fixed = TRUE)),
               label = "[derived] and [output] must not appear as unknown sections")
  expect_true(any(grepl("derived", out, ignore.case = TRUE)),
              label = "[derived] section should be printed as present")
  expect_true(any(grepl("output", out, ignore.case = TRUE)),
              label = "[output] section should be printed as present")
})
test_that("bundled igd_inverse_gaussian example validates (igd() parses against the engine)", {
  # Guards that the bundled inverse-Gaussian example stays parseable by the
  # pinned ferx-core: validation runs the model through the Rust engine, so a
  # valid result confirms the built-in `igd(mat, cv2)` input rate is understood
  # (the data/alias are covered by test-example.R; a full fit is skipped — the
  # anchor data is mis-specified for IG, so it stalls rather than converges).
  ex  <- ferx_example("igd_inverse_gaussian")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled biphasic_igd_absorption example validates (FR1*igd + FR2*igd parses)", {
  # Guards the shared pathway-fraction multiplier (ferx-core #388): a valid result
  # confirms the engine accepts a declared-parameter fraction (`FR*igd(...)`) and
  # two input-rate terms on one compartment, and that the fit-init fraction checks
  # (0 < FR <= 1, Σ FR ≈ 1) pass for FR1 / FR2 = 1 - FR1. The matched, recovering
  # fit is exercised by the ferx-core NONMEM anchor (freijer_biphasic_ig), not here.
  ex  <- ferx_example("biphasic_igd_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled weibull_absorption example validates (weibull() parses against the engine)", {
  # Guards that the bundled Weibull example stays parseable by the pinned
  # ferx-core: a valid result confirms the built-in `weibull(td, beta)` input
  # rate is understood by the engine (data/alias covered by test-example.R; a
  # full fit is exercised by the ferx-core NONMEM anchor, not here).
  ex  <- ferx_example("weibull_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled zero_order_absorption example validates (zero_order() parses against the engine)", {
  # Guards that the bundled zero-order example stays parseable by the pinned
  # ferx-core: a valid result confirms the built-in `zero_order(dur)` input rate
  # is understood by the engine (data/alias covered by test-example.R; the value
  # path is anchored against a NONMEM ADVAN1 modeled-duration run in ferx-core).
  ex  <- ferx_example("zero_order_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled parallel_absorption example validates (FR1*first_order + FR2*first_order parses)", {
  # Guards parallel dual-first-order absorption (ferx-core #505): a valid result
  # confirms the engine accepts the built-in `first_order(ka)` input rate composed
  # as two fraction-weighted pathways on one compartment, and that the fit-init
  # fraction checks (0 < FR <= 1, sum FR ~ 1) pass for FR1 / FR2 = 1 - FR1. The
  # matched, recovering fit is exercised by the ferx-core NONMEM anchor, not here.
  ex  <- ferx_example("parallel_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled mixed_absorption example validates (zero_order + first_order compose with fractions)", {
  # Guards mixed zero-order + first-order absorption (ferx-core #505): a valid
  # result confirms a fraction-weighted `zero_order(dur)` and `first_order(ka)` on
  # one compartment parse against the pinned engine -- the per-segment zero-order
  # channel now carries the pathway fraction (anchored to NONMEM in ferx-core).
  ex  <- ferx_example("mixed_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled per_route_lag_absorption example validates (per-route lag= parses on a parallel route)", {
  # Guards per-route absorption lag (ferx-core #856): a valid result confirms the
  # pinned engine accepts the optional `lag=` argument on a fraction-weighted
  # `first_order(ka, lag=LAG2)` pathway (IR + delayed-release), a strict superset of
  # the parallel first_order composition. The matched, recovering fit is the
  # ferx-core NONMEM $DES anchor (per_route_lag.ctl), not exercised here.
  ex  <- ferx_example("per_route_lag_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled sequential_absorption example validates (zero_order + first-order compose)", {
  # Guards the sequential (zero-order depot fill -> first-order ka) composition:
  # a valid result confirms zero_order(dur) into a depot plus a hand-written
  # `- KA*depot` transfer parses against the engine (its own simulated dataset is
  # covered by test-example.R; a full fit recovers the simulation truth offline).
  ex  <- ferx_example("sequential_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})

# -- Start-side bound checks (ferx-core #1251, pinned at 944cbf1e) --
#
# These pin the breaking half of the engine bump: a `theta` whose initial
# estimate is *strictly* outside its own declared range used to be clamped onto
# the bound and fitted from there, silently and on every run. Nothing in this
# package asserted the new refusal, so a pin move that reverted or reworded it
# would leave NEWS.md describing behaviour ferx no longer has, with green CI.
#
# `bad_theta()` swaps only the theta line of the known-good model above, so a
# failure here is about the bound check and not about the rest of the file.
bad_theta <- function(line) {
  s <- VALID_WARFARIN_SECTIONS
  s[grep("^  theta ", s)[1]] <- line
  write_ferx(s)
}

test_that("a theta starting outside its declared range is refused, either side", {
  for (line in c("  theta TVCL(0.05, 0.1, 10.0)",     # below the lower bound
                 "  theta TVCL(50.0, 0.1, 10.0)")) {  # above the upper bound
    path <- bad_theta(line)
    on.exit(unlink(path), add = TRUE)
    capture.output(res <- ferx_model_validate(path))
    expect_false(isTRUE(res$ok), info = line)
    expect_true("E_THETA_INIT_OUTSIDE_BOUNDS" %in% res$diagnostics$code,
                info = line)
  }
})

test_that("a theta starting exactly ON a bound is left alone", {
  # The geometry that separates the implemented rule from the obvious
  # alternative: the check is "strictly outside", not "outside or equal". On the
  # bound the clamp is a no-op, so there is nothing to report and the model is
  # valid. Without this case a check that also rejected `==` would satisfy the
  # test above just as well.
  for (line in c("  theta TVCL(0.1, 0.1, 10.0)",
                 "  theta TVCL(10.0, 0.1, 10.0)")) {
    path <- bad_theta(line)
    on.exit(unlink(path), add = TRUE)
    capture.output(res <- ferx_model_validate(path))
    expect_true(isTRUE(res$ok), info = line)
    # `ok` alone is the baseline model's answer too, so it cannot tell "left
    # alone" from "reported at a severity that still leaves ok TRUE". #1251
    # already splits this behaviour across an error and a warning by severity,
    # so that drift is not hypothetical - assert the silence directly.
    codes <- if (is.data.frame(res$diagnostics)) res$diagnostics$code else character(0)
    expect_false(any(grepl("INIT_OUTSIDE_BOUNDS", codes, fixed = TRUE)), info = line)
  }
})

test_that("a theta whose declared range is empty is reported, not a panic", {
  # Bounds swapped: an empty optimizer box. This used to panic in the bound
  # clamp; it is now E_INIT_BOUNDS_INVERTED and, unlike the check above, carries
  # no `maxiter = 0` exemption.
  path <- bad_theta("  theta TVCL(1.0, 5.0, 2.0)")
  on.exit(unlink(path))
  capture.output(res <- ferx_model_validate(path))
  expect_false(isTRUE(res$ok))
  expect_true("E_INIT_BOUNDS_INVERTED" %in% res$diagnostics$code)
})

# -- The refused fit carries the same code validation prints (#367) -----------

test_that("ferx_fit() attaches the engine's diagnostic code to a refused fit", {
  # `ferx_model_validate()` returned `E_UNKNOWN_BLOCK` while fitting the same
  # file raised the same prose without it, so a script could branch on the code
  # only on the path most users reach second.
  ex    <- ferx_example("warfarin")
  lines <- c(readLines(ex$model), "", "[fit_option]", "  method = focei")
  path  <- withr::local_tempfile(fileext = ".ferx")
  writeLines(lines, path)

  err <- tryCatch(
    suppressWarnings(ferx_fit(path, ex$data, covariance = FALSE, verbose = FALSE)),
    error = function(e) e
  )
  expect_s3_class(err, "ferx_engine_error")
  expect_identical(err$code, "E_UNKNOWN_BLOCK")
  expect_identical(err$block, "fit_option")
  expect_true(grepl("fit_options", err$suggestion, fixed = TRUE))
  # The engine's own prose survives at the front, so existing matches still
  # work, and the code is readable from the console message too.
  expect_true(grepl("Unknown block", conditionMessage(err), fixed = TRUE))
  expect_true(grepl("[E_UNKNOWN_BLOCK]", conditionMessage(err), fixed = TRUE))

  # And it is the same identifier the validator reports for that file.
  invisible(capture.output(res <- ferx_model_validate(path)))
  expect_true(err$code %in% res$diagnostics$code)
})

test_that(".ferx_engine_error() leaves a condition it cannot attribute alone", {
  # The fallback matters as much as the happy path: an engine failure the
  # validation pass does not explain must reach the caller as the engine wrote
  # it, not wearing a code taken from an unrelated finding. A clean model has
  # no error diagnostics at all, so nothing can be attributed.
  ex <- ferx_example("warfarin")
  e  <- simpleError("Fit error: something the validator does not check")
  out <- ferx:::.ferx_engine_error(e, ex$model, ex$data)
  expect_identical(out, e)
  expect_false(inherits(out, "ferx_engine_error"))
  expect_null(out$code)
})

# A refusal of the call itself is not a model or data failure (#517). On a
# model and data with exactly one validation error - warfarin_iov with neither
# an `iov_column` nor an OCC column, E_IOV_MISSING_OCC - the single-error
# fallback used to label a bad `settings` entry with that error's code.
# Mutations that redden this: the glue raising these refusals without
# `argument_refusal()` (the fallback attaches E_IOV_MISSING_OCC), or
# `.ferx_engine_coded_error()` without its empty-code return (the condition
# becomes a `ferx_engine_error` with an empty code).
test_that("ferx_fit() leaves its own argument refusals uncoded (#517)", {
  ex <- ferx_example("warfarin_iov")
  model <- withr::local_tempfile(fileext = ".ferx")
  text <- readLines(ex$model)
  writeLines(text[!grepl("iov_column", text, fixed = TRUE)], model)
  rows <- utils::read.csv(ex$data, na.strings = ".")
  rows$OCC <- NULL
  data <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
  invisible(capture.output(res <- ferx_model_validate(model, data)))
  expect_identical(res$diagnostics$code[res$diagnostics$severity == "error"],
                   "E_IOV_MISSING_OCC")

  plain <- function(expr, pattern) {
    e <- tryCatch({
      suppressWarnings(expr)
      NULL
    }, error = function(e) e)
    expect_s3_class(e, "error")
    expect_false(inherits(e, "ferx_engine_error"), info = pattern)
    expect_null(e$code)
    msg <- if (inherits(e, "condition")) conditionMessage(e) else ""
    expect_match(msg, pattern, fixed = TRUE)
    expect_no_match(msg, "[E_", fixed = TRUE)
  }
  fit <- function(...) ferx_fit(model, data, covariance = FALSE, verbose = FALSE, ...)
  plain(fit(settings = list(iov_occasion = "dose", sir = TRUE)),
        "setting `sir` conflicts with a dedicated ferx_fit() argument")
  plain(fit(settings = list(iov_occasion = "dose", not_a_setting = 1)),
        "unknown fit setting `not_a_setting`")
  plain(fit(settings = list(iov_occasion = "weekly")),
        "fit option `iov_occasion`: unknown value `weekly`")

  # Two refusals of the call only the raw glue can reach: a method token the
  # engine's parser refuses, and `settings` keys without values (R builds the
  # two vectors together). Each record is taken straight after its raise: the
  # next raise overwrites the slot.
  glue_refusal <- function(method, keys, values, pattern) {
    raised <- tryCatch(
      ferx:::ferx_rust_fit(model, data, method, "false", "false", "", 0L, "", "",
                           "", keys, values),
      error = function(e) conditionMessage(e)
    )
    d <- ferx:::ferx_rust_take_engine_diagnostic()
    expect_match(raised, pattern, fixed = TRUE)
    expect_identical(d$code, "", info = pattern)
    expect_identical(charToRaw(d$text %||% ""), charToRaw(raised), info = pattern)
  }
  glue_refusal("agq", character(), character(), "method = agq has been removed")
  glue_refusal(character(), "maxiter", character(),
               "settings keys/values length mismatch (1 vs 0)")

  # The model's own refusal still gets its code: without a rule anywhere the
  # fit is refused for want of occasion labels.
  e <- tryCatch(suppressWarnings(fit()), error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_IOV_MISSING_OCC")
})

# R truncates every error the glue raises at `warning.length` (`Rf_error()`
# formats into a buffer that size), so a refusal longer than that arrives as
# a prefix of the text the glue recorded (#542).
test_that("the raised text matches its record when R truncated it (#542)", {
  match <- ferx:::.ferx_is_raised_text
  text <- paste0(strrep("x", 1200L), " [tail]")
  expect_true(match(text, text))
  expect_true(match(text, substr(text, 1L, 999L)))
  # In a multibyte locale R backs off to a character boundary, up to three
  # bytes short of the limit.
  expect_true(match(text, substr(text, 1L, 996L)))
  # Anything else that is a prefix is a different message.
  expect_false(match(text, substr(text, 1L, 995L)))
  expect_false(match(text, substr(text, 1L, 1000L)))
  expect_false(match(text, "x"))
  expect_false(match("short", "short message"))
  expect_false(match(text, paste0("y", substr(text, 2L, 999L))))
  # Byte-wise: an em-dash cut in half is still a prefix of the record.
  dash <- paste0(strrep("x", 998L), "\u2014", strrep("x", 100L))
  cut <- rawToChar(charToRaw(dash)[1:999])
  expect_true(match(dash, cut))
  # The limit is the session's `warning.length`.
  withr::local_options(warning.length = 2000L)
  expect_false(match(text, substr(text, 1L, 999L)))
  expect_true(match(paste0(text, strrep("z", 1000L)),
                    substr(paste0(text, strrep("z", 1000L)), 1L, 1999L)))
})

test_that("a truncated coded refusal keeps its code and its full text (#542)", {
  text <- paste0("Error simulating: ", strrep("long remediation ", 70L), "end.")
  expect_gt(nchar(text, "bytes"), 1000L)
  testthat::local_mocked_bindings(
    ferx_rust_take_engine_diagnostic = function() {
      list(text = text, message = text, code = "E_TEST", block = "",
           line = 0L, suggestion = "")
    },
    .package = "ferx"
  )
  raised <- substr(text, 1L, 999L)
  e <- ferx:::.ferx_engine_coded_error(simpleError(raised), raised)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_TEST")
  expect_identical(conditionMessage(e), paste0(text, " [E_TEST]"))
})

test_that("a truncated argument refusal stays uncoded (#542, #517)", {
  # The glue's refusal of an unknown setting has an empty code. Truncated, it
  # used to miss its record and reach the single-error fallback, which labelled
  # it with the model's one unrelated finding (E_IOV_MISSING_OCC).
  ex <- ferx_example("warfarin_iov")
  model <- withr::local_tempfile(fileext = ".ferx")
  text <- readLines(ex$model)
  writeLines(text[!grepl("iov_column", text, fixed = TRUE)], model)
  rows <- utils::read.csv(ex$data, na.strings = ".")
  rows$OCC <- NULL
  data <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")

  key <- paste0("not_a_setting_", strrep("x", 1200L))
  settings <- list(iov_occasion = "dose", 1)
  names(settings)[2L] <- key
  e <- tryCatch({
    suppressWarnings(ferx_fit(model, data, covariance = FALSE, verbose = FALSE,
                              settings = settings))
    NULL
  }, error = function(e) e)
  expect_s3_class(e, "error")
  expect_false(inherits(e, "ferx_engine_error"))
  expect_null(e$code)
  msg <- conditionMessage(e)
  expect_match(msg, "unknown fit setting `not_a_setting_", fixed = TRUE)
  expect_no_match(msg, "[E_", fixed = TRUE)
})

test_that("a class-blind tryCatch still sees an ordinary error", {
  ex    <- ferx_example("warfarin")
  lines <- c(readLines(ex$model), "", "[fit_option]", "  method = focei")
  path  <- withr::local_tempfile(fileext = ".ferx")
  writeLines(lines, path)
  expect_error(
    suppressWarnings(ferx_fit(path, ex$data, covariance = FALSE, verbose = FALSE)),
    "Unknown block"
  )
})

test_that("a model file gone before the data bind is bound from the one read (ferx-core #1752)", {
  skip_on_os("windows")  # needs a FIFO
  skip_if(!nzchar(Sys.which("mkfifo")), "mkfifo not available")
  ex    <- ferx_example("warfarin")
  dir   <- withr::local_tempdir()
  model <- file.path(dir, "m.ferx")
  file.copy(ex$model, model)
  # The data file is a FIFO, so the writer decides when the read finishes: it
  # deletes the model while the validator blocks on the pipe, after the parse
  # and before the re-read the data bind does, then hands over the rows.
  # Assumption: the model parse finishes within the writer's 2 s sleep (it
  # takes milliseconds); a slower parse would see the model already gone.
  fifo <- file.path(dir, "data.csv")
  system2("mkfifo", shQuote(fifo))
  system(sprintf("(sleep 2; rm %s; cat %s > %s) &",
                 shQuote(model), shQuote(ex$data), shQuote(fifo)))
  # Bounded unblockers, so a broken validator fails this test rather than
  # hanging the suite or outliving it. A reader that opens the FIFO a second
  # time gets EOF after 30 s; on exit, a non-blocking read-open releases a
  # writer still waiting in open() because the validator never read the pipe.
  system(sprintf("(sleep 30; : <> %s) 2>/dev/null &", shQuote(fifo)))
  withr::defer(close(fifo(fifo, "r", blocking = FALSE)))
  utils::capture.output(res <- ferx_model_validate(model, fifo))
  errs <- res$diagnostics[res$diagnostics$severity == "error", ]
  # The check reads the model once and binds the text it parsed, so a model
  # deleted after that read no longer matters. Before ferx-core 712cd47b this
  # was E_THETA_LEVEL_BINDING on `parameters`, at 712cd47b E_MODEL_REREAD:
  # the data bind read the file a second time.
  expect_identical(errs$code, character(0))
  expect_true(isTRUE(res$ok))
})

test_that("ferx_fit() refuses a theta starting outside its declared range", {
  # ferx_model_validate() is advisory; this is the assertion that pins the
  # *breaking* change, because a fit is what a user's script actually runs.
  ex   <- ferx_example("warfarin")
  path <- bad_theta("  theta TVCL(0.05, 0.1, 10.0)")
  on.exit(unlink(path))
  expect_error(
    ferx_fit(path, ex$data, covariance = FALSE, verbose = FALSE),
    "declared lower bound"
  )
})

# -- Failures that are not the model or the data (#521) ----------------------

# warfarin_iov with neither an `iov_column` nor an OCC column: validation finds
# exactly one error, E_IOV_MISSING_OCC, so the single-error fallback in
# `.ferx_engine_error()` has a code to hand any failure it reaches.
iov_without_occasions <- function(env = parent.frame()) {
  ex <- ferx_example("warfarin_iov")
  model <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  text <- readLines(ex$model)
  writeLines(text[!grepl("iov_column", text, fixed = TRUE)], model)
  rows <- utils::read.csv(ex$data, na.strings = ".")
  rows$OCC <- NULL
  data <- withr::local_tempfile(fileext = ".csv", .local_envir = env)
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
  list(model = model, data = data)
}

fit_failure <- function(expr) {
  tryCatch({
    suppressWarnings(expr)
    NULL
  }, error = function(e) e)
}

# Ctrl-C and a worker panic cannot be raised on demand, so `ferx_rust_fit()`
# raises their text here and the record the glue keeps beside it. Mutations
# that redden this: the glue recording either without its kind (it reaches the
# fallback and comes back `... [E_IOV_MISSING_OCC]`), or
# `.ferx_failure_condition()` returning the condition unclassed.
test_that("a cancelled or panicked fit is not labelled with a validation code (#521)", {
  f <- iov_without_occasions()
  for (case in list(
    list(kind = "cancelled", text = "ferx_fit: cancelled by user",
         class = "ferx_cancelled"),
    list(kind = "panic", text = "Fit error: worker thread panicked",
         class = "ferx_panic")
  )) {
    testthat::local_mocked_bindings(
      ferx_rust_fit = function(...) stop(case$text, call. = FALSE),
      ferx_rust_take_engine_diagnostic = function() {
        list(kind = case$kind, text = case$text, message = case$text,
             code = "", block = "", line = 0L, suggestion = "")
      },
      ferx_rust_validate_model = function(...) stop("re-validation was called"),
      .package = "ferx"
    )
    e <- fit_failure(ferx_fit(f$model, f$data, covariance = FALSE,
                              verbose = FALSE))
    expect_s3_class(e, case$class)
    expect_s3_class(e, "error")
    expect_false(inherits(e, "ferx_engine_error"), info = case$kind)
    expect_null(e$code)
    expect_identical(conditionMessage(e), case$text)
  }
})

# The glue side of the test above, with nothing mocked: the arms `ferx_rust_fit`
# takes for a cancelled run and a panicked worker (`fit_outcome`, driven by a
# test hook because neither can be caused on demand), and a panic `entry`
# catches in any entry point. Each must leave a record of its kind, and R must
# turn that record into the classed condition. Mutation that reddens this: any
# of the three recording no kind (or `Call`).
test_that("the glue records a cancel and a panic with their kind (#521)", {
  for (case in list(
    list(raise = function() ferx:::ferx_rust_test_fit_outcome("cancelled"),
         kind = "cancelled", text = "ferx_fit: cancelled by user",
         class = "ferx_cancelled"),
    list(raise = function() ferx:::ferx_rust_test_fit_outcome("panicked"),
         kind = "panic", text = "Fit error: worker thread panicked",
         class = "ferx_panic"),
    list(raise = function() ferx:::ferx_rust_test_panic("engine blew up"),
         kind = "panic", text = "engine blew up", class = "ferx_panic")
  )) {
    raised <- tryCatch(case$raise(), error = function(e) e)
    expect_identical(conditionMessage(raised), case$text)
    d <- ferx:::ferx_rust_take_engine_diagnostic()
    expect_identical(d$kind, case$kind, info = case$text)
    expect_identical(d$code, "", info = case$text)
    expect_identical(charToRaw(d$text), charToRaw(case$text))

    # Through `.ferx_engine_error()`, on the one-error model: classed, uncoded,
    # never re-validated.
    f <- iov_without_occasions()
    raised <- tryCatch(case$raise(), error = function(e) e)
    e <- ferx:::.ferx_engine_error(raised, f$model, f$data)
    expect_s3_class(e, case$class)
    expect_false(inherits(e, "ferx_engine_error"), info = case$text)
    expect_null(e$code)
  }
})

test_that("every glue record names its kind (#521)", {
  f <- iov_without_occasions()
  tryCatch(
    ferx:::ferx_rust_fit(f$model, f$data, character(), "false", "false", "",
                         0L, "", "", "", "maxiter", character()),
    error = function(e) NULL
  )
  expect_identical(ferx:::ferx_rust_take_engine_diagnostic()$kind, "call")

  ex <- ferx_example("warfarin")
  tryCatch(ferx:::ferx_rust_predict(ex$model, infusion_into_cmt0_data()),
           error = function(e) NULL)
  d <- ferx:::ferx_rust_take_engine_diagnostic()
  expect_identical(d$code, "E_DOSE_CMT_NOT_INFUSABLE")
  expect_identical(d$kind, "engine")
})

# The second half: the fallback re-validates the model file and the data, and
# sees none of the call's settings. With `settings = list(iov_occasion =
# "dose")` the model fits, so a failure of that fit has nothing to do with the
# one finding validation reports without the setting.
test_that("the fallback does not label a fit whose settings validation cannot see (#521)", {
  f <- iov_without_occasions()
  testthat::local_mocked_bindings(
    ferx_rust_fit = function(...) stop("Fit error: an uncoded failure", call. = FALSE),
    ferx_rust_take_engine_diagnostic = function() NULL,
    .package = "ferx"
  )
  fit <- function(...) {
    fit_failure(ferx_fit(f$model, f$data, verbose = FALSE, ...))
  }
  # No override: the fallback still labels the failure, as before.
  e <- fit()
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_IOV_MISSING_OCC")
  # Anything the re-validation does not see: no code, never a wrong one.
  for (args in list(
    list(settings = list(iov_occasion = "dose")),
    list(ignore = "TIME > 100"),
    list(ignore_ids = 1),
    list(method = "foce"),
    list(covariance = FALSE),
    list(inits_from_nca = TRUE)
  )) {
    e <- do.call(fit, args)
    expect_s3_class(e, "error")
    expect_false(inherits(e, "ferx_engine_error"), info = names(args))
    expect_no_match(conditionMessage(e), "[E_", fixed = TRUE)
  }
  # `verbose` and `threads` change nothing a check looks at.
  e <- fit(threads = 1L)
  expect_identical(e$code, "E_IOV_MISSING_OCC")
})

test_that("the issue's model fits with the occasion rule it was refused without (#521)", {
  f <- iov_without_occasions()
  # A selection that removes every observation, on top of the setting: the
  # refusal is the selection's, not E_IOV_MISSING_OCC.
  e <- fit_failure(ferx_fit(f$model, f$data, covariance = FALSE, verbose = FALSE,
                            settings = list(iov_occasion = "dose"),
                            ignore = "EVID == 0"))
  expect_s3_class(e, "ferx_engine_error")
  expect_match(conditionMessage(e), "nothing the likelihood can score",
               fixed = TRUE)
  expect_identical(e$code, "E_NO_SCORED_OBSERVATIONS")
  expect_no_match(conditionMessage(e), "E_IOV_MISSING_OCC", fixed = TRUE)
})

test_that("a model that does not parse is still labelled with call-time settings (#521)", {
  ex    <- ferx_example("warfarin")
  lines <- c(readLines(ex$model), "", "[not_a_block]", "  x = 1")
  path  <- withr::local_tempfile(fileext = ".ferx")
  writeLines(lines, path)
  e <- fit_failure(ferx_fit(path, ex$data, verbose = FALSE, method = "foce",
                            settings = list(maxiter = 5)))
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_UNKNOWN_BLOCK")
})

test_that("an info diagnostic prints as info, not as a warning (#660)", {
  # ferx-core #2427 added `Severity::Info` (W_LEVEL_SUMMARY is the first). The
  # engine result is mocked so every mix of severities is covered; the test
  # below checks the real engine's rows.
  ex <- ferx_example("warfarin")
  validate_with <- function(severity, code) {
    n <- length(severity)
    testthat::local_mocked_bindings(
      ferx_rust_validate_model = function(model_path, data_path) {
        list(ok = TRUE, model = "warfarin", data = "",
             severity = severity, code = code,
             message = paste("finding", code), block = rep("", n),
             line = rep(0L, n), suggestion = rep("", n))
      },
      .package = "ferx"
    )
    out <- capture.output(res <- ferx_model_validate(ex$model))
    list(out = out, res = res)
  }

  v <- validate_with("info", "W_LEVEL_SUMMARY")
  expect_true(v$res$ok)
  expect_identical(v$res$diagnostics$severity, "info")
  expect_true(any(grepl("Result: VALID (with notes)", v$out, fixed = TRUE)))
  expect_false(any(grepl("with warnings", v$out, fixed = TRUE)))
  expect_true(any(grepl("* info W_LEVEL_SUMMARY", v$out, fixed = TRUE)))
  expect_false(any(grepl("* warning W_LEVEL_SUMMARY", v$out, fixed = TRUE)))

  v <- validate_with(c("warning", "info"), c("W_SOMETHING", "W_LEVEL_SUMMARY"))
  expect_true(any(grepl("Result: VALID (with warnings and notes)", v$out,
                        fixed = TRUE)))
  expect_true(any(grepl("* warning W_SOMETHING", v$out, fixed = TRUE)))
  expect_true(any(grepl("* info W_LEVEL_SUMMARY", v$out, fixed = TRUE)))

  v <- validate_with("warning", "W_SOMETHING")
  expect_true(any(grepl("Result: VALID (with warnings)", v$out, fixed = TRUE)))
})

test_that("the engine's W_LEVEL_SUMMARY rows arrive as info notes (ferx-core #2427)", {
  ex <- ferx_example("nested_levels")
  out <- capture.output(v <- ferx_model_validate(ex$model, data = ex$data))
  info <- v$diagnostics[v$diagnostics$code == "W_LEVEL_SUMMARY", , drop = FALSE]
  expect_gt(nrow(info), 0L)
  expect_true(all(info$severity == "info"))
  expect_true(v$ok)
  expect_false(any(v$diagnostics$severity == "warning"))
  expect_true("Result: VALID (with notes)" %in% out)
  expect_true(any(grepl("^  \\* info W_LEVEL_SUMMARY", out)))
})
