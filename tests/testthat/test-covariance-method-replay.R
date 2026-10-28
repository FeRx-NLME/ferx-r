# ferx_covariance(fit) replays the fit's own covariance_method (#659).
#
# Until #659 the argument defaulted to "r", so a fit made with `rsr` or `s` (in
# the call's `settings` or in [fit_options]) was replayed as `r` unless the
# caller repeated it: on mbma_placebo the replayed SEs sat 17-19% from the
# inline rsr ones (TVE0 1.013 -> 0.952, ED50 8.196 -> 9.244). The default is
# now the fit's request, `fit$call_settings` ahead of the model file, and the
# estimator that ran is recorded as `fit$covariance_method`.

q <- function(x) suppressWarnings(suppressMessages(x))
se <- function(f) lapply(f[c("se_theta", "se_omega", "se_sigma", "se_kappa")], unname)

# ferx_covariance() with the conditions it raised about the method kept apart
# from the engine's: those are the ones this file asserts on.
replay <- function(...) {
  msgs <- character()
  warns <- character()
  out <- withCallingHandlers(
    ferx_covariance(...),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  mine <- function(x) x[grepl("covariance_method", x, fixed = TRUE)]
  list(fit = out, msgs = mine(msgs), warns = mine(warns))
}

round_trip <- function(f) {
  path <- withr::local_tempfile(fileext = ".fitrx", .local_envir = parent.frame())
  q(ferx_save_fit(f, path))
  q(ferx_load_fit(path))
}

# warfarin under each estimator, from the call. Its 32 subjects are more than
# its free parameters, so S is invertible and `s` computes too.
warfarin_by_method <- local({
  out <- NULL
  function() {
    if (is.null(out)) {
      ex <- ferx_example("warfarin")
      out <<- lapply(c(r = "r", rsr = "rsr", s = "s"), function(m) {
        q(ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                   covariance = TRUE,
                   settings = list(maxiter = 30L, covariance_method = m)))
      })
    }
    out
  }
})

# mbma_placebo with `covariance_method = <m>` in its [fit_options], profiled
# (its PLACEBO block, unbounded) or not: the case #659 measured.
mbma_file_model <- function(m, profiled = FALSE) {
  ex <- ferx_example("mbma_placebo")
  txt <- readLines(ex$model)
  if (profiled) {
    at <- grep("^  theta PLACEBO\\[", txt)
    stopifnot(length(at) == 1L)
    txt[at] <- paste(sub("(0.0, -50.0, 50.0)", "(0.0)", txt[at], fixed = TRUE),
                     "profile")
  }
  at <- grep("^\\[fit_options\\]", txt)
  stopifnot(length(at) == 1L)
  txt <- append(txt, paste0("  covariance_method = ", m), after = at)
  model <- withr::local_tempfile(fileext = ".ferx", .local_envir = parent.frame())
  writeLines(txt, model)
  list(model = model, data = ex$data)
}

test_that("ferx_covariance(fit) repeats each estimator the call asked for, live and reloaded", {
  skip_on_cran()
  fits <- warfarin_by_method()
  for (m in names(fits)) {
    f <- fits[[m]]
    expect_identical(f$covariance_status, "computed", info = m)
    # The fit records the estimator that ran.
    expect_identical(f$covariance_method, m, info = m)
    got <- replay(f)
    expect_identical(se(got$fit), se(f), info = m)
    expect_identical(got$fit$covariance_method, m, info = m)
    expect_identical(got$msgs, character(), info = m)
    expect_identical(got$warns, character(), info = m)

    loaded <- round_trip(f)
    expect_identical(loaded$covariance_method, m, info = m)
    got <- replay(loaded)
    expect_identical(se(got$fit), se(f), info = m)
    expect_identical(got$fit$covariance_method, m, info = m)
    expect_identical(c(got$msgs, got$warns), character(), info = m)
  }
  # Three estimators, three sets of standard errors: the replays above are not
  # one estimate three times.
  expect_false(isTRUE(all.equal(se(fits$r), se(fits$rsr))))
  expect_false(isTRUE(all.equal(se(fits$r), se(fits$s))))
})

test_that("an explicit covariance_method still overrides the fit's", {
  skip_on_cran()
  fits <- warfarin_by_method()
  # The three fits differ only in the covariance step, so the estimates are
  # one: the rsr fit replayed under an explicit "r" is the r fit's step.
  expect_identical(fits$rsr$theta, fits$r$theta)
  for (m in c("r", "s")) {
    got <- replay(fits$rsr, covariance_method = m)
    expect_identical(se(got$fit), se(fits[[m]]), info = m)
    expect_identical(got$fit$covariance_method, m, info = m)
    # Asked for, so nothing to say about the fit's.
    expect_identical(c(got$msgs, got$warns), character(), info = m)
  }
  got <- replay(fits$r, covariance_method = "SANDWICH")
  expect_identical(se(got$fit), se(fits$rsr))
  expect_identical(got$fit$covariance_method, "rsr")
})

test_that("[fit_options] covariance_method is replayed on mbma_placebo, profiled and not", {
  skip_on_cran()
  for (profiled in c(FALSE, TRUE)) {
    info <- if (profiled) "profiled" else "unprofiled"
    m <- mbma_file_model("rsr", profiled)
    f <- q(ferx_fit(m$model, m$data, verbose = FALSE))
    expect_identical(f$covariance_status, "computed", info = info)
    expect_identical(f$covariance_method, "rsr", info = info)
    expect_null(f$call_settings$covariance_method, info = info)
    got <- replay(f)
    expect_identical(se(got$fit), se(f), info = info)
    expect_identical(got$fit$covariance_method, "rsr", info = info)
    expect_identical(c(got$msgs, got$warns), character(), info = info)
    expect_identical(se(replay(round_trip(f))$fit), se(f), info = info)
    # The glue reads it from the hash-checked model file, so a fit without R's
    # records of the call and the file (a bundle ferx-core wrote) keeps it.
    bare <- f
    bare$call_settings <- list()
    bare$model_file_settings <- list()
    expect_identical(se(replay(bare)$fit), se(f), info = info)
    # The old default, "r", is another answer here (17-19% on the structural
    # SEs at ferx-core 32dbb309).
    r <- replay(f, covariance_method = "r")$fit
    expect_identical(r$covariance_method, "r", info = info)
    expect_gt(abs(r$se_theta[["ED50"]] / f$se_theta[["ED50"]] - 1), 0.05)
  }
})

test_that("an `s` the fit could not form is not formed by the replay either", {
  skip_on_cran()
  # mbma_placebo has six subjects, fewer than its free parameters: S is
  # singular, the fit's step fails, and so does the replay of it.
  m <- mbma_file_model("s")
  f <- q(ferx_fit(m$model, m$data, verbose = FALSE))
  expect_identical(f$covariance_status, "failed")
  expect_identical(f$covariance_method, NA_character_)
  got <- replay(f)
  expect_identical(got$fit$covariance_status, "failed")
  expect_identical(got$fit$covariance_method, NA_character_)
  # No matrix on either side: nothing to label, nothing to note.
  expect_identical(c(got$msgs, got$warns), character())
  # An explicit "r" does form one, so the failure was the method's.
  r <- replay(f, covariance_method = "r")$fit
  expect_identical(r$covariance_status, "computed")
  expect_identical(r$covariance_method, "r")
})

test_that("a fit saved before #659 loads, and falls back to its call, the model file, or r, with a note", {
  skip_on_cran()
  fits <- warfarin_by_method()
  legacy <- function(f) {
    f$covariance_method <- NULL  # what a bundle saved before #659 lacks
    round_trip(f)
  }
  # The call settings travelled in the bundle before #659 too.
  old <- legacy(fits$rsr)
  expect_identical(old$covariance_method, NA_character_)
  got <- replay(old)
  expect_identical(se(got$fit), se(fits$rsr))
  expect_identical(got$fit$covariance_method, "rsr")
  expect_length(got$msgs, 1L)
  expect_match(got$msgs, "does not record which covariance_method", fixed = TRUE)
  expect_match(got$msgs, "before ferx-r #659", fixed = TRUE)
  expect_match(got$msgs, "replaying \"rsr\", from the fit's call settings", fixed = TRUE)
  expect_identical(got$warns, character())

  # With no record of the call either, the default.
  old$call_settings <- list()
  got <- replay(old)
  expect_identical(se(got$fit), se(fits$r))
  expect_match(got$msgs, "replaying \"r\", from the default", fixed = TRUE)
  expect_identical(got$warns, character())

  # The model file's [fit_options] in between.
  m <- mbma_file_model("rsr")
  f <- q(ferx_fit(m$model, m$data, verbose = FALSE))
  old <- legacy(f)
  old$call_settings <- list()
  got <- replay(old)
  expect_identical(se(got$fit), se(f))
  expect_match(got$msgs, "replaying \"rsr\", from the model file's [fit_options]",
               fixed = TRUE)

  # A fit with no matrix records no estimator and needs no note.
  got <- replay(warfarin_fit())
  expect_identical(c(got$msgs, got$warns), character())
})

test_that("a replay that runs another estimator than the fit recorded warns", {
  skip_on_cran()
  fits <- warfarin_by_method()
  # The label survives, the record of the call that asked for it does not.
  f <- fits$rsr
  f$call_settings <- list()
  got <- replay(f)
  expect_identical(se(got$fit), se(fits$r))
  expect_identical(got$fit$covariance_method, "r")
  expect_length(got$warns, 1L)
  expect_match(got$warns, "the fit's covariance step ran covariance_method = \"rsr\"",
               fixed = TRUE)
  expect_match(got$warns, "this replay ran \"r\" (requested \"r\", from the default",
               fixed = TRUE)
  expect_match(got$warns, "its standard errors are not the fit's", fixed = TRUE)
  expect_match(got$warns, "Pass covariance_method = \"rsr\" to repeat the fit's",
               fixed = TRUE)
  expect_identical(got$msgs, character())
})
