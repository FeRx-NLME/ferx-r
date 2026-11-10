# ferx_covariance(fit) and ferx_sir(fit) repeat the in-fit steps (#511, #472).
#
# A fit records how it was scored: `fit$scoring_settings` (the inner-loop and
# ODE settings of the stage that produced the estimates, ferx-core #1805),
# `fit$sir_settings` (what its SIR ran under, ferx-core #1758) and, in memory,
# `fit$packed_estimate` (the optimizer's exact packed vector). The standalone
# steps hand all three back, so with no arguments they are identical() to the
# in-fit steps, and an explicit argument overrides by editing the record.
#
# Every row asserts its premise first - dropping the field it is about moves the
# result on that fixture - so a setting that does not move its fixture cannot
# pass by accident (#511's `inner_tol = 1e-5` is the engine default and
# `mu_referencing = FALSE` is inert on warfarin; both rows moved, plan section 0b).
#
# Engine: the R glue is built FERX_NO_AUTODIFF=1 in the measurements quoted
# here; ferx-core's inner gradients on these fixtures are analytic either way.

ps_scratch <- function(ext, env = parent.frame()) {
  withr::local_tempfile(fileext = ext, .local_envir = env)
}

# SIR fits draw 300, keep 100, seed 1 - the in-fit run each S4/S5 row repeats.
ps_sir_settings <- list(sir_samples = 300L, sir_resamples = 100L, sir_seed = 1L)

# Fits are shared between tests, so each is made once per file run.
ps_cache <- new.env(parent = emptyenv())
ps_fit <- function(key, example, settings = list(), sir = FALSE, ...) {
  if (is.null(ps_cache[[key]])) {
    ex <- if (is.list(example)) example else ferx_example(example)
    if (sir) settings <- c(settings, ps_sir_settings)
    ps_cache[[key]] <- suppressWarnings(ferx_fit(
      ex$model, ex$data, covariance = TRUE, verbose = FALSE, sir = sir,
      settings = if (length(settings)) settings else NULL, ...
    ))
  }
  fit <- ps_cache[[key]]
  skip_if(is.null(fit$cov_matrix), paste(key, "covariance step did not converge"))
  fit
}

ps_round_trip <- function(fit, env = parent.frame()) {
  path <- ps_scratch(".fitrx", env)
  unlink(path)
  ferx_save_fit(fit, path)
  ferx_load_fit(path)
}

# `fit` with the fields `keys` gone from its bytes, which is what the steps read.
# "sir_settings" is the `sir.settings` record.
ps_without <- function(fit, ...) {
  keys <- c(...)
  fitrx_edit_fit(fit, function(w) {
    for (k in setdiff(keys, "sir_settings")) w[[k]] <- NULL
    if ("sir_settings" %in% keys && !is.null(w$sir)) w$sir$settings <- NULL
    w
  })
}

# `fit` with one key of its scoring record set in its bytes (NULL drops it).
ps_scoring <- function(fit, key, value) {
  fitrx_edit_fit(fit, function(w) {
    w$scoring_settings[[key]] <- value
    w
  })
}

ps_sir_fields <- function(fit) {
  fit[c("sir_ess", "sir_ci_theta", "sir_ci_omega", "sir_ci_sigma", "sir_ci_kappa")]
}

ps_scoring_keys <- c(
  "inner_maxiter", "inner_tol", "inner_restarts", "inner_screen",
  "inner_twin_seed", "mu_referencing", "n_agq",
  "inner_optimizer", "ebe_warm_start", "ode_reltol", "ode_abstol",
  "ode_max_steps", "ode_method", "ode_stiff_abort_after", "ode_auto_switch",
  "gradient_method", "nn_l2_lambda", "nn_smooth_lambda"
)

test_that("S11: a fit shows its scoring record, and its SIR record only after SIR", {
  fit <- ps_fit("w_default", "warfarin")
  expect_identical(names(fit$scoring_settings), ps_scoring_keys)
  expect_identical(fit$scoring_settings$inner_optimizer, "auto")
  expect_true(is.numeric(fit$packed_estimate) && length(fit$packed_estimate) > 0L)
  expect_null(fit$sir_settings)

  sir <- ps_fit("w_sir_df3", "warfarin", list(sir_df = 3), sir = TRUE)
  expect_identical(
    names(sir$sir_settings),
    c("samples", "resamples", "seed", "df", "scale", "keep_samples", "scoring")
  )
  expect_identical(names(sir$sir_settings$scoring), ps_scoring_keys)
  expect_identical(sir$sir_settings[c("samples", "resamples", "seed", "df")],
                   list(samples = 300L, resamples = 100L, seed = 1, df = 3))
})

test_that("S12: a gradient = \"fd\" fit stays on finite differences in ferx_covariance(), also after a reload", {
  # ferx-core #1835: the fit records its gradient route, and the post-hoc step
  # takes it as a union with the caller's (either side's fd wins).
  fd <- ps_fit("w_fd", "warfarin", gradient = "fd")
  expect_identical(fd$scoring_settings$gradient_method, "fd")
  expect_identical(ps_fit("w_default", "warfarin")$scoring_settings$gradient_method,
                   "auto")
  # Premise: the analytic route moves this fit's covariance.
  as_auto <- ps_scoring(fd, "gradient_method", "auto")
  expect_false(identical(ferx_covariance(as_auto)$cov_matrix, fd$cov_matrix))

  cov <- ferx_covariance(fd)
  expect_identical(cov$cov_matrix, fd$cov_matrix)
  loaded <- ps_round_trip(fd)
  expect_identical(loaded$scoring_settings, fd$scoring_settings)
  expect_identical(ferx_covariance(loaded)$cov_matrix, fd$cov_matrix)
})

test_that("S13: a scoring record without gradient_method (saved before ferx-core #1835) reads it as auto", {
  # In the SIR record's scoring half, read through the record decoder.
  rec <- ferx:::ferx_rust_settings_record("sir_settings", NULL, "test", "rec")
  rec$scoring$gradient_method <- NULL
  back <- ferx:::ferx_rust_settings_record("sir_settings", rec, "test", "rec")
  expect_identical(back$scoring$gradient_method, "auto")
})

test_that("S1: on a default fit ferx_covariance(fit) is identical() to the in-fit step", {
  # The original 1e-7 case (#511): warfarin_iov SE(TVKA) 0.76168056 vs
  # 0.76168066 was the packed estimate, not the settings (plan section 0a).
  fit <- ps_fit("iov_default", "warfarin_iov")
  # On the finite-difference route the step reads the packed estimate; see S1b.
  fd <- ps_fit("iov_fd", "warfarin_iov", gradient = "fd")
  expect_identical(ferx_covariance(fd)$cov_matrix, fd$cov_matrix)

  cov <- ferx_covariance(fit)
  expect_identical(cov$cov_matrix, fit$cov_matrix)
  expect_identical(cov$se_theta, fit$se_theta)
})

test_that("S1b: on gradient = \"fd\" the packed estimate is what makes S1 identical()", {
  # Premise for S1's finite-difference half. Since ferx-core #2093 the default
  # fit's covariance step is the exact R-matrix, which does not read the packed
  # estimate at all. The FD route does: without the vector, core re-packs the
  # reported estimates, and the FD Hessian amplifies the last-bit difference of
  # that repack (2.8e-5 to 5.8e-5 relative on macOS arm64). Whether the repack
  # differs at all is round-off, so at some optima on some platforms it comes
  # back bit-equal (ubuntu at ferx-core 497a7a9) and the premise cannot be shown
  # there.
  fd <- ps_fit("iov_fd", "warfarin_iov", gradient = "fd")
  no_packed <- ps_without(fd, "packed_estimate")
  moved <- !identical(ferx_covariance(no_packed)$cov_matrix, fd$cov_matrix)
  skip_if(!moved, "the repack of the reported estimates is bit-equal here")
  expect_true(moved)
})

# Fits made without SIR, so `run_covariance`'s fallback to the SIR record's
# scoring half cannot stand in for the record under test.
#
# w_lbfgs and ode_rtol3 run on gradient = "fd". Since ferx-core #2093 their
# default fits take the exact R-matrix, which neither setting moves on these
# fixtures (identical() with and without the record); on the finite-difference
# stencil both still do (measured 2.8e-2 and 3.2e-2 relative). mm_rtol3 keeps
# a premise on the exact R-matrix itself.
ps_cov_rows <- list(
  list(key = "w_lbfgs", example = "warfarin", setting = "inner_optimizer",
       settings = list(inner_optimizer = "lbfgs"), args = list(gradient = "fd"),
       value = "lbfgs"),
  list(key = "w_tol3", example = "warfarin", setting = "inner_tol",
       settings = list(inner_tol = 1e-3), value = 1e-3),
  list(key = "iov_mu", example = "warfarin_iov", setting = "mu_referencing",
       args = list(mu_referencing = FALSE), value = FALSE),
  list(key = "ode_rtol3", example = "warfarin_ode", setting = "ode_reltol",
       settings = list(ode_reltol = 1e-3), args = list(gradient = "fd"), value = 1e-3),
  # The same setting on the default route, where it does reach the exact
  # R-matrix: mm_oral is a nonlinear ODE, so ode_reltol reaches its solve
  # (warfarin_ode takes the exact linear path and ignores it). Dropping the
  # record moved this covariance by 35.9 relative.
  list(key = "mm_rtol3", example = "mm_oral", setting = "ode_reltol",
       settings = list(ode_reltol = 1e-3), value = 1e-3)
)

ps_cov_fit <- function(row) {
  do.call(ps_fit, c(list(row$key, row$example, row$settings %||% list()), row$args))
}

test_that("S2: ferx_covariance(fit) follows the fit's scoring record", {
  for (row in ps_cov_rows) {
    fit <- ps_cov_fit(row)
    # The record holds the setting as the stage ran it.
    expect_identical(fit$scoring_settings[[row$setting]], row$value, info = row$key)
    # Premise: on this fixture the setting moves the covariance.
    no_record <- ps_without(fit, "scoring_settings")
    expect_false(identical(ferx_covariance(no_record)$cov_matrix, fit$cov_matrix),
                 info = row$key)
    cov <- ferx_covariance(fit)
    expect_identical(cov$cov_matrix, fit$cov_matrix, info = row$key)
    expect_identical(cov$se_theta, fit$se_theta, info = row$key)
  }
})

# In-fit SIR at 300 / 100 / seed 1, each under one non-default setting.
ps_sir_rows <- list(
  list(key = "w_sir_df3", example = "warfarin", settings = list(sir_df = 3)),
  list(key = "w_sir_maxit5", example = "warfarin", settings = list(inner_maxiter = 5L)),
  list(key = "iov_sir_natural", example = "warfarin_iov",
       settings = list(sir_scale = "natural")),
  # mm_oral, not warfarin_ode: since SIR re-solves each EBE as the fit did
  # (ferx-core #2279), ode_reltol no longer reaches warfarin_ode's exact linear
  # path, and dropping the record scored the same draws identically.
  list(key = "mm_sir_rtol3", example = "mm_oral",
       settings = list(ode_reltol = 1e-3)),
  # A multivariate-normal proposal: `.fitrx` writes its df as "normal" since
  # ferx-core #1819, so S5's reload repeats it too.
  list(key = "w_sir_normal", example = "warfarin", settings = list(sir_df = Inf))
)

ps_sir_fit <- function(row) {
  fit <- ps_fit(row$key, row$example, row$settings, sir = TRUE)
  skip_if(is.null(fit$sir_ess), paste(row$key, "in-fit SIR did not run"))
  fit
}

test_that("S4: ferx_sir(fit) with no arguments repeats the in-fit SIR", {
  for (row in ps_sir_rows) {
    fit <- ps_sir_fit(row)
    # Premise: the same draws without the records score differently. Both go:
    # without a SIR record, ferx_sir() scores under the fit's stage record,
    # which carries the scoring rows' setting too.
    no_record <- ps_without(fit, "sir_settings", "scoring_settings")
    expect_false(identical(
      ps_sir_fields(ferx_sir(no_record, 300L, 100L, sir_seed = 1L)),
      ps_sir_fields(fit)
    ), info = row$key)
    expect_identical(ps_sir_fields(ferx_sir(fit)), ps_sir_fields(fit), info = row$key)
  }
})

test_that("S3: the scoring record survives ferx_save_fit() / ferx_load_fit()", {
  for (row in ps_cov_rows) {
    fit <- ps_cov_fit(row)
    loaded <- ps_round_trip(fit)
    expect_identical(loaded$scoring_settings, fit$scoring_settings, info = row$key)
    # The bundle carries the packed estimate bit for bit (ferx-core #1815), so
    # the loaded fit repeats the in-fit step exactly.
    expect_identical(loaded$packed_estimate, fit$packed_estimate, info = row$key)
    cov <- ferx_covariance(loaded)
    expect_identical(cov$cov_matrix, fit$cov_matrix, info = row$key)
    expect_identical(cov$se_theta, fit$se_theta, info = row$key)
    # Premise: the record still moves the loaded fit's covariance.
    no_record <- ps_without(loaded, "scoring_settings")
    expect_false(identical(ferx_covariance(no_record)$cov_matrix, cov$cov_matrix),
                 info = row$key)
  }
})

ps_bundle_wire <- function(path) {
  con <- unz(path, "fit.json")
  on.exit(close(con))
  jsonlite::parse_json(paste(readLines(con, warn = FALSE), collapse = "\n"))
}

test_that("S5: the SIR record survives the round trip, so ferx_sir(loaded) repeats the in-fit SIR", {
  for (row in ps_sir_rows) {
    fit <- ps_sir_fit(row)
    path <- ps_scratch(".fitrx")
    unlink(path)
    ferx_save_fit(fit, path)
    # ferx-core refuses a bundle whose top-level `sir_seed` is not the record's.
    wire <- ps_bundle_wire(path)
    expect_identical(wire$sir_seed, 1L, info = row$key)
    expect_identical(wire$sir$settings$seed, 1L, info = row$key)

    loaded <- ferx_load_fit(path)
    expect_identical(loaded$sir_settings, fit$sir_settings, info = row$key)
    expect_identical(loaded$sir_seed_used, 1, info = row$key)
    expect_identical(ps_sir_fields(ferx_sir(loaded)), ps_sir_fields(fit), info = row$key)
    # Premise: without the records the loaded fit resamples differently.
    no_record <- ps_without(loaded, "sir_settings", "scoring_settings")
    expect_false(identical(
      ps_sir_fields(ferx_sir(no_record, 300L, 100L, sir_seed = 1L)),
      ps_sir_fields(fit)
    ), info = row$key)
  }
})

test_that("S10: the .fitrx spelling of both records is ferx-core's", {
  # ferx-core's ScoringSettingsWire keys, enums as their [fit_options] tokens.
  w <- fitrx_fit_json(ps_cov_fit(ps_cov_rows[[1L]])$bundle)$scoring_settings
  expect_setequal(names(w), ps_scoring_keys)
  expect_identical(w$inner_optimizer, "lbfgs")
  expect_null(w$ode_stiff_abort_after)

  # The covariate-NN penalty strengths (ferx-core #1352) are optional, and the
  # default record leaves them unrecorded: R edits it into a fit with no record,
  # and a recorded 0 would win over the model file's strength.
  sc <- ferx:::.ferx_default_settings("scoring_settings")
  expect_identical(sc$nn_l2_lambda, NA_real_)
  expect_identical(sc$nn_smooth_lambda, NA_real_)
  expect_identical(ferx:::.ferx_default_settings("sir_settings")$scoring$nn_l2_lambda,
                   NA_real_)

  # SirSettingsWire: the scoring half flattened beside the SIR keys.
  ws <- fitrx_fit_json(ps_sir_fit(ps_sir_rows[[3L]])$bundle)$sir$settings
  expect_setequal(names(ws), c("samples", "resamples", "seed", "df", "scale",
                               "keep_samples", ps_scoring_keys))
  expect_identical(ws$scale, "natural")
  # A normal proposal (`df = Inf`) is `"df": "normal"` on the wire, ferx-core's
  # own spelling (ferx-core#1819), and reads back as `Inf`.
  normal <- ps_sir_fit(ps_sir_rows[[5L]])
  expect_identical(fitrx_fit_json(normal$bundle)$sir$settings$df, "normal")
  expect_identical(normal$sir_settings$df, Inf)
})

test_that("S12: a record the engine would refuse is refused by name", {
  sc <- ferx:::.ferx_default_settings("scoring_settings")
  decode <- function(record) {
    ferx:::ferx_rust_settings_record("scoring_settings", record, "ferx_covariance",
                                     "fit$scoring_settings")
  }
  extra <- c(sc, list(inner_maxit = 5L))
  expect_error(decode(extra),
               'ferx_covariance: `fit$scoring_settings`: unknown key "inner_maxit"',
               fixed = TRUE)
  missing <- sc[names(sc) != "n_agq"]
  expect_error(decode(missing),
               'ferx_covariance: `fit$scoring_settings`: missing key "n_agq"',
               fixed = TRUE)
  bad <- sc
  bad$inner_optimizer <- "newton"
  expect_error(decode(bad),
               "ferx_covariance: `fit$scoring_settings`: fit option `inner_optimizer`: unknown value `newton`",
               fixed = TRUE)
  bad <- sc
  bad$ode_reltol <- -1
  expect_error(decode(bad),
               "ferx_covariance: `fit$scoring_settings`: ode_reltol must be a positive finite value, got -1",
               fixed = TRUE)
  # A bundle holding such a record does not load at all.
  expect_error(fitrx_edit_fit(ps_fit("w_default", "warfarin"), function(b) {
    b$scoring_settings$inner_optimizer <- "newton"
    b
  }), 'unknown scoring_settings.inner_optimizer "newton"', fixed = TRUE)
})

test_that("S6: an explicit argument overrides the fit's record, even at its default value", {
  # The engine resolves a record by value, so an explicit default passed as an
  # option would lose to it (plan section 0c); the argument edits the record.
  fit <- ps_cov_fit(ps_cov_rows[[3L]])  # warfarin_iov, recorded mu_referencing = FALSE
  no_record <- ps_without(fit, "scoring_settings")
  over <- ferx_covariance(fit, mu_referencing = TRUE)$cov_matrix
  expect_identical(over, ferx_covariance(no_record)$cov_matrix)
  expect_false(identical(over, fit$cov_matrix))
  expect_identical(ferx_covariance(fit, mu_referencing = FALSE)$cov_matrix,
                   fit$cov_matrix)

  nat <- ps_sir_fit(ps_sir_rows[[3L]])  # warfarin_iov, recorded sir_scale = "natural"
  packed <- ferx_sir(nat, sir_scale = "packed")
  ref <- ps_without(nat, "sir_settings")
  expect_identical(ps_sir_fields(packed),
                   ps_sir_fields(ferx_sir(ref, 300L, 100L, sir_seed = 1L)))
  expect_false(identical(ps_sir_fields(packed), ps_sir_fields(nat)))
  # The run records what it did; the settings not passed stay the fit's.
  expect_identical(packed$sir_settings$scale, "packed")
  expect_identical(packed$sir_settings[c("samples", "resamples", "seed")],
                   nat$sir_settings[c("samples", "resamples", "seed")])

  # A fit carrying only a SIR record (a bundle from between ferx-core #1758 and
  # #426) is scored under that record's scoring half; an explicit argument
  # edits that half, not the engine default (review round 1, finding 4).
  # inner_maxiter = 5 no longer binds the covariance step since ferx-core
  # #2280 (its near-mode inner BFGS converges sooner; the record moved the
  # matrix 8.8e-10 relative before, 0 after), so the record also carries a
  # loose inner_tol, which does (1.9e-8 relative).
  sir_only <- fitrx_edit_fit(ps_sir_fit(ps_sir_rows[[2L]]), function(w) {  # warfarin, inner_maxiter = 5
    w$scoring_settings <- NULL
    w$sir$settings$inner_tol <- 1e-2
    w
  })
  expect_identical(ferx_covariance(sir_only, mu_referencing = TRUE)$cov_matrix,
                   ferx_covariance(sir_only)$cov_matrix)
  # Premise: the SIR record's scoring half moves this fixture's covariance.
  neither <- ps_without(sir_only, "sir_settings")
  expect_false(identical(ferx_covariance(neither)$cov_matrix,
                         ferx_covariance(sir_only)$cov_matrix))
})

test_that("S13: a resamples > samples refusal names where each value came from", {
  fit <- ps_sir_fit(ps_sir_rows[[1L]])  # recorded 300 / 100
  expect_error(ferx_sir(fit, sir_samples = 50L),
               "`sir_resamples` (100, recorded on the fit) must be <= `sir_samples` (50)",
               fixed = TRUE)
  bare <- ps_without(fit, "sir_settings")
  expect_error(ferx_sir(bare, sir_samples = 100L),
               "`sir_resamples` (250, the default) must be <= `sir_samples` (100)",
               fixed = TRUE)
  expect_error(ferx_sir(bare, sir_samples = 10L, sir_resamples = 20L),
               "`sir_resamples` (20) must be <= `sir_samples` (10)", fixed = TRUE)
})

test_that("S7: a fit with no record is scored as it always was, and says nothing", {
  bare <- function(fit) {
    ps_without(fit, "scoring_settings", "sir_settings", "sir_seed", "packed_estimate")
  }
  f0 <- bare(ps_cov_fit(ps_cov_rows[[3L]]))
  expect_no_warning(cov <- ferx_covariance(f0))
  # An explicit argument on a record-less fit edits the engine's default record:
  # at the default value it is the same run.
  expect_identical(cov$cov_matrix,
                   ferx_covariance(f0, mu_referencing = TRUE)$cov_matrix)

  s0 <- bare(ps_sir_fit(ps_sir_rows[[1L]]))
  expect_no_warning(sir <- ferx_sir(s0))
  expect_identical(sir$sir_settings[c("samples", "resamples", "seed", "df", "scale")],
                   list(samples = 1000L, resamples = 250L, seed = 12345, df = 5,
                        scale = "packed"))
  expect_identical(ps_sir_fields(sir),
                   ps_sir_fields(ferx_sir(s0, 1000L, 250L, sir_scale = "packed")))
})

test_that("S8: a fit whose estimates were edited does not keep its packed estimate", {
  fit <- ps_fit("w_default", "warfarin")
  edited <- fitrx_edit_fit(fit, function(w) {
    w$theta$estimates[[1L]] <- w$theta$estimates[[1L]] * 1.01
    w
  })
  edited_np <- ps_without(edited, "packed_estimate")
  # With the packed vector the step would differentiate at the unedited point,
  # which is the in-fit covariance.
  cov <- ferx_covariance(edited)$cov_matrix
  expect_identical(cov, ferx_covariance(edited_np)$cov_matrix)
  expect_false(identical(cov, fit$cov_matrix))
})

test_that("S9: a [mixture] per-class override fit is scored in memory and after a reload", {
  ex <- list(model = test_path("fixtures", "mixture_override.ferx"),
             data = test_path("fixtures", "mixture_iv.csv"))
  fit <- ps_fit("mix_override", ex, method = "focei")
  # The fit records the fitted class-2 variance on its own (ferx-core #1765),
  # which the engine reads ahead of the packed estimate, so a fit without the
  # latter is still scored as the fit was.
  expect_length(fit$mixture_overrides, 1L)
  no_packed <- ps_without(fit, "packed_estimate")
  expect_identical(ferx_covariance(no_packed)$cov_matrix, fit$cov_matrix)
  # The descriptors and SEs go with the values: the loader refuses them alone
  # (ferx-core #2833).
  neither <- ps_without(no_packed, "mixture_overrides", "mixture_override_params",
                        "se_mixture_overrides")
  expect_error(ferx_covariance(neither), "records neither their fitted values",
               fixed = TRUE)

  cov <- ferx_covariance(fit)
  expect_identical(cov$cov_matrix, fit$cov_matrix)
  expect_identical(cov$se_theta, fit$se_theta)

  # The bundle carries both the packed estimate (ferx-core #1815) and the
  # overrides (#1765), bit for bit, so the reloaded fit is scored too.
  loaded <- ps_round_trip(fit)
  expect_identical(loaded$packed_estimate, fit$packed_estimate)
  expect_identical(loaded$mixture_overrides, fit$mixture_overrides)
  cov <- ferx_covariance(loaded)
  expect_identical(cov$cov_matrix, fit$cov_matrix)
  expect_identical(cov$se_theta, fit$se_theta)
})

test_that("S9b: an SAEM [mixture] override fit is scored in memory and after a reload", {
  # SAEM carries no packed estimate, so before ferx-core #1765 the post-hoc
  # steps refused this fit outright; its recorded overrides now let
  # ferx_covariance(fit) repeat the in-fit step.
  ex <- list(model = test_path("fixtures", "mixture_override.ferx"),
             data = test_path("fixtures", "mixture_iv.csv"))
  fit <- ps_fit("mix_override_saem", ex, method = "saem")
  expect_null(fit$packed_estimate)
  expect_length(fit$mixture_overrides, 1L)
  expect_identical(ferx_covariance(fit)$cov_matrix, fit$cov_matrix)
  loaded <- ps_round_trip(fit)
  expect_identical(loaded$mixture_overrides, fit$mixture_overrides)
  expect_identical(ferx_covariance(loaded)$cov_matrix, fit$cov_matrix)
})

# Since ferx-core #1824 `run_covariance` / `run_sir` score the LOQ method the
# fit reports, so the skeleton must carry the fit's (#558 review 1): a
# hard-coded "drop" re-scored an M3 fit without its censored rows (SE off by
# up to 37% on warfarin_bloq).
ps_bloq_drop_model <- function(env = parent.frame()) {
  # warfarin_bloq's model without its `bloq_method = m3` line: a drop-file
  # model, so M3 reaches the fit only through `ferx_fit(bloq_method =)`.
  path <- ps_scratch(".ferx", env)
  lines <- readLines(ferx_example("warfarin_bloq")$model)
  writeLines(lines[!grepl("^\\s*bloq_method\\s*=", lines)], path)
  path
}

test_that("S14: ferx_covariance(fit) on an M3 fit scores M3, in memory and after a reload", {
  file_m3 <- ps_fit("bloq_file_m3", "warfarin_bloq")
  drop_model <- ps_bloq_drop_model()
  arg_m3 <- ps_fit("bloq_arg_m3",
                   list(model = drop_model, data = ferx_example("warfarin_bloq")$data),
                   bloq_method = "m3")
  for (fit in list(file_m3, arg_m3)) {
    expect_identical(fit$bloq_method_label, "m3")
    # Premise: scored as drop, the censored rows move the covariance.
    as_drop <- fitrx_edit_fit(fit, function(w) {
      w$bloq_method <- "drop"
      w
    })
    expect_false(identical(ferx_covariance(as_drop)$se_theta, fit$se_theta))

    cov <- ferx_covariance(fit)
    expect_identical(cov$cov_matrix, fit$cov_matrix)
    expect_identical(cov$se_theta, fit$se_theta)
    loaded <- ps_round_trip(fit)
    expect_identical(loaded$bloq_method_label, "m3")
    expect_identical(ferx_covariance(loaded)$se_theta, fit$se_theta)
  }
})

test_that("S15: ferx_sir(fit) on an M3 fit repeats the in-fit SIR", {
  fit <- ps_fit("bloq_file_m3_sir", "warfarin_bloq", sir = TRUE)
  skip_if(is.null(fit$sir_ess), "bloq_file_m3_sir in-fit SIR did not run")
  expect_identical(ps_sir_fields(ferx_sir(fit)), ps_sir_fields(fit))
})

test_that("S13: inner_screen is recorded, and reads 0 when absent (#589)", {
  sc <- ferx:::.ferx_default_settings("scoring_settings")
  expect_identical(sc$inner_screen, 100L)
  mem <- sc[names(sc) != "inner_screen"]
  expect_identical(
    ferx:::ferx_rust_settings_record("scoring_settings", mem, "", "")$inner_screen, 0L
  )
})

test_that("S13b: inner_twin_seed is recorded, and reads FALSE when absent (ferx-core #2095)", {
  sc <- ferx:::.ferx_default_settings("scoring_settings")
  expect_identical(sc$inner_twin_seed, TRUE)
  mem <- sc[names(sc) != "inner_twin_seed"]
  expect_identical(
    ferx:::ferx_rust_settings_record("scoring_settings", mem, "", "")$inner_twin_seed, FALSE
  )
  # settings = list(inner_twin_seed =) reaches the engine's record.
  ex <- ferx_example("warfarin")
  fit_r <- function(settings) {
    suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                              verbose = FALSE, settings = settings))
  }
  expect_identical(fit_r(list(inner_twin_seed = FALSE))$scoring_settings$inner_twin_seed, FALSE)
  expect_identical(fit_r(list())$scoring_settings$inner_twin_seed, TRUE)
})

test_that("S14: settings = list(inner_screen =) reaches the engine as the model file's does (#589)", {
  ex <- ferx_example("warfarin")
  fit_r <- function(settings) {
    suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                              verbose = FALSE, settings = settings))
  }
  off <- fit_r(list(inner_screen = 0L))
  expect_identical(off$scoring_settings$inner_screen, 0L)
  expect_identical(fit_r(list(inner_screen = 7L))$scoring_settings$inner_screen, 7L)

  # The same fit with the option in [fit_options], as the CLI reads it.
  model <- ps_scratch(".ferx")
  text <- readLines(ex$model)
  at <- grep("^\\[fit_options\\]", text)
  writeLines(append(text, "  inner_screen = 0", after = at), model)
  from_file <- suppressWarnings(ferx_fit(model, ex$data, covariance = FALSE,
                                         verbose = FALSE))
  expect_identical(from_file$scoring_settings$inner_screen, 0L)
  expect_identical(from_file$ofv, off$ofv)
  expect_identical(from_file$theta, off$theta)
})

# Migration policy for the two records (#535). A key a record lacks is filled
# with what the engine that wrote it ran under, as ferx-core's
# `ScoringSettingsWire` fills it; a key this ferx does not know drops a
# bundle's record on load and is refused, with the way out, on a fit in memory.

# The keys a record may lack, and the value each is read as. Every other key
# is refused when missing. A key added to the records must be added here or
# to the refusals on purpose: the test walks the engine's own record.
ps_missing_fills <- list(inner_restarts = 1L, inner_screen = 0L, inner_twin_seed = FALSE,
                         gradient_method = "auto",
                         nn_l2_lambda = NA_real_, nn_smooth_lambda = NA_real_)

test_that("S16: a missing key is filled as ferx-core's wire fills it, or refused with the way out (#535)", {
  decode <- function(kind, record, field) {
    ferx:::ferx_rust_settings_record(kind, record, "ferx_covariance", field)
  }
  sc <- ferx:::.ferx_default_settings("scoring_settings")
  # A non-default value everywhere a fill could be confused with the default.
  sc$inner_restarts <- 4L
  sc$inner_screen <- 7L
  sc$inner_twin_seed <- TRUE
  sc$gradient_method <- "fd"
  sc$nn_l2_lambda <- 0.5
  sc$nn_smooth_lambda <- 0.25
  for (name in names(sc)) {
    rec <- sc[names(sc) != name]
    if (name %in% names(ps_missing_fills)) {
      back <- decode("scoring_settings", rec, "fit$scoring_settings")
      expect_identical(back[[name]], ps_missing_fills[[name]], info = name)
      expect_identical(back[names(back) != name], rec, info = name)
    } else {
      err <- tryCatch(decode("scoring_settings", rec, "fit$scoring_settings"),
                      error = conditionMessage)
      expect_match(err, sprintf('`fit$scoring_settings`: missing key "%s"', name),
                   fixed = TRUE, info = name)
      expect_match(err, "set `fit$scoring_settings <- NULL`", fixed = TRUE, info = name)
    }
  }
  # The scoring half of the SIR record is filled from the same table, and
  # every SIR key has been on the record since it existed.
  sr <- ferx:::.ferx_default_settings("sir_settings")
  sr$scoring[c("inner_restarts", "inner_screen")] <- NULL
  back <- decode("sir_settings", sr, "fit$sir_settings")$scoring
  expect_identical(back$inner_restarts, ps_missing_fills$inner_restarts)
  expect_identical(back$inner_screen, ps_missing_fills$inner_screen)
  for (name in setdiff(names(sr), "scoring")) {
    err <- tryCatch(decode("sir_settings", sr[names(sr) != name], "fit$sir_settings"),
                    error = conditionMessage)
    expect_match(err, sprintf('missing key "%s"', name), fixed = TRUE, info = name)
    expect_match(err, "set `fit$sir_settings <- NULL`", fixed = TRUE, info = name)
  }
})

test_that("S17: ferx_sir(fit, args) on a fit with no SIR record scores under its stage record", {
  skip_on_cran()
  # `run_sir` scores under the SIR record's scoring half before
  # fit$scoring_settings, so the record R builds for an explicit argument must
  # carry the fit's stage settings, not the engine default (#646 review: at
  # inner_maxiter 5 the default 200 moved the ESS from 119.406 to 119.575).
  fit <- ps_fit("w_maxit5_nosir", "warfarin", list(inner_maxiter = 5L))
  expect_null(fit$sir_settings)
  ran <- ferx_sir(fit, 300L, 100L, sir_seed = 1L)
  expect_identical(ran$sir_settings$scoring$inner_maxiter, 5L)
  # The record it wrote repeats the same draws.
  expect_identical(ps_sir_fields(ran), ps_sir_fields(ferx_sir(ran)))
  # Premise: under the engine default the draws score differently.
  bare <- ps_without(fit, "scoring_settings")
  expect_false(identical(ps_sir_fields(ferx_sir(bare, 300L, 100L, sir_seed = 1L)),
                         ps_sir_fields(ran)))
})

test_that("S18: a regularized [covariate_nn] fit records its NN penalty strength, and post-hoc records keep it", {
  skip_on_cran()
  # warfarin_dcm with `nn_l2 = 0.1`, evaluated at its initial estimates. A
  # recorded strength wins over the caller's in ferx-core (#1352), so every
  # record R hands a post-hoc step must carry the fit's, and a record R builds
  # for a fit with none must leave it unrecorded (NA) so the model file's
  # strength stands rather than 0.
  ex <- ferx_example("warfarin_dcm")
  t <- readLines(ex$model)
  t <- append(t, "  nn_l2 = 0.1", grep("^  method = focei", t))
  model <- ps_scratch(".ferx")
  writeLines(t, model)
  fit <- suppressWarnings(ferx_fit(model, ex$data, covariance = FALSE, verbose = FALSE,
                                   settings = list(maxiter = 0L)))
  expect_identical(fit$scoring_settings$nn_l2_lambda, 0.1)
  expect_identical(fit$scoring_settings$nn_smooth_lambda, 0)
  args <- ferx:::.ferx_fit_record_args(fit, scoring = list(mu_referencing = TRUE),
                                       sir = list(samples = 300L))
  expect_identical(args$scoring_settings$nn_l2_lambda, 0.1)
  expect_identical(args$sir_settings$scoring$nn_l2_lambda, 0.1)
  bare <- fit
  bare$scoring_settings <- NULL
  args <- ferx:::.ferx_fit_record_args(bare, scoring = list(mu_referencing = TRUE),
                                       sir = list(samples = 300L))
  expect_identical(args$scoring_settings$nn_l2_lambda, NA_real_)
  expect_identical(args$sir_settings$scoring$nn_l2_lambda, NA_real_)
  # And the strength survives the bundle.
  expect_identical(ps_round_trip(fit)$scoring_settings, fit$scoring_settings)
})

test_that("S18: a bundle record naming an unknown key is dropped with a warning; the rest loads (#535)", {
  fit <- ps_sir_fit(ps_sir_rows[[1L]])
  path <- ps_scratch(".fitrx")
  unlink(path)
  ferx_save_fit(fit, path, include_data = TRUE)
  newer <- list(
    scoring_settings = fitrx_edit_bundle(path, function(w) {
      w$scoring_settings$inner_future <- 7L
      w
    }),
    sir_settings = fitrx_edit_bundle(path, function(w) {
      w$sir$settings$inner_future <- 7L
      w
    })
  )
  without <- list(
    scoring_settings = fitrx_edit_bundle(path, function(w) {
      w$scoring_settings <- NULL
      w
    }),
    sir_settings = fitrx_edit_bundle(path, function(w) {
      w$sir$settings <- NULL
      w
    })
  )
  for (kind in names(newer)) {
    w <- expect_warning(loaded <- ferx_load_fit(newer[[kind]]),
                        class = "ferx_settings_record_dropped_warning")
    expect_match(conditionMessage(w), '"inner_future"', fixed = TRUE)
    expect_identical(w$keys, "inner_future")
    expect_identical(w$record, kind)
    expect_null(loaded[[kind]])
    other <- setdiff(names(newer), kind)
    expect_identical(loaded[[other]], fit[[other]], info = kind)
    expect_identical(loaded$theta, fit$theta, info = kind)
    # Exactly the fit a bundle without the record loads as (bar the paths each
    # load stages its bundled files to), so ferx_covariance() and ferx_sir()
    # run on it as on that one.
    staged <- c("bundle", "model_path", "data_path")
    expect_identical(loaded[setdiff(names(loaded), staged)],
                     ferx_load_fit(without[[kind]])[setdiff(names(loaded), staged)],
                     info = kind)
  }
  loaded <- suppressWarnings(ferx_load_fit(newer$scoring_settings))
  expect_identical(ferx_covariance(loaded)$cov_matrix,
                   ferx_covariance(ferx_load_fit(without$scoring_settings))$cov_matrix)
})
