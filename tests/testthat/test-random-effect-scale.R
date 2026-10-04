# A random effect's variance is reported on the scale the effect enters its
# parameter: CV% for a log-normal effect, SD for an additive or logit one,
# nothing for any other shape (ferx-core #1643). The ETA side reads each ETA's
# `eta_param_info` entry *by name* (#438): the engine lists that in
# [individual_parameters] statement order, not in `eta_names` order. The kappa
# side reads `kappa_param_types`, which the engine keeps parallel to
# `kappa_names`.

ferx_rust_eta_info_by_name <- getFromNamespace("ferx_rust_eta_info_by_name", "ferx")

# One mbma_placebo fit per file: KAPPA_ARM is the additive, NARM-weighted kappa
# (`BASE = TVE0 + PLACEBO + ETA_E0 + KAPPA_ARM`) that read as CV% = 1249.1.
.mbma_fit_cache <- new.env(parent = emptyenv())
mbma_fit <- function() {
  if (is.null(.mbma_fit_cache$fit)) {
    ex <- ferx_example("mbma_placebo")
    .mbma_fit_cache$fit <- ferx_fit(ex$model, ex$data, verbose = FALSE)
  }
  .mbma_fit_cache$fit
}

# #438's fixture: ETAs declared ETA_V, ETA_CL, statements in the other order,
# and one ETA (ETA_CL) sharing its statement with a kappa.
order_fixture_fit <- function(env = parent.frame()) {
  model <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  theta TVV(8.1, 0.1, 500.0)",
    "  theta TVKA(1.0, 0.01, 50.0)",
    "  omega ETA_V  ~ 1.0",
    "  omega ETA_CL ~ 0.07",
    "  kappa KAPPA_CL ~ 0.04",
    "  sigma PROP_ERR ~ 0.01 (sd)",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL + KAPPA_CL)",
    "  V  = TVV + ETA_V",
    "  KA = TVKA",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)",
    "",
    "[fit_options]",
    "  method     = foce",
    "  iov_column = OCC",
    "  covariance = false"
  ), model)
  ferx_fit(model, ferx_example("warfarin_iov")$data, verbose = FALSE,
           settings = list(maxiter = 5L))
}

# The printed row of one parameter.
printed_row <- function(fit, pattern) {
  out <- capture.output(print(fit))
  row <- grep(pattern, out, value = TRUE)
  expect_length(row, 1L)
  row
}

# Rewrite one bundle's fit.json in place.
edit_bundle <- function(path, edit) {
  staging <- withr::local_tempdir()
  utils::unzip(path, exdir = staging)
  json <- file.path(staging, "fit.json")
  wire <- jsonlite::read_json(json, simplifyVector = FALSE)
  jsonlite::write_json(edit(wire), json, auto_unbox = TRUE, digits = NA,
                       null = "null", na = "null")
  unlink(path)
  withr::with_dir(staging, utils::zip(path, list.files(staging), flags = "-q"))
  invisible(path)
}

# -- kappa rows (ferx-core #1643) ---------------------------------------------

test_that("mbma_placebo: KAPPA_ARM prints an SD at weight 1, not a CV%", {
  skip_on_cran()
  fit <- mbma_fit()
  expect_identical(fit$kappa_param_types, c(KAPPA_ARM = "additive"))
  row <- printed_row(fit, "^  KAPPA_ARM = ")
  var <- fit$omega_iov[1L, 1L]
  expect_true(grepl(sprintf("(SD = %.4f at weight 1)", sqrt(var)), row, fixed = TRUE),
              info = row)
  expect_false(grepl("CV%", row, fixed = TRUE), info = row)
})

test_that("mbma_placebo: ferx_estimates() gives KAPPA_ARM's scale as additive", {
  skip_on_cran()
  est <- mbma_fit()$estimates
  expect_identical(est["KAPPA_ARM", "scale"], "additive")
  # An omega row carries its ETA's engine classification.
  expect_identical(est["ETA_E0", "scale"], mbma_fit()$eta_param_types[[1L]])
  # The scale is a separate column: the estimate is still a variance, and its
  # interval is still the variance-scale Wald interval.
  expect_identical(est["KAPPA_ARM", "transform"], "variance")
  expect_true(is.na(est["TVE0", "scale"]))
  expect_true(is.na(est["ADD_ERR", "scale"]))
})

test_that("a kappa row reads on its scale, both sides of every gate", {
  note <- getFromNamespace(".ferx_kappa_variance_note", "ferx")
  # log-normal and unknown: the exact CV% every row printed before.
  expect_identical(note("log_normal", 0.2, FALSE), "CV% = 47.1")
  expect_identical(note("log_normal", 0.2, TRUE), "CV% = 47.1")
  expect_identical(note(NA_character_, 0.2, FALSE), "CV% = 47.1")
  expect_identical(note("log_normal", -1, FALSE), "CV% = 0.0")
  # additive: an SD, labelled at weight 1 only when the kappa is weighted.
  expect_identical(note("additive", 0.25, FALSE), "SD = 0.5000")
  expect_identical(note("additive", 0.25, TRUE), "SD = 0.5000 at weight 1")
  expect_identical(note("additive", -1, FALSE), "SD = 0.0000")
  # logit: an SD on the logit scale.
  expect_identical(note("logit", 0.25, FALSE), "SD = 0.5000, logit scale")
  expect_identical(note("logit", 0.25, TRUE), "SD = 0.5000 at weight 1, logit scale")
  # custom: nothing to say.
  expect_null(note("custom", 0.25, FALSE))
  expect_null(note("custom", 0.25, TRUE))
})

test_that("a custom kappa prints no parenthetical at all", {
  fit <- make_fake_fit(
    omega = matrix(0.10, 1, 1),
    omega_iov = matrix(0.04, 1, 1), kappa_names = "KAPPA_V",
    se_kappa = 0.01, shrinkage_kappa = 0.1, kappa_param_types = "custom"
  )
  row <- printed_row(fit, "^  KAPPA_V = ")
  expect_identical(row, "  KAPPA_V = 0.040000  SE = 0.010000  Shrinkage = 10.0%")
})

test_that("print labels an additive kappa's SD 'at weight 1' only when it is weighted", {
  base <- list(
    omega = matrix(0.10, 1, 1),
    omega_iov = matrix(0.04, 1, 1), kappa_names = "KAPPA_V",
    se_kappa = 0.01, shrinkage_kappa = 0.1, kappa_param_types = "additive"
  )
  plain <- do.call(make_fake_fit, base)
  expect_identical(printed_row(plain, "^  KAPPA_V = "),
                   "  KAPPA_V = 0.040000  (SD = 0.2000)  SE = 0.010000  Shrinkage = 10.0%")
  weighted <- do.call(make_fake_fit, c(base, list(
    kappa_weights = c(KAPPA_V = "NARM"), kappa_weight_typical = c(KAPPA_V = 4)
  )))
  expect_identical(printed_row(weighted, "^  KAPPA_V = "),
                   "  KAPPA_V = 0.040000  (SD = 0.2000 at weight 1)  SE = 0.010000  Shrinkage = 10.0%")
})

# -- ETA rows by name (#438) --------------------------------------------------

test_that("eta types follow the ETA's name, not eta_param_info's position", {
  skip_on_cran()
  fit <- order_fixture_fit()
  expect_identical(fit$eta_names, c("ETA_V", "ETA_CL"))
  expect_identical(fit$eta_param_types, c("additive", "log_normal"))
  expect_identical(fit$eta_linked_theta, c("TVV", "TVCL"))
  expect_true(grepl("[additive]", printed_row(fit, "^  ETA_V "), fixed = TRUE))
  expect_true(grepl("[log-normal]", printed_row(fit, "^  ETA_CL "), fixed = TRUE))
  # The kappa shares ETA_CL's statement, so it is log-normal too.
  expect_identical(fit$kappa_param_types, c(KAPPA_CL = "log_normal"))
  expect_identical(fit$estimates["ETA_V", "scale"], "additive")
  expect_identical(fit$estimates["ETA_CL", "scale"], "log_normal")
})

test_that("the by-name lookup handles a skipped, a repeated and a disagreeing ETA", {
  out <- ferx_rust_eta_info_by_name(
    c("ETA_A", "ETA_B", "ETA_C", "ETA_D"),
    c("ETA_C", "ETA_A", "ETA_C", "ETA_D", "ETA_D"),
    c("additive", "logit", "additive", "log_normal", "additive"),
    c("TVC", "TVA", "TVC2", "TVD", "TVD")
  )
  # ETA_A: one entry. ETA_B: none, so the log-normal default and no theta.
  # ETA_C: twice, types agree, thetas do not. ETA_D: twice, types disagree.
  expect_identical(out$eta_param_types, c("logit", "log_normal", "additive", "custom"))
  expect_identical(out$eta_linked_theta, c("TVA", "", "", "TVD"))
  # No classification at all (an old engine, a skeleton result): empty, so the
  # R-side default applies exactly as before.
  none <- ferx_rust_eta_info_by_name(c("ETA_A"), character(), character(), character())
  expect_identical(none$eta_param_types, character())
  expect_identical(none$eta_linked_theta, character())
})

test_that("a bundle with statement-ordered eta_param_info loads by name", {
  skip_on_cran()
  fit <- order_fixture_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  # What ferx-core itself writes: statement order, ETA_CL first.
  edit_bundle(path, function(w) {
    w$eta_param_info <- list(
      list(eta_name = "ETA_CL", param_type = "log_normal", linked_theta = "TVCL",
           individual_param_name = "CL"),
      list(eta_name = "ETA_V", param_type = "additive", linked_theta = "TVV",
           individual_param_name = "V")
    )
    w
  })
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$eta_param_types, c("additive", "log_normal"))
  expect_identical(loaded$eta_linked_theta, c("TVV", "TVCL"))
  expect_identical(unname(loaded$kappa_param_types), "log_normal")
})

# -- persistence ----------------------------------------------------------------

test_that("kappa_param_types survives a save/load round trip", {
  skip_on_cran()
  fit <- mbma_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$kappa_param_types, fit$kappa_param_types)
  expect_identical(printed_row(loaded, "^  KAPPA_ARM = "),
                   printed_row(fit, "^  KAPPA_ARM = "))
})

test_that("a .fitrx saved before kappa_param_types loads and prints as before", {
  skip_on_cran()
  fit <- mbma_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  # The writer before this change never emitted the key.
  edit_bundle(path, function(w) {
    w$iov$kappa_param_types <- NULL
    w
  })
  loaded <- ferx_load_fit(path)
  expect_length(loaded$kappa_param_types, 0L)
  # The reloaded variance, not the fit's: the JSON round trip can move its
  # last digit, and this CV% (of an additive kappa) is ~1e35.
  var <- loaded$omega_iov[1L, 1L]
  row <- printed_row(loaded, "^  KAPPA_ARM = ")
  # The row exactly as print.ferx_fit() wrote it before: the log-normal CV%.
  expect_true(startsWith(row, sprintf(
    "  KAPPA_ARM = %.6f  (CV%% = %.1f)  SE = ", var, sqrt(exp(var) - 1) * 100
  )), info = row)
  expect_true(is.na(loaded$estimates["KAPPA_ARM", "scale"]))
  expect_identical(loaded$estimates["KAPPA_ARM", "transform"], "variance")
})
