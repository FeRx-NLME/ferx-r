# Helpers ---------------------------------------------------------------

make_fake_fit <- function(...) {
  defaults <- list(
    model_name = "test_model",
    data_name = "test_data",
    converged = TRUE,
    method = "focei",
    method_chain = NULL,
    gradient = "auto",
    ofv = -100.0,
    aic = -90.0,
    bic = -80.0,
    n_subjects = 10L,
    n_obs = 100L,
    n_parameters = 5L,
    n_iterations = 50L,
    theta = c(CL = 1.0, V = 10.0),
    se_theta = NULL,
    omega = NULL,
    se_omega = NULL,
    sigma = NULL,
    se_sigma = NULL,
    shrinkage_eta = NULL,
    shrinkage_eps = NULL,
    covariance_status = "not_requested",
    wall_time_secs = 5.0,
    ferx_version = "0.1.0",
    ebe_convergence_warnings = 0L,
    max_unconverged_subjects = NULL,
    total_ebe_fallbacks = 0L,
    call_settings = list(),
    sir_ess = NULL,
    saem_n_subjects_hmc = NULL,
    warnings = character(0)
  )
  structure(modifyList(defaults, list(...)), class = "ferx_fit")
}

# summary.ferx_fit structure -------------------------------------------

test_that("summary.ferx_fit returns a ferx_summary with expected fields", {
  fit <- make_fake_fit()
  s <- summary(fit)

  expect_s3_class(s, "ferx_summary")
  expect_named(s, c(
    "model_name", "data_name", "gradient", "gradient_used",
    "method", "method_chain",
    "converged", "ofv", "aic", "bic", "n_subjects", "n_obs",
    "n_parameters", "n_iterations", "theta", "se_theta", "omega",
    "se_omega", "sigma", "se_sigma", "shrinkage_eta", "shrinkage_eps",
    "covariance_status", "eigenvalues", "condition_number",
    "wall_time_secs", "ferx_version",
    "ebe_convergence_warnings", "max_unconverged_subjects",
    "total_ebe_fallbacks", "model_structure", "call_settings",
    "model_file_settings", "sir_ess", "importance_sampling", "warnings",
    "uses_sde", "dw_statistic", "iwres_lag1_r", "saem_n_subjects_hmc"
  ), ignore.order = TRUE)
})

test_that("summary.ferx_fit passes through call_settings", {
  fit <- make_fake_fit(call_settings = list(optimizer = "slsqp", max_iter = 200L))
  s <- summary(fit)

  expect_equal(s$call_settings$optimizer, "slsqp")
  expect_equal(s$call_settings$max_iter, 200L)
})

test_that("summary.ferx_fit defaults call_settings to empty list when absent", {
  fit <- make_fake_fit()
  fit$call_settings <- NULL  # simulate pre-PR fit object
  s <- summary(fit)

  expect_equal(s$call_settings, list())
})

test_that("summary.ferx_fit passes through warnings", {
  fit <- make_fake_fit(warnings = c("ETA1 may be non-normal", "ETA2 boundary"))
  s <- summary(fit)

  expect_equal(s$warnings, c("ETA1 may be non-normal", "ETA2 boundary"))
})

test_that("summary.ferx_fit passes through sir_ess", {
  fit <- make_fake_fit(sir_ess = 847.2)
  s <- summary(fit)

  expect_equal(s$sir_ess, 847.2)
})

test_that("summary.ferx_fit passes through method_chain", {
  fit <- make_fake_fit(method = "focei", method_chain = c("saem", "focei"))
  s <- summary(fit)

  expect_equal(s$method_chain, c("saem", "focei"))
})

# print.ferx_summary output --------------------------------------------

test_that("print.ferx_summary shows method in uppercase for a single method", {
  s <- summary(make_fake_fit(method = "focei", method_chain = NULL))
  out <- capture.output(print(s))

  expect_true(any(grepl("Method:\\s+FOCEI", out)))
})

test_that("print.ferx_summary shows full chain in uppercase when method_chain > 1", {
  s <- summary(make_fake_fit(method = "focei", method_chain = c("saem", "focei")))
  out <- capture.output(print(s))

  expect_true(any(grepl("Method:\\s+SAEM -> FOCEI", out)))
  expect_false(any(grepl("saem|focei", out)))
})

test_that("print.ferx_summary omits Settings block when both settings sources are empty", {
  s <- summary(make_fake_fit(call_settings = list()))
  out <- capture.output(print(s))

  expect_false(any(grepl("^Settings", out)))
})

test_that("print.ferx_summary shows Settings block when call_settings is non-empty", {
  s <- summary(make_fake_fit(call_settings = list(optimizer = "slsqp", max_iter = 200)))
  out <- capture.output(print(s))

  expect_true(any(grepl("^Settings", out)))
  expect_true(any(grepl("optimizer", out)))
  expect_true(any(grepl("slsqp", out)))
  expect_true(any(grepl("max_iter", out)))
})

test_that("print.ferx_summary omits Warnings block when no warnings", {
  s <- summary(make_fake_fit(
    warnings = character(0), ebe_convergence_warnings = 0L, total_ebe_fallbacks = 0L
  ))
  out <- capture.output(print(s))

  expect_false(any(grepl("^Warnings:", out)))
})

test_that("print.ferx_summary shows model warnings", {
  s <- summary(make_fake_fit(warnings = "ETA1 Shapiro-Wilk p=0.03"))
  out <- capture.output(print(s))

  expect_true(any(grepl("Warnings:", out)))
  expect_true(any(grepl("ETA1 Shapiro-Wilk", out)))
})

test_that("print.ferx_summary aggregates EBE convergence warnings", {
  s <- summary(make_fake_fit(
    ebe_convergence_warnings = 3L,
    max_unconverged_subjects = 2L
  ))
  out <- capture.output(print(s))

  expect_true(any(grepl("Warnings:", out)))
  expect_true(any(grepl("3 EBE convergence", out)))
  expect_true(any(grepl("max unconverged subjects: 2", out)))
})

test_that("print.ferx_summary aggregates EBE fallback warnings", {
  s <- summary(make_fake_fit(total_ebe_fallbacks = 5L))
  out <- capture.output(print(s))

  expect_true(any(grepl("5 EBE fallback", out)))
})

test_that("print.ferx_summary omits SIR ESS line when sir_ess is NULL", {
  s <- summary(make_fake_fit(sir_ess = NULL))
  out <- capture.output(print(s))

  expect_false(any(grepl("SIR ESS", out)))
})

test_that("print.ferx_summary shows SIR ESS when sir_ess is present", {
  s <- summary(make_fake_fit(sir_ess = 512.3))
  out <- capture.output(print(s))

  expect_true(any(grepl("SIR ESS.*512\\.3", out)))
})

test_that("print.ferx_summary shows R package version with core version in parens", {
  s <- summary(make_fake_fit(ferx_version = "0.1.0"))
  out <- capture.output(print(s))
  pkg_v <- as.character(utils::packageVersion("ferx"))
  expect_true(any(grepl(
    paste0("ferx v", pkg_v, " \\(core v0\\.1\\.0\\)"), out, fixed = FALSE
  )))
})

test_that("print.ferx_fit shows R package version with core version in parens", {
  om <- matrix(0.04, 1, 1, dimnames = list("ETA_CL", "ETA_CL"))
  fit <- make_fake_fit(
    ferx_version = "0.1.0",
    omega = om,
    sigma = 0.01
  )
  out <- capture.output(suppressWarnings(print(fit)))
  pkg_v <- as.character(utils::packageVersion("ferx"))
  expect_true(any(grepl(
    paste0("ferx v", pkg_v, " \\(core v0\\.1\\.0\\)"), out, fixed = FALSE
  )))
})


test_that("print.ferx_summary wraps output in dashed borders", {
  s <- summary(make_fake_fit())
  out <- capture.output(print(s))

  bars <- grep("^-{10,}$", out)
  expect_gte(length(bars), 2L)
})

test_that("print.ferx_summary shows Structure line when model_structure is non-NULL", {
  s <- summary(make_fake_fit(model_structure = list(
    theta_names = c("TVCL", "TVV"),
    model_type  = "1-cpt oral",
    iiv         = c("ETA_CL", "ETA_V"),
    iov         = character(0),
    residual    = "proportional"
  )))
  out <- capture.output(print(s))

  expect_true(any(grepl("^Structure:", out)))
  expect_true(any(grepl("1-cpt oral", out)))
  expect_true(any(grepl("ETA_CL", out)))
  expect_true(any(grepl("proportional", out)))
})

test_that("print.ferx_summary omits Structure line when model_structure is NULL", {
  s <- summary(make_fake_fit())
  s$model_structure <- NULL
  out <- capture.output(print(s))

  expect_false(any(grepl("^Structure:", out)))
})

# gradient_used display ------------------------------------------------

test_that(".ferx_short_gradient_label maps engine labels to short tokens", {
  expect_identical(ferx:::.ferx_short_gradient_label("Enzyme AD"), "ad")
  expect_identical(ferx:::.ferx_short_gradient_label("finite differences"), "fd")
  expect_identical(ferx:::.ferx_short_gradient_label("N/A"), "N/A")
  expect_identical(ferx:::.ferx_short_gradient_label(""), NA_character_)
  expect_identical(ferx:::.ferx_short_gradient_label(NULL), NA_character_)
})

test_that("summary.ferx_fit carries gradient_used through", {
  fit <- make_fake_fit(gradient = "auto", gradient_used = "ad")
  s <- summary(fit)
  expect_identical(s$gradient_used, "ad")
})

test_that("print.ferx_summary shows requested -> used arrow when used is known", {
  s <- summary(make_fake_fit(gradient = "auto", gradient_used = "ad"))
  out <- capture.output(print(s))
  expect_true(any(grepl("Gradient:\\s+auto \\(requested\\) -> ad \\(used\\)", out)))
})

test_that("print.ferx_summary falls back to (requested) only when used is NA", {
  s <- summary(make_fake_fit(gradient = "auto", gradient_used = NA_character_))
  out <- capture.output(print(s))
  expect_true(any(grepl("Gradient:\\s+auto \\(requested\\)$", out)))
  expect_false(any(grepl("\\(used\\)", out)))
})

# .ferx_parse_setting_value --------------------------------------------

test_that(".ferx_parse_setting_value converts booleans", {
  expect_identical(ferx:::.ferx_parse_setting_value("true"), TRUE)
  expect_identical(ferx:::.ferx_parse_setting_value("false"), FALSE)
})

test_that(".ferx_parse_setting_value converts numerics", {
  expect_equal(ferx:::.ferx_parse_setting_value("200"), 200)
  expect_equal(ferx:::.ferx_parse_setting_value("3.14"), 3.14)
})

test_that(".ferx_parse_setting_value converts null to NA", {
  expect_identical(ferx:::.ferx_parse_setting_value("null"), NA)
})

test_that(".ferx_parse_setting_value leaves strings as character", {
  expect_identical(ferx:::.ferx_parse_setting_value("slsqp"), "slsqp")
})

# call_settings round-trip type fidelity --------------------------------

test_that("call_settings stores typed values after settings are serialised", {
  # Simulate what ferx_fit() does: serialise then parse back.
  parts <- ferx:::.ferx_settings_to_strings(list(
    optimizer = "slsqp",
    max_iter = 200L,
    scale_params = TRUE,
    stagnation_guard = FALSE,
    reconverge_gradient_interval = 5L
  ))
  typed <- setNames(
    lapply(parts$values, ferx:::.ferx_parse_setting_value),
    parts$keys
  )

  expect_identical(typed$optimizer, "slsqp")
  expect_equal(typed$max_iter, 200)
  expect_identical(typed$scale_params, TRUE)
  expect_identical(typed$stagnation_guard, FALSE)
  expect_equal(typed$reconverge_gradient_interval, 5)
})

test_that("reconverge_gradient_interval serialises as a plain integer string", {
  # Integer-valued settings must not pick up scientific notation or decimals
  # when forwarded to the Rust fit-option parser (which expects a usize).
  parts <- ferx:::.ferx_settings_to_strings(list(reconverge_gradient_interval = 10L))
  expect_identical(parts$keys, "reconverge_gradient_interval")
  expect_identical(parts$values, "10")

  # 0 is the documented "off" value and must round-trip, not become "" or NA.
  parts0 <- ferx:::.ferx_settings_to_strings(list(reconverge_gradient_interval = 0L))
  expect_identical(parts0$values, "0")
})
