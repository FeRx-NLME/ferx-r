# Unit tests for the pure-R diagnostic/formatting helpers in diagnostics.R:
# the warning-guidance lookup and the cli styling shims. No model fit required.

# ---------------------------------------------------------------------------
# .ferx_warning_guidance
# ---------------------------------------------------------------------------
test_that(".ferx_warning_guidance toggles the SDE hint for positive autocorrelation", {
  pos_ode <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive", uses_sde = FALSE)
  expect_match(pos_ode, "Positive IWRES autocorrelation")
  expect_match(pos_ode, "SDE process noise")            # hint shown when not already using SDE
  pos_sde <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive", uses_sde = TRUE)
  expect_match(pos_sde, "Positive IWRES autocorrelation")
  expect_false(grepl("SDE process noise", pos_sde))     # hint suppressed
})

test_that(".ferx_warning_guidance returns negative-autocorrelation guidance", {
  neg <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "Negative autocorrelation")
  expect_match(neg, "Negative IWRES autocorrelation")
})

test_that(".ferx_warning_guidance maps every known category to non-empty text", {
  cats <- c("convergence", "covariance_step", "condition_number", "optimizer_health",
            "eta_normality", "bloq_method", "sir", "importance_sampling", "data_quality",
            "omega_structure", "ebe_convergence", "gradient_fallback", "mu_referencing",
            "optimizer_config", "multi_start", "threads", "cancelled", "unused_parameter")
  for (cat in cats) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
})

test_that(".ferx_warning_guidance returns NULL for an unknown category", {
  expect_null(ferx:::.ferx_warning_guidance("not_a_real_category"))
})

# ---------------------------------------------------------------------------
# .ferx_use_cli / .ferx_style
# ---------------------------------------------------------------------------
test_that(".ferx_use_cli follows cli's colour detection", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 1)
  expect_false(ferx:::.ferx_use_cli())
  withr::local_options(cli.num_colors = 256)
  expect_true(ferx:::.ferx_use_cli())
})

test_that(".ferx_style returns text unchanged when cli is off", {
  for (s in c("bold", "green", "red", "yellow", "dim", "unknownstyle")) {
    expect_identical(ferx:::.ferx_style("hi", s, use_cli = FALSE), "hi")
  }
})

test_that(".ferx_style applies ANSI for each known style when cli is on", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 256)
  for (s in c("bold", "green", "red", "yellow", "dim")) {
    styled <- ferx:::.ferx_style("hi", s, use_cli = TRUE)
    expect_true(cli::ansi_has_any(styled), info = s)   # colour actually applied
    expect_identical(cli::ansi_strip(styled), "hi", info = s) # text intact
  }
  # Unknown style falls through to the raw text even with cli on.
  expect_identical(ferx:::.ferx_style("hi", "nope", use_cli = TRUE), "hi")
})
