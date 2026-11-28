# A [mixture] fit that records its per-class override values only in
# `fit$packed_estimate` (#675): a `.fitrx` saved before #655 carries no
# `mixture_overrides`. ferx_sir() / ferx_covariance() read the packed estimate
# from the bytes (test-posthoc-settings.R S9), and so does every from-fit
# reader that rebuilds the classes, instead of refusing the fit with
# `E_FIT_INPUT_MISSING`.
#
# The oracle is the full fit's own output: the engine takes the packed
# estimate only when it unpacks bit for bit to the fit's theta / omega / sigma,
# and its override segment holds the values `mixture_overrides` does, so the
# same seed gives identical output.

mp_model <- function() test_path("fixtures", "mixture_override.ferx")
mp_data <- function() test_path("fixtures", "mixture_iv.csv")
# The fixture fitted with its covariance step, for the uncertainty paths.
mp_cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- muffle_mixl_lint(ferx_fit(mp_model(), mp_data(), method = "focei",
                                        covariance = TRUE, verbose = FALSE))
    }
    fit
  }
})

# `run(fit)` and `run()` on the fit without its override values are identical.
expect_same_without_override_values <- function(fit, run) {
  legacy <- mixture_fit_without_override_values(fit)
  # Premise: the values are gone and the packed estimate is not.
  expect_length(fit$mixture_overrides, 1L)
  expect_null(legacy$mixture_overrides)
  expect_gt(length(legacy$packed_estimate), 0L)
  expect_identical(legacy$packed_estimate, fit$packed_estimate)
  full <- muffle_mixl_lint(run(fit))
  expect_identical(muffle_mixl_lint(run(legacy)), full)
  invisible(full)
}

test_that("ferx_simulate(fit = ) reads the overrides from the packed estimate (#675)", {
  full <- expect_same_without_override_values(mixture_override_fit(), function(f) {
    ferx_simulate(mp_model(), mp_data(), n_sim = 20L, seed = 5L, fit = f)
  })
  expect_setequal(unique(full$MIXNUM), 1:2)
})

test_that("ferx_calc_npde() reads the overrides from the packed estimate (#675)", {
  full <- expect_same_without_override_values(mixture_override_fit(), function(f) {
    ferx_calc_npde(f, nsim = 100L, seed = 4L, model = mp_model(), data = mp_data())$sdtab
  })
  expect_true(all(is.finite(full$NPD)))
})

test_that("ferx_simulate_adaptive(fit = ) reads the overrides from the packed estimate (#675)", {
  # Same parameters and bounds as the fit's model, so the same packed layout.
  path <- .write_mixture_adaptive_model()
  full <- expect_same_without_override_values(mixture_override_fit(), function(f) {
    ferx_simulate_adaptive(path, mp_data(), n_sim = 3L, seed = 1L, fit = f)
  })
  expect_gt(nrow(full$trajectories), 0L)
  # With `uncertainty =` the draws are taken around a second rebuild of the
  # fit, which needs the packed estimate as well (#803 review, finding 1).
  fit <- mp_cov_fit()
  skip_if(is.null(fit$cov_matrix), "mixture covariance step did not converge")
  expect_same_without_override_values(fit, function(f) {
    ferx_simulate_adaptive(path, mp_data(), n_sim = 2L, seed = 1L, fit = f,
                           uncertainty = list(n_draws = 2L))
  })
})

test_that("ferx_simulate_with_uncertainty() reads the overrides from the packed estimate (#675)", {
  fit <- mp_cov_fit()
  skip_if(is.null(fit$cov_matrix), "mixture covariance step did not converge")
  full <- expect_same_without_override_values(fit, function(f) {
    ferx_simulate_with_uncertainty(mp_model(), mp_data(), f, n_uncertainty_draws = 3L,
                                   n_sim_per_draw = 10L, seed = 3L)
  })
  expect_setequal(unique(full$MIXNUM), 1:2)
})

# The four TTE readers, on the fixture with two competing TTE endpoints
# (cause A's hazard carries ETA_CL). The VPCs simulate their replicates, so
# they rebuild the classes and refused such a fit; ferx_predict_survival()
# reads no class matrix and must not change a row (#803 review).
mp_tte_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- muffle_mixl_lint(ferx_fit(
        test_path("fixtures", "mixture_tte.ferx"), test_path("fixtures", "mixture_tte.csv"),
        method = "focei", verbose = FALSE
      ))
    }
    fit
  }
})

test_that("the TTE VPCs and ferx_predict_survival() read the overrides from the packed estimate (#675)", {
  mod <- test_path("fixtures", "mixture_tte.ferx")
  dat <- test_path("fixtures", "mixture_tte.csv")
  fit <- mp_tte_fit()
  expect_true(fit$converged)
  expect_same_without_override_values(fit, function(f) {
    ferx_vpc_tte(mod, dat, n_sim = 20L, seed = 1L, fit = f, cmt = 2L)
  })
  expect_same_without_override_values(fit, function(f) {
    ferx_vpc_kmmc(mod, dat, covariate = "WT", n_sim = 20L, seed = 1L, fit = f, cmt = 2L)
  })
  expect_same_without_override_values(fit, function(f) {
    ferx_vpc_cif(mod, dat, n_sim = 20L, seed = 1L, fit = f)
  })
  expect_same_without_override_values(fit, function(f) {
    ferx_predict_survival(mod, dat, times = c(0, 6, 12), fit = f)
  })
})

test_that("ferx_predict(fit = ) predicts as the full fit (#675)", {
  # PRED reads no class Omega/Sigma, so this path never refused such a fit; the
  # packed estimate now reaches it too and must not change a row.
  expect_same_without_override_values(mixture_override_fit(), function(f) {
    ferx_predict(mp_model(), mp_data(), fit = f)
  })
})

test_that("a [mixture] fit with neither override values nor a packed estimate is refused, coded (#675)", {
  neither <- fitrx_edit_fit(mixture_fit_without_override_values(mixture_override_fit()),
                            function(w) {
                              w$packed_estimate <- NULL
                              w
                            })
  expect_null(neither$packed_estimate)
  probe <- engine_error_probe(muffle_mixl_lint(ferx_simulate(
    mp_model(), mp_data(), n_sim = 1L, seed = 1L, fit = neither
  )))
  # Coded, as the same refusal from ferx_simulate_with_uncertainty() is (#672).
  expect_coded_refusal(probe, "omega(2) ETA_CL", "E_FIT_INPUT_MISSING")
  expect_null(probe$value)
})

test_that("every glue binding that takes mixture_overrides takes packed_estimate too (#675)", {
  # The two travel together from validate_fit_for_params(): a new from-fit
  # binding given only the first would refuse an old [mixture] fit again.
  ns <- asNamespace("ferx")
  takes <- function(arg) {
    Filter(function(f) arg %in% names(formals(get(f, ns))),
           grep("^ferx_rust_", ls(ns), value = TRUE))
  }
  readers <- takes("mixture_overrides")
  expect_gte(length(readers), 10L)
  expect_setequal(takes("packed_estimate"), readers)
  # And the R side hands both to each of them.
  pieces <- validate_fit_for_params(mixture_override_fit())
  expect_identical(pieces$packed_estimate, mixture_override_fit()$packed_estimate)
})
