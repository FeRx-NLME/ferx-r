# Tests for vine-copula SAEM output (omega_dist = vine), surfaced as
# `fit$vine_copula` (marginals / pairs / params data frames) and the scalar
# `fit$vine_corrected_ofv`. See ferx-core "Vine-Copula SAEM".

# Cached short vine SAEM fit on the bundled warfarin_vine example. Iteration
# counts are overridden low: these tests check output shape, not convergence.
vine_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin_vine")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method = "saem", verbose = FALSE, covariance = FALSE,
        settings = list(
          n_exploration = 40L,
          n_convergence = 20L,
          omega_burnin  = 5L,
          seed          = 12345L
        )
      )
    }
    fit
  }
})

test_that("vine_copula is a list of marginals/pairs/params data frames", {
  fit <- vine_fit()

  expect_type(fit$vine_copula, "list")
  expect_named(fit$vine_copula, c("marginals", "pairs", "params"))

  expect_s3_class(fit$vine_copula$marginals, "data.frame")
  expect_identical(names(fit$vine_copula$marginals), c("eta", "mean", "sd"))
  expect_s3_class(fit$vine_copula$pairs, "data.frame")
  expect_identical(
    names(fit$vine_copula$pairs),
    c("tree", "pair", "family", "kendall_tau", "tail_dep_lower", "tail_dep_upper")
  )
  expect_s3_class(fit$vine_copula$params, "data.frame")
  expect_identical(
    names(fit$vine_copula$params),
    c("tree", "pair", "parameter", "estimate", "se")
  )
})

test_that("vine marginals carry one finite-sd Gaussian per random effect", {
  fit <- vine_fit()
  m <- fit$vine_copula$marginals

  # The warfarin_vine example declares two etas: ETA_CL, ETA_V.
  expect_setequal(m$eta, c("ETA_CL", "ETA_V"))
  expect_true(all(is.finite(m$sd)))
  expect_true(all(m$sd > 0))
})

test_that("vine pairs report a recognised family with a finite Kendall's tau", {
  fit <- vine_fit()
  p <- fit$vine_copula$pairs

  expect_gte(nrow(p), 1L)
  known <- c("gaussian", "student_t", "clayton", "gumbel", "frank")
  expect_true(all(p$family %in% known))
  expect_true(all(is.finite(p$kendall_tau)))
  expect_true(all(p$kendall_tau >= -1 & p$kendall_tau <= 1))
  # Tail dependence is a probability in [0, 1].
  expect_true(all(p$tail_dep_lower >= 0 & p$tail_dep_lower <= 1))
  expect_true(all(p$tail_dep_upper >= 0 & p$tail_dep_upper <= 1))
})

test_that("vine params list each copula parameter with a finite estimate", {
  fit <- vine_fit()
  pr <- fit$vine_copula$params

  expect_gte(nrow(pr), 1L)
  expect_true(all(nzchar(pr$parameter)))
  expect_true(all(is.finite(pr$estimate)))
})

test_that("vine_corrected_ofv is a finite scalar comparable to a Gaussian OFV", {
  fit <- vine_fit()
  expect_type(fit$vine_corrected_ofv, "double")
  expect_length(fit$vine_corrected_ofv, 1L)
  expect_true(is.finite(fit$vine_corrected_ofv))
})

test_that("a Gaussian SAEM fit carries no vine output", {
  ex <- ferx_example("warfarin_saem")
  fit <- ferx_fit(
    ex$model, ex$data,
    method = "saem", verbose = FALSE, covariance = FALSE,
    settings = list(n_exploration = 20L, n_convergence = 10L, omega_burnin = 5L)
  )
  expect_null(fit$vine_copula)
  expect_null(fit$vine_corrected_ofv)
})

test_that("vine output survives a .fitrx save/load round-trip", {
  fit <- vine_fit()

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_equal(loaded$vine_copula$marginals, fit$vine_copula$marginals)
  expect_equal(loaded$vine_copula$pairs, fit$vine_copula$pairs)
  expect_equal(loaded$vine_copula$params, fit$vine_copula$params)
  expect_equal(loaded$vine_corrected_ofv, fit$vine_corrected_ofv)
})

test_that("print.ferx_fit shows a VINE COPULA section for a vine fit", {
  fit <- vine_fit()
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "VINE COPULA")
  expect_match(out, "Pair-copulas")
})

test_that("ferx_vine_density returns a bivariate grid named after the etas", {
  fit <- vine_fit()
  d <- ferx_vine_density(fit, n = 10L)

  expect_s3_class(d, "data.frame")
  expect_identical(names(d), c("ETA_CL", "ETA_V", "density"))
  expect_equal(nrow(d), 100L) # n * n
  expect_true(all(is.finite(d$density)))
  expect_true(all(d$density >= 0))
})

test_that("ferx_vine_density evaluates a proper density (integrates to ~1)", {
  fit <- vine_fit()
  d <- ferx_vine_density(fit, n = 80L, range_sd = 5)

  ux <- sort(unique(d[[1]]))
  uy <- sort(unique(d[[2]]))
  cell <- (ux[2] - ux[1]) * (uy[2] - uy[1])
  integral <- sum(d$density) * cell
  expect_equal(integral, 1, tolerance = 0.05)
})

test_that("ferx_vine_density supports a 1-D marginal slice", {
  fit <- vine_fit()
  d <- ferx_vine_density(fit, etas = "ETA_CL", n = 15L)

  expect_identical(names(d), c("ETA_CL", "density"))
  expect_equal(nrow(d), 15L)
  expect_true(all(is.finite(d$density)))
})

test_that("ferx_vine_density errors clearly on a non-vine fit and unknown etas", {
  fit <- vine_fit()
  gauss <- fit
  gauss$vine_copula <- NULL
  expect_error(ferx_vine_density(gauss), "no vine-copula")
  expect_error(ferx_vine_density(fit, etas = "ETA_NOPE"), "Unknown eta")
})

test_that("ferx_simulate draws ETAs from the vine when the fit carries one", {
  fit <- vine_fit()
  ex <- ferx_example("warfarin_vine")

  s_vine <- suppressWarnings(
    ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 7L, fit = fit)
  )
  expect_s3_class(s_vine, "data.frame")
  expect_true(all(is.finite(s_vine$DV_SIM)))

  # Same fitted theta/omega/sigma, but with the vine stripped -> Gaussian eta
  # draws. With the vine present the simulated values must differ, proving the
  # copula path is actually used rather than the multivariate normal.
  gauss <- fit
  gauss$vine_copula <- NULL
  s_gauss <- suppressWarnings(
    ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 7L, fit = gauss)
  )
  expect_false(isTRUE(all.equal(s_vine$DV_SIM, s_gauss$DV_SIM)))
})
