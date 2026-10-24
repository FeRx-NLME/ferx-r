# plot_param_dist() (#478): the implied population density of each
# ETA-bearing parameter, by change of variables from N(0, omega) through the
# parameter's scale and shape, and the panels drawn from it.

param_density <- getFromNamespace(".ferx_param_density", "ferx")
shape_h <- getFromNamespace(".ferx_shape_h", "ferx")
shape_dh <- getFromNamespace(".ferx_shape_dh", "ferx")

# Trapezoid integral of f over x.
trapz <- function(x, f) sum(diff(x) * (head(f, -1L) + tail(f, -1L)) / 2)

# -- the density ----------------------------------------------------------------

test_that("the implied density integrates to 1 on every scale and shape", {
  cases <- list(
    list(type = "log_normal", base = log(5), x = seq(1e-4, 60, length.out = 40001L)),
    list(type = "additive", base = 5, x = seq(-5, 15, length.out = 40001L)),
    list(type = "logit", base = qlogis(0.3), x = seq(1e-6, 1 - 1e-6, length.out = 40001L)),
    list(type = "log_normal", base = log(5), shape = list(kind = "boxcox", value = 0.5),
         x = seq(1e-4, 60, length.out = 40001L)),
    list(type = "log_normal", base = log(5), shape = list(kind = "boxcox", value = -0.5),
         x = seq(1e-4, 60, length.out = 40001L)),
    list(type = "log_normal", base = log(5), shape = list(kind = "tdist", value = 5),
         x = seq(1e-4, 200, length.out = 80001L)),
    list(type = "log_normal", base = log(5), shape = list(kind = "johndraper", value = 0.5),
         x = seq(1e-4, 60, length.out = 40001L))
  )
  for (cs in cases) {
    f <- param_density(cs$x, cs$base, 0.25, cs$type, cs$shape)
    expect_equal(trapz(cs$x, f), 1, tolerance = 2e-3,
                 label = paste(cs$type, cs$shape$kind %||% "", cs$shape$value %||% ""))
  }
})

test_that("the unshaped scales are the textbook densities", {
  x <- seq(0.5, 20, length.out = 200L)
  expect_equal(param_density(x, log(5), 0.1, "log_normal"),
               dlnorm(x, log(5), sqrt(0.1)), tolerance = 1e-5)
  xa <- seq(3, 7, length.out = 200L)
  expect_equal(param_density(xa, 5, 0.4, "additive"), dnorm(xa, 5, sqrt(0.4)),
               tolerance = 1e-5)
  xl <- seq(0.05, 0.95, length.out = 200L)
  logitnorm <- dnorm(qlogis(xl), qlogis(0.3), sqrt(0.5)) / (xl * (1 - xl))
  expect_equal(param_density(xl, qlogis(0.3), 0.5, "logit"), logitnorm, tolerance = 1e-5)
})

test_that("box-cox at lambda -> 0 is the log-normal density", {
  x <- seq(0.5, 20, length.out = 200L)
  ln <- dlnorm(x, log(5), sqrt(0.1))
  for (lam in c(1e-3, 1e-6, 0)) {
    f <- param_density(x, log(5), 0.1, "log_normal", list(kind = "boxcox", value = lam))
    expect_equal(f, ln, tolerance = 5e-3, label = paste("lambda", lam))
  }
})

test_that("a subject mixture averages the per-subject densities", {
  x <- seq(0.5, 30, length.out = 300L)
  mix <- param_density(x, log(c(4, 8)), 0.1, "log_normal")
  expect_equal(mix, (dlnorm(x, log(4), sqrt(0.1)) + dlnorm(x, log(8), sqrt(0.1))) / 2,
               tolerance = 1e-5)
})

test_that("a custom or unknown scale has no implied density", {
  expect_true(all(is.na(param_density(1:3, 1, 0.1, "custom"))))
  expect_true(all(is.na(param_density(1:3, 1, 0.1, "unknown"))))
})

# -- the shape kernels ----------------------------------------------------------

test_that("each shape reduces to the identity at its identity value", {
  eta <- seq(-2, 2, by = 0.25)
  expect_equal(shape_h(eta, "boxcox", 0), eta)
  expect_equal(shape_h(eta, "johndraper", 1), eta)
  expect_equal(shape_h(eta, "tdist", 1e8), eta, tolerance = 1e-7)
  expect_identical(shape_h(eta), eta)
})

test_that("the shape kernels match ferx-core's closed forms", {
  # boxcox: (exp(eta)^lambda - 1) / lambda; johndraper:
  # sign(eta) * ((|eta| + 1)^lambda - 1) / lambda (ferx-core eta_shape.rs).
  eta <- c(-1.3, -0.2, 0.4, 1.7)
  expect_equal(shape_h(eta, "boxcox", 0.6), (exp(eta)^0.6 - 1) / 0.6)
  expect_equal(shape_h(eta, "johndraper", 0.6), sign(eta) * ((abs(eta) + 1)^0.6 - 1) / 0.6)
  nu <- 5
  td <- eta * (1 + (eta^2 + 1) / (4 * nu) + (5 * eta^4 + 16 * eta^2 + 3) / (96 * nu^2) +
                 (3 * eta^6 + 19 * eta^4 + 17 * eta^2 - 15) / (384 * nu^3))
  expect_equal(shape_h(eta, "tdist", nu), td)
})

test_that("the kernel derivatives match finite differences", {
  eta <- c(-1.3, -0.2, 0.4, 1.7)
  d <- 1e-6
  for (sh in list(list("boxcox", 0.6), list("boxcox", -1), list("tdist", 5),
                  list("johndraper", 0.6))) {
    fd <- (shape_h(eta + d, sh[[1]], sh[[2]]) - shape_h(eta - d, sh[[1]], sh[[2]])) / (2 * d)
    expect_equal(shape_dh(eta, sh[[1]], sh[[2]]), fd, tolerance = 1e-6, label = sh[[1]])
  }
})

# -- the plot -------------------------------------------------------------------

# A warfarin model whose V is additive (ETA_V) and CL log-normal (ETA_CL).
additive_v_fit <- function(env = parent.frame()) {
  model <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  theta TVV(8.1, 0.1, 500.0)",
    "  theta TVKA(1.0, 0.01, 50.0)",
    "  omega ETA_CL ~ 0.07",
    "  omega ETA_V  ~ 1.0",
    "  sigma PROP_ERR ~ 0.01 (sd)",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
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
    "  covariance = false"
  ), model)
  ferx_fit(model, ferx_example("warfarin")$data, verbose = FALSE,
           settings = list(maxiter = 5L))
}

test_that("plot_param_dist() draws a warfarin fit without error", {
  skip_on_cran()
  fit <- warfarin_fit()
  expect_identical(fit$eta_param_info$individual_param, c("CL", "V", "KA"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_param_dist(fit)
  expect_named(out, c("CL", "V", "KA"))
  expect_identical(out$CL$eta, "ETA_CL")
  expect_identical(out$CL$param_type, "log_normal")
  # Plainly log-normal: no reference curve, and the density covers its range.
  expect_true(all(is.na(out$CL$density$reference)))
  expect_gt(trapz(out$CL$density$x, out$CL$density$density), 0.9)
  expect_silent(plot_param_dist(fit, param = "V", qq = TRUE))
  expect_error(plot_param_dist(fit, param = "Q"), "No ETA-bearing individual parameter")
})

test_that("an additive ETA gets the log-normal reference", {
  skip_on_cran()
  fit <- additive_v_fit()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_param_dist(fit, qq = TRUE)
  expect_setequal(names(out), c("CL", "V"))
  expect_identical(out$V$param_type, "additive")
  expect_true(any(is.finite(out$V$density$reference)))
})

test_that("eta_param_info survives a .fitrx round trip", {
  skip_on_cran()
  fit <- warfarin_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$eta_param_info, fit$eta_param_info)
  # The engine's own reader and writer keep the parameter names too.
  trip <- fitrx_engine_trip(fit)
  expect_identical(ferx_load_fit(trip$core_path)$eta_param_info$individual_param,
                   c("CL", "V", "KA"))
  # A loaded fit has no individual_estimates, and says so.
  expect_error(plot_param_dist(loaded), "carries no `individual_estimates`")
})

test_that("a bundle saved before eta_param_info loads with no parameter names", {
  skip_on_cran()
  fit <- warfarin_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  old <- fitrx_edit_bundle(path, function(w) {
    w$eta_param_info <- lapply(w$eta_param_info, function(e) {
      e$individual_param_name <- NULL
      e
    })
    w
  })
  loaded <- ferx_load_fit(old)
  expect_identical(loaded$eta_param_info$individual_param, c("", "", ""))
  expect_identical(loaded$eta_param_types, fit$eta_param_types)
})

test_that("without parameter names the plot pairs ETA_X with X", {
  skip_on_cran()
  fit <- warfarin_fit()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  unnamed <- fit
  unnamed$eta_param_info$individual_param <- ""
  expect_named(plot_param_dist(unnamed), c("CL", "V", "KA"))
  bare <- fit
  bare$eta_param_info <- NULL
  out <- plot_param_dist(bare)
  expect_named(out, c("CL", "V", "KA"))
  expect_identical(out$V$param_type, "log_normal")
})

test_that("a shape the engine reports is drawn and tagged", {
  skip_on_cran()
  fit <- warfarin_fit()
  # What ferx-core #359 P1 part B will supply; built by hand until it lands.
  fit$eta_param_info$shape_kind <- c("boxcox", NA, NA)
  fit$eta_param_info$shape_value <- c(0.42, NA, NA)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_param_dist(fit, param = "CL")
  expect_identical(out$CL$shape, list(kind = "boxcox", value = 0.42))
  expect_true(any(is.finite(out$CL$density$reference)))
  tag <- getFromNamespace(".ferx_param_dist_tag", "ferx")
  expect_identical(tag(list(param_type = "log_normal", shape = out$CL$shape)),
                   "log-normal, box-cox lambda=0.42")
})

test_that("conddist = TRUE plots conditional-distribution draws", {
  skip_on_cran()
  ex <- ferx_example("warfarin_saem")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, method = "saem",
                                   covariance = FALSE, verbose = FALSE,
                                   settings = list(conddist = TRUE)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_param_dist(fit, conddist = TRUE, n_draws = 5L, qq = TRUE)
  expect_length(out$CL$values, 5L * nrow(fit$individual_estimates))
  expect_true(all(is.finite(out$CL$values) & out$CL$values > 0))
  expect_error(plot_param_dist(warfarin_fit(), conddist = TRUE))
})
