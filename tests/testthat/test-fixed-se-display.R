# FIX parameters in print.ferx_fit() and fit$estimates (#451), and the
# missing-flags warning in ferx_save_fit() (#452).
#
# The engine reports SE 0 for a FIX entry. print() must say FIXED, and
# fit$estimates must carry NA SE / RSE with fixed = TRUE. The trigger is the
# fit's *_fixed flags, never SE == 0: a free parameter whose covariance step
# was skipped, failed or came back exactly 0 must keep printing its SE.

.compute_estimates <- getFromNamespace(".ferx_compute_estimates", "ferx")

# The print row for one label: a THETA table row starts with the name, the
# OMEGA / SIGMA / OMEGA_IOV rows are indented by two spaces. Whitespace then
# anything but `~` after the label skips the shrinkage line (`ETA_CL: 3%`) and
# an off-diagonal row that starts with the same eta (`ETA_V ~ ETA_CL`).
print_line <- function(out, label) {
  hit <- grep(sprintf("^\\s*%s\\s+[^~\\s]", label), out, value = TRUE, perl = TRUE)
  expect_length(hit, 1L)
  hit
}

# The SE field of a FIX row. KAPPA rows are ferx-core's own console rows
# (#470), which label a FIX kappa `NAME [FIX]` with `SE = ---`; every other
# block is R's and says FIXED.
fixed_se <- function(label) if (startsWith(label, "KAPPA")) "SE = ---" else "SE = FIXED"

# -- Crafted fits: every cell of class x fixed/free x covariance state --------

# One FIX and one free entry per class, all with SE 0 from the "engine": the
# pair only renders apart if the flags, not the SE, decide. Block omega
# (ETA_CL, ETA_V) is the FIX block, ETA_KA is free; the off-diagonal pair of
# the FIX block is FIX too.
crafted_fit <- function(se = c("zero", "real", "none"), flags = TRUE) {
  se <- match.arg(se)
  theta <- c(TVCL = 0.134, TVKA = 1.0)
  omega <- matrix(c(0.07, 0.02, 0, 0.02, 0.02, 0, 0, 0, 0.4), 3, 3,
                  dimnames = rep(list(c("ETA_CL", "ETA_V", "ETA_KA")), 2))
  sigma <- c(PROP_ERR = 0.1, ADD_ERR = 0.5)
  omega_iov <- matrix(c(0.04, 0, 0, 0.02), 2, 2,
                      dimnames = rep(list(c("KAPPA_CL", "KAPPA_V")), 2))
  # Lower triangle, column-major: (1,1) (2,1) (3,1) (2,2) (3,2) (3,3).
  se_omega_real <- c(0.011, 0.005, 0, 0.004, 0, 0.09)
  fit <- list(
    theta = theta,
    theta_names = names(theta),
    theta_transforms = c(TVCL = "identity", TVKA = "identity"),
    omega = omega, eta_names = colnames(omega),
    eta_param_types = rep("log_normal", 3L),
    sigma = sigma, sigma_names = names(sigma),
    sigma_types = c("proportional", "additive"),
    omega_iov = omega_iov, kappa_names = colnames(omega_iov),
    se_theta = switch(se, zero = c(TVCL = 0, TVKA = 0), real = c(TVCL = 0.02, TVKA = 0), none = NULL),
    se_omega = switch(se, zero = rep(0, 6L), real = replace(se_omega_real, c(1, 2, 4), 0), none = NULL),
    se_sigma = switch(se, zero = c(0, 0), real = c(0.01, 0), none = NULL),
    se_kappa = switch(se, zero = c(0, 0), real = c(0.008, 0), none = NULL)
  )
  if (flags) {
    fit$theta_fixed <- c(FALSE, TRUE)
    fit$omega_fixed <- c(TRUE, TRUE, FALSE)
    fit$sigma_fixed <- c(FALSE, TRUE)
    fit$kappa_fixed <- c(FALSE, TRUE)
  }
  class(fit) <- "ferx_fit"
  fit
}

fixed_labels <- c("TVKA", "ETA_CL", "ETA_V", "ADD_ERR", "KAPPA_V")
free_labels  <- c("TVCL", "ETA_KA", "PROP_ERR", "KAPPA_CL")

test_that("print: the FIX entry of every class says FIXED, its free twin keeps SE 0", {
  out <- capture.output(print(crafted_fit("zero")))
  # theta: SE and %RSE columns both read FIXED.
  expect_match(print_line(out, "TVKA"), "^TVKA\\s+1\\.000000\\s+FIXED\\s+FIXED\\s*$")
  for (lbl in setdiff(fixed_labels, "TVKA")) {
    ln <- print_line(out, lbl)
    expect_match(ln, fixed_se(lbl), fixed = TRUE, info = lbl)
    expect_no_match(ln, "SE = 0", fixed = TRUE)
  }
  # Free twins with an SE of exactly 0: printed as 0, not FIXED (kills an
  # SE == 0 shortcut).
  expect_match(print_line(out, "TVCL"), "^TVCL\\s+0\\.134000\\s+0\\.000000\\s+0\\.0\\s*$")
  for (lbl in setdiff(free_labels, "TVCL")) {
    ln <- print_line(out, lbl)
    expect_match(ln, "SE = 0.000000", fixed = TRUE, info = lbl)
    expect_no_match(ln, "FIXED", fixed = TRUE)
  }
  # The off-diagonal of the FIX block.
  expect_match(print_line(out, "ETA_V ~ ETA_CL"), "SE = FIXED\\s*$")
})

test_that("print: covariance step ok, failed or skipped - FIXED regardless, free per its SE", {
  real <- capture.output(print(crafted_fit("real")))
  expect_match(print_line(real, "TVCL"), "^TVCL\\s+0\\.134000\\s+0\\.020000\\s+14\\.9\\s*$")
  expect_match(print_line(real, "ETA_KA"), "SE = 0.090000", fixed = TRUE)
  expect_match(print_line(real, "PROP_ERR"), "SE = 0.010000", fixed = TRUE)
  expect_match(print_line(real, "KAPPA_CL"), "SE = 0.008000", fixed = TRUE)
  # Skipped (covariance = FALSE) or failed: no SE at all. FIX is still known,
  # the free entries say N/A.
  none <- capture.output(print(crafted_fit("none")))
  expect_match(print_line(none, "TVKA"), "FIXED\\s+FIXED\\s*$")
  expect_match(print_line(none, "TVCL"), "N/A\\s+N/A\\s*$")
  for (lbl in setdiff(fixed_labels, "TVKA")) {
    expect_match(print_line(none, lbl), fixed_se(lbl), fixed = TRUE, info = lbl)
  }
  for (lbl in setdiff(free_labels, "TVCL")) {
    expect_match(print_line(none, lbl), "SE = N/A", fixed = TRUE, info = lbl)
  }
  expect_match(print_line(none, "ETA_V ~ ETA_CL"), "SE = FIXED\\s*$")
})

test_that("print: a fit with no *_fixed fields does not crash and claims no FIXED", {
  out <- capture.output(print(crafted_fit("zero", flags = FALSE)))
  expect_false(any(grepl("FIXED", out, fixed = TRUE)))
  expect_match(print_line(out, "TVKA"), "^TVKA\\s+1\\.000000\\s+0\\.000000\\s+0\\.0\\s*$")
  expect_match(print_line(out, "ADD_ERR"), "SE = 0.000000", fixed = TRUE)
})

test_that("estimates: fixed column, NA SE / RSE / CI for FIX, free SE 0 kept", {
  est <- .compute_estimates(crafted_fit("zero"))
  expect_true("fixed" %in% names(est))
  expect_identical(est[fixed_labels, "fixed"], rep(TRUE, length(fixed_labels)))
  expect_identical(est[free_labels, "fixed"], rep(FALSE, length(free_labels)))
  for (col in c("se", "rse_pct", "lower_95", "upper_95")) {
    expect_true(all(is.na(est[fixed_labels, col])), info = col)
  }
  expect_identical(est[free_labels, "se"], rep(0, length(free_labels)))
  # The point estimate of a FIX row is untouched.
  expect_identical(est["TVKA", "estimate"], 1.0)
  # No flags: fixed = FALSE everywhere, SE as reported.
  old <- .compute_estimates(crafted_fit("zero", flags = FALSE))
  expect_identical(old$fixed, rep(FALSE, nrow(old)))
  expect_identical(old["TVKA", "se"], 0)
})

test_that("estimates: a FIX log or logit theta keeps its natural-scale estimate, without an interval", {
  for (tf in c("log", "logit")) {
    fit <- crafted_fit("zero")
    fit$theta_transforms[["TVKA"]] <- tf
    est <- .compute_estimates(fit)
    want <- if (tf == "log") exp(1.0) else stats::plogis(1.0)
    expect_equal(est["TVKA", "estimate_natural"], want, info = tf)
    expect_true(is.na(est["TVKA", "lower_95_natural"]), info = tf)
    expect_true(is.na(est["TVKA", "upper_95_natural"]), info = tf)
  }
})

test_that("ferx_se(): NA for a FIX parameter, without the no-standard-errors warning", {
  fit <- crafted_fit("zero")
  fit$estimates <- .compute_estimates(fit)
  expect_silent(se <- ferx_se(fit, c("TVKA", "ADD_ERR")))
  expect_identical(se, c(TVKA = NA_real_, ADD_ERR = NA_real_))
  # A free parameter with no SE still warns.
  none <- crafted_fit("none")
  none$estimates <- .compute_estimates(none)
  expect_warning(ferx_se(none, c("TVCL", "TVKA")), "no standard errors")
})

# -- ferx_save_fit() with flags missing (#452) ---------------------------------

test_that("ferx_save_fit() warns, per class, when a fit carries no FIX flags", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_iov")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE))
  f <- tempfile(fileext = ".fitrx")
  on.exit(unlink(f), add = TRUE)
  # With every flag present: quiet.
  expect_no_warning(ferx_save_fit(fit, f), message = "FIX flags")
  old <- fit
  old[c("theta_fixed", "omega_fixed", "sigma_fixed", "kappa_fixed")] <- NULL
  w <- character()
  withCallingHandlers(ferx_save_fit(old, f),
                      warning = function(cnd) {
                        w <<- c(w, conditionMessage(cnd))
                        invokeRestart("muffleWarning")
                      })
  w <- grep("FIX flags", w, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "theta, omega, sigma, kappa", fixed = TRUE)
  expect_match(w, "made before ferx recorded them", fixed = TRUE)
  expect_match(w, "recorded as estimated", fixed = TRUE)
  expect_match(w, "reload as free", fixed = TRUE)
  expect_match(w, "Refit", fixed = TRUE)
  # One class missing: only that one is named.
  part <- fit
  part$sigma_fixed <- NULL
  expect_warning(ferx_save_fit(part, f), "FIX flags for sigma are unknown")
})

# -- Live fits -----------------------------------------------------------------

mbma_live <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      ex <- ferx_example("mbma_placebo")
      cached <<- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE))
    }
    cached
  }
})

test_that("mbma_placebo: ADD_ERR (sigma FIX) prints FIXED; a free theta keeps its SE", {
  skip_on_cran()
  fit <- mbma_live()
  skip_if(is.null(fit$se_theta), "covariance step did not run - skipping")
  expect_identical(fit$sigma_fixed, TRUE)
  out <- capture.output(print(fit))
  ln <- print_line(out, "ADD_ERR")
  expect_match(ln, "SE = FIXED", fixed = TRUE)
  expect_no_match(ln, "SE = 0", fixed = TRUE)
  est <- fit$estimates
  expect_true(is.na(est["ADD_ERR", "se"]))
  expect_true(is.na(est["ADD_ERR", "rse_pct"]))
  expect_true(est["ADD_ERR", "fixed"])
  # TVE0 is free and estimated.
  expect_false(est["TVE0", "fixed"])
  expect_gt(est["TVE0", "se"], 0)
  expect_equal(est["TVE0", "se"], unname(fit$se_theta[["TVE0"]]))
  expect_match(print_line(out, "TVE0"),
               sprintf("^TVE0\\s+\\S+\\s+%s\\s", sprintf("%.6f", fit$se_theta[["TVE0"]])))
})

# The per-class live fixture: warfarin_iov with one FIX and one free entry per
# class, a FIX block omega among them.
per_class_model <- function() {
  ex <- ferx_example("warfarin_iov")
  m <- c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  theta TVV(8.1, 0.1, 500.0)",
    "  theta TVKA(1.0, 0.01, 50.0) FIX",
    "  block_omega (ETA_CL, ETA_V) = [0.07, 0.01, 0.02] FIX",
    "  omega ETA_KA ~ 0.40",
    "  kappa KAPPA_CL ~ 0.04",
    "  kappa KAPPA_V ~ 0.02 FIX",
    "  sigma PROP_ERR ~ 0.1 (sd)",
    "  sigma ADD_ERR ~ 0.1 (sd) FIX",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL + KAPPA_CL)",
    "  V  = TVV  * exp(ETA_V + KAPPA_V)",
    "  KA = TVKA * exp(ETA_KA)",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ combined(PROP_ERR, ADD_ERR)",
    "[fit_options]",
    "  method     = foce",
    "  iov_column = OCC",
    "  covariance = true"
  )
  p <- tempfile(fileext = ".ferx")
  writeLines(m, p)
  list(model = p, data = ex$data)
}

per_class_live <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      x <- per_class_model()
      cached <<- suppressWarnings(ferx_fit(x$model, x$data, verbose = FALSE))
    }
    cached
  }
})

test_that("live fit, one FIX and one free per class: flags, print and estimates agree", {
  skip_on_cran()
  fit <- per_class_live()
  expect_identical(fit$theta_fixed, c(FALSE, FALSE, TRUE))
  expect_identical(fit$omega_fixed, c(TRUE, TRUE, FALSE))
  expect_identical(fit$sigma_fixed, c(FALSE, TRUE))
  expect_identical(fit$kappa_fixed, c(FALSE, TRUE))
  out <- capture.output(print(fit))
  expect_match(print_line(out, "TVKA"), "FIXED\\s+FIXED\\s*$")
  for (lbl in c("ETA_CL", "ETA_V", "ADD_ERR", "KAPPA_V")) {
    expect_match(print_line(out, lbl), fixed_se(lbl), fixed = TRUE, info = lbl)
  }
  expect_match(print_line(out, "ETA_V ~ ETA_CL"), "SE = FIXED\\s*$")
  est <- fit$estimates
  expect_identical(est[c("TVKA", "ETA_CL", "ETA_V", "ADD_ERR", "KAPPA_V"), "fixed"], rep(TRUE, 5L))
  expect_identical(est[c("TVCL", "TVV", "ETA_KA", "PROP_ERR", "KAPPA_CL"), "fixed"], rep(FALSE, 5L))
  expect_true(all(is.na(est[est$fixed, "se"])))
  skip_if(is.null(fit$se_theta), "covariance step did not run - skipping")
  for (lbl in c("TVCL", "ETA_KA", "PROP_ERR", "KAPPA_CL")) {
    ln <- print_line(out, lbl)
    expect_no_match(ln, "FIXED", fixed = TRUE)
    expect_false(is.na(est[lbl, "se"]), info = lbl)
  }
})

test_that("live fit with covariance = FALSE: FIX still FIXED, free entries N/A", {
  skip_on_cran()
  x <- per_class_model()
  fit <- suppressWarnings(ferx_fit(x$model, x$data, covariance = FALSE, verbose = FALSE))
  expect_null(fit$se_theta)
  out <- capture.output(print(fit))
  expect_match(print_line(out, "TVKA"), "FIXED\\s+FIXED\\s*$")
  expect_match(print_line(out, "TVCL"), "N/A\\s+N/A\\s*$")
  expect_match(print_line(out, "ADD_ERR"), "SE = FIXED", fixed = TRUE)
  expect_match(print_line(out, "PROP_ERR"), "SE = N/A", fixed = TRUE)
})

test_that(".fitrx round trip prints the same SE / FIXED as the live fit", {
  skip_on_cran()
  cases <- list(
    per_class = list(fit = per_class_live(),
                     theta = c("TVCL", "TVV", "TVKA"),
                     other = c("ETA_CL", "ETA_V", "ETA_KA", "ETA_V ~ ETA_CL",
                               "KAPPA_CL", "KAPPA_V", "PROP_ERR", "ADD_ERR"),
                     fixed = c("TVKA", "ETA_CL", "ETA_V", "ETA_V ~ ETA_CL",
                               "KAPPA_V", "ADD_ERR")),
    mbma = list(fit = mbma_live(), theta = c("TVE0", "EMAX", "ED50", "ET50"),
                other = c("ETA_E0", "KAPPA_ARM", "ADD_ERR"), fixed = "ADD_ERR")
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    f <- tempfile(fileext = ".fitrx")
    expect_no_warning(ferx_save_fit(cs$fit, f), message = "FIX flags")
    fit2 <- suppressWarnings(ferx_load_fit(f))
    unlink(f)
    a <- capture.output(print(cs$fit))
    b <- capture.output(print(fit2))
    # A THETA row is the whole story: estimate, SE and %RSE.
    for (lbl in cs$theta) {
      expect_identical(print_line(b, lbl), print_line(a, lbl), info = paste(nm, lbl))
    }
    # Omega / sigma / kappa rows: compare the SE field. The rest of a sigma
    # row (its type tag) does not survive a reload today - unrelated to FIX.
    se_field <- function(out, lbl) regmatches(print_line(out, lbl),
                                              regexpr("SE = \\S+", print_line(out, lbl)))
    for (lbl in cs$other) {
      expect_identical(se_field(b, lbl), se_field(a, lbl), info = paste(nm, lbl))
    }
    for (lbl in setdiff(cs$fixed, cs$theta)) {
      expect_identical(se_field(b, lbl), fixed_se(lbl), info = paste(nm, lbl))
    }
    expect_identical(fit2$estimates$fixed, cs$fit$estimates$fixed)
    expect_identical(fit2$estimates$se, cs$fit$estimates$se)
  }
})
