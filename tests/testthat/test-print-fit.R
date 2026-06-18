# make_fake_fit() lives in helper-trace.R and is auto-loaded by testthat.

# print.ferx_fit CV% formula --------------------------------------------

test_that("print.ferx_fit uses exact log-normal CV% when eta_param_types is absent (fallback)", {
  # No eta_param_types field -> defaults to log_normal.
  # omega = 0.40: exact CV% = sqrt(exp(0.40) - 1) * 100 ≈ 70.1 (not 63.2 from sqrt(0.40)*100)
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(length(omega_line) == 1L)
  expect_false(grepl("63\\.2", omega_line))  # old approximate value
  expect_true(grepl("70\\.1", omega_line))   # exact log-normal CV%
})

test_that("print.ferx_fit uses exact log-normal CV% when eta_param_types is 'log_normal'", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "log_normal")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("70\\.1", omega_line))
})

test_that("print.ferx_fit shows SD_logit for logit eta_param_types", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "logit")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("logit", omega_line))
  expect_true(grepl("SD_logit", omega_line))
  # CV% is not meaningful for logit ETAs — must not appear
  expect_false(grepl("CV% =", omega_line))
})

test_that("print.ferx_fit CV% equals zero when omega diagonal is zero", {
  fit <- make_fake_fit(omega = matrix(0, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("CV% = 0\\.0", omega_line))
})

# SIGMA display (ferx#59 / ferx-core#57): sigma is on the SD scale, so
# print() must show variance = sigma^2 always and CV% = sigma * 100 only
# for proportional components.

test_that("print.ferx_fit shows variance + CV% for proportional sigma using its declared name", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.1,
    sigma_names = "PROP_ERR",
    sigma_types = "proportional"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("PROP_ERR", out)]
  expect_length(sigma_line, 1L)
  # variance = 0.1^2 = 0.01; CV% = 0.1 * 100 = 10.0.
  expect_true(grepl("var = 0\\.010000", sigma_line))
  expect_true(grepl("CV% = 10\\.0", sigma_line))
  # No legacy SIGMA(1) label when a name is supplied.
  expect_false(any(grepl("SIGMA\\(1\\)", out)))
})

test_that("print.ferx_fit shows variance but no CV% for additive sigma", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.5,
    sigma_names = "ADD_ERR",
    sigma_types = "additive"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("ADD_ERR", out)]
  expect_length(sigma_line, 1L)
  # variance = 0.5^2 = 0.25; no CV% on observation-unit scale.
  expect_true(grepl("var = 0\\.250000", sigma_line))
  expect_false(grepl("CV%", sigma_line))
})

test_that("print.ferx_fit handles combined error: CV% on prop component only", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = c(0.2, 0.5),
    sigma_names = c("PROP_ERR", "ADD_ERR"),
    sigma_types = c("proportional", "additive")
  )
  out <- capture.output(print(fit))
  prop_line <- out[grepl("PROP_ERR", out)]
  add_line  <- out[grepl("ADD_ERR",  out)]
  expect_length(prop_line, 1L)
  expect_length(add_line,  1L)
  expect_true(grepl("CV% = 20\\.0", prop_line))
  expect_false(grepl("CV%", add_line))
  expect_true(grepl("var = 0\\.040000", prop_line))
  expect_true(grepl("var = 0\\.250000", add_line))
})

test_that("print.ferx_fit falls back to SIGMA(i) when sigma_names is missing", {
  # sigma_names absent (older Rust binary or unit test that didn't pass them).
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.1,
    sigma_names = NULL,
    sigma_types = NULL
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("SIGMA\\(1\\)", out)))
})

# Block-omega correlations (#60): print uses omega_param_corr when the engine
# provides it (bivariate-lognormal formula), and falls back to the eta-level
# Pearson formula when it is absent.

test_that("print.ferx_fit uses omega_param_corr value and 'param corr' label when present", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  pc <- matrix(c(1.0, 0.5227, 0.5227, 1.0), 2, 2)
  fit <- make_fake_fit(
    omega            = om,
    omega_param_corr = pc,
    eta_param_types  = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("param corr = 0\\.5227", corr_line))
})

test_that("print.ferx_fit falls back to eta-level corr label when omega_param_corr is NULL", {
  # cov = 0.025, vars = 0.10 -> Pearson corr = 0.025 / 0.10 = 0.25
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega            = om,
    omega_param_corr = NULL,
    eta_param_types  = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("\\(corr = 0\\.2500\\)", corr_line))
  expect_false(grepl("param corr", corr_line))
})

test_that("print.ferx_fit uses omega_iov_param_corr when present for IOV correlations", {
  om     <- matrix(0.10, 1, 1)
  iov    <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  iov_pc <- matrix(c(1.0, 0.5227, 0.5227, 1.0), 2, 2)
  fit <- make_fake_fit(
    omega                 = om,
    omega_iov             = iov,
    omega_iov_param_corr  = iov_pc,
    kappa_names           = c("KAPPA1", "KAPPA2"),
    se_kappa              = NULL,
    shrinkage_kappa       = NULL
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("KAPPA2 ~ KAPPA1", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("param corr = 0\\.5227", corr_line))
})

test_that("print.ferx_fit uses exact log-normal CV% for OMEGA_IOV when kappa_param_types absent", {
  # kappa = 0.20 -> exact CV% = sqrt(exp(0.20) - 1) * 100 ≈ 47.1 (not 44.7 from sqrt(0.20)*100)
  fit <- make_fake_fit(
    omega           = matrix(0.10, 1, 1),
    omega_iov       = matrix(0.20, 1, 1),
    kappa_names     = "KAPPA1",
    se_kappa        = NULL,
    shrinkage_kappa = NULL
  )
  out <- capture.output(print(fit))
  kappa_line <- out[grepl("KAPPA1", out)]
  expect_true(length(kappa_line) >= 1L)
  expect_false(grepl("44\\.7", kappa_line[1]))  # old approximate value
  expect_true(grepl("47\\.1", kappa_line[1]))   # exact log-normal CV%
})

# NN-aware output (Option E, ferx-core#52): per-weight thetas are hidden from
# the THETA table and a compact NEURAL NETWORKS block summarises them.

test_that("print.ferx_fit skips NN-weight thetas and emits NEURAL NETWORKS block", {
  # 2 baseline thetas (TVCL, TVV1) + 4 NN weights at offset=2.
  theta <- c(TVCL = 4.0, TVV1 = 40.0, W1 = 0.1, W2 = -0.2, W3 = 0.3, W4 = 0.0)
  nn <- list(list(
    name              = "TYPICAL_PK",
    shape             = c(2L, 2L, 2L),
    hidden_activation = "tanh",
    output_activation = "exp",
    n_weights         = 4L,
    weights_offset    = 2L,
    input_names       = c("WT", "CRCL"),
    output_names      = c("CL", "V1"),
    weights           = c(0.1, -0.2, 0.3, 0.0)
  ))
  fit <- make_fake_fit(
    omega           = matrix(0.10, 1, 1),
    theta           = theta,
    neural_networks = nn
  )
  out <- capture.output(print(fit))

  # Baseline thetas appear in the THETA table.
  expect_true(any(grepl("^TVCL\\s", out)))
  expect_true(any(grepl("^TVV1\\s", out)))
  # NN-weight thetas are skipped from the THETA table (no `W1`/`W2`/... rows).
  theta_block <- out[seq(
    which(grepl("^\\s*THETA\\s*$|--- THETA", out))[1L],
    (which(grepl("NEURAL NETWORKS", out)) - 1L)[1L]
  )]
  expect_false(any(grepl("^W[1-4]\\s", theta_block)))

  # NEURAL NETWORKS block contains the network metadata and weight summary.
  expect_true(any(grepl("NEURAL NETWORKS", out)))
  expect_true(any(grepl("TYPICAL_PK.*shape=\\[2, 2, 2\\].*activation=tanh/exp.*n_weights=4", out)))
  expect_true(any(grepl("inputs:\\s+\\[WT, CRCL\\]", out)))
  expect_true(any(grepl("outputs:\\s+\\[CL, V1\\]", out)))
  expect_true(any(grepl("weights: min -0\\.2000\\s+max 0\\.3000", out)))
})

test_that("print.ferx_fit omits NEURAL NETWORKS block when neural_networks is empty/NULL", {
  fit <- make_fake_fit(omega = matrix(0.10, 1, 1))
  out <- capture.output(print(fit))
  expect_false(any(grepl("NEURAL NETWORKS", out)))
})

# --- init_as_sd annotation tests ---

test_that("print.ferx_fit annotates omega with [initial specified as SD] when flag is TRUE", {
  fit <- make_fake_fit(
    omega           = matrix(0.09, 1, 1),
    eta_names       = "ETA_CL",
    omega_init_as_sd = TRUE
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("ETA_CL", out)]
  expect_true(length(omega_line) > 0L)
  expect_true(any(grepl("\\[initial specified as SD\\]", omega_line)))
})

test_that("print.ferx_fit does NOT annotate omega when flag is FALSE", {
  fit <- make_fake_fit(
    omega           = matrix(0.09, 1, 1),
    eta_names       = "ETA_CL",
    omega_init_as_sd = FALSE
  )
  out <- capture.output(print(fit))
  expect_false(any(grepl("\\[initial specified as SD\\]", out)))
})

test_that("print.ferx_fit annotates sigma with [initial specified as SD] when flag is TRUE", {
  fit <- make_fake_fit(
    omega            = matrix(0.09, 1, 1),
    sigma            = c(PROP_ERR = 0.01),
    sigma_names      = "PROP_ERR",
    sigma_types      = "proportional",
    sigma_init_as_sd = TRUE
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("PROP_ERR", out)]
  expect_true(length(sigma_line) > 0L)
  expect_true(any(grepl("\\[initial specified as SD\\]", sigma_line)))
})

test_that("print.ferx_fit handles missing omega_init_as_sd gracefully (old fit object)", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1), eta_names = "ETA_CL")
  # omega_init_as_sd not set at all — should not error, no annotation
  out <- capture.output(print(fit))
  expect_false(any(grepl("\\[initial specified as SD\\]", out)))
})

# STATUS line surfaces failing critical category --------------------------

test_that("print.ferx_fit appends critical category to STATUS when fit did not converge", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = FALSE,
    warnings_structured = data.frame(
      severity = c("critical", "info"),
      category = c("convergence", "mu_referencing"),
      message  = c("Outer optimization did not converge",
                   "mu-ref: CL, V"),
      source_method = c("", ""),
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_length(status_line, 1L)
  expect_true(grepl("NOT CONVERGED", status_line))
  # Inline category appears between the label and the iteration/time tail.
  expect_true(grepl("\\(convergence\\)", status_line))
})

test_that("print.ferx_fit deduplicates and joins multiple critical categories", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = FALSE,
    warnings_structured = data.frame(
      severity = c("critical", "critical", "critical"),
      category = c("convergence", "covariance_step", "convergence"),
      message  = c("a", "b", "c"),
      source_method = c("", "", ""),
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_true(grepl("\\(convergence, covariance_step\\)", status_line))
})

test_that("print.ferx_fit STATUS has no inline category when there are no critical warnings", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = TRUE,
    warnings_structured = data.frame(
      severity = "info",
      category = "mu_referencing",
      message  = "mu-ref: CL",
      source_method = "",
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_true(grepl("CONVERGED", status_line))
  # No parenthesised category snippet because nothing is critical.
  expect_false(grepl("\\(mu_referencing\\)|\\(convergence\\)", status_line))
})

test_that("print.ferx_fit STATUS has no inline category when non-converged fit has only info warnings", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = FALSE,
    warnings_structured = data.frame(
      severity = "info",
      category = "mu_referencing",
      message  = "mu-ref: CL",
      source_method = "",
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_true(grepl("NOT CONVERGED", status_line))
  expect_false(grepl("\\(", status_line))
})

# se_omega display (#143): print() shows "SE = ..." for omega diagonal and
# block off-diagonal elements, reading via .omega_se_at(). Block omega carries
# a full lower-triangle se_omega; diagonal omega carries one SE per variance.

test_that("print.ferx_fit shows SE for omega diagonal when se_omega is present", {
  fit <- make_fake_fit(
    omega    = matrix(0.40, 1, 1),
    se_omega = 0.05
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_length(omega_line, 1L)
  expect_true(grepl("SE = 0\\.050000", omega_line))
})

test_that("print.ferx_fit shows off-diagonal SE for block omega (full lower-triangle se_omega)", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  # Full lower-triangle, column-major: [(1,1), (2,1), (2,2)]
  fit <- make_fake_fit(
    omega    = om,
    se_omega = c(0.011, 0.022, 0.033),
    eta_param_types = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  # Diagonal SE on the OMEGA(1,1) line.
  diag_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("SE = 0\\.011000", diag_line))
  # Off-diagonal covariance line carries the (2,1) SE.
  cov_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(cov_line, 1L)
  expect_true(grepl("SE = 0\\.022000", cov_line))
})

test_that("print.ferx_fit omits off-diagonal SE when se_omega is diagonal-only", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega    = om,
    se_omega = c(0.011, 0.033),  # length n_eta -> diagonal-only, no cov SE
    eta_param_types = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  cov_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(cov_line, 1L)
  expect_false(grepl("SE =", cov_line))
})

test_that("print.ferx_fit shows off-diagonal SE with named etas (block omega)", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega     = om,
    se_omega  = c(0.011, 0.022, 0.033),
    eta_names = c("ETA_CL", "ETA_V"),
    eta_param_types = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  # Loop emits the higher-index eta first (i outer, j < i), so the (2,1)
  # pair prints as "ETA_V ~ ETA_CL".
  cov_line <- out[grepl("ETA_V ~ ETA_CL : cov =", out)]
  expect_length(cov_line, 1L)
  expect_true(grepl("SE = 0\\.022000", cov_line))
})
