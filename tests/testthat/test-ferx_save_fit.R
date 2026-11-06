
# ---- header from test-bayes.R ----
# Tier-1 tests for the Bayes (method = "bayes") R-side surface: the R-hat
# convergence flag, the print/summary rendering of `$bayes`, and the persist
# wire builder. These exercise the wrapper logic with a mocked posterior so
# they run in milliseconds — the real MCMC engine round-trip is covered by the
# (slow) integration tests in test-fit.R / test-persist.R.

# A minimal but complete `$bayes` posterior summary, shaped exactly as
# `fit_result_to_list()` (Rust) and `ferx_load_fit()` (persist) produce it.
fake_bayes <- function(max_rhat = 1.004) {
  list(
    n_chains          = 2,
    n_warmup          = 300,
    n_draws_per_chain = 600,
    n_divergent       = 0,
    max_rhat          = max_rhat,
    param_names       = c("TVCL", "TVV", "OMEGA_CL", "SIGMA_PROP"),
    mean              = c(0.133, 8.0, 0.09, 0.10),
    sd                = c(0.01, 0.4, 0.02, 0.01),
    q025              = c(0.115, 7.2, 0.05, 0.08),
    median            = c(0.133, 8.0, 0.09, 0.10),
    q975              = c(0.151, 8.8, 0.13, 0.12),
    rhat              = c(1.001, 1.002, 1.004, 1.000),
    ess_bulk          = c(820, 760, 540, 910),
    ess_tail          = c(800, 740, 520, 900),
    mcse              = c(3e-4, 0.01, 8e-4, 3e-4)
  )
}

# .ferx_bayes_rhat_flag ------------------------------------------------------



# print.ferx_fit — BAYES posterior block -------------------------------------




# summary.ferx_fit / print.ferx_summary --------------------------------------



# .fitrx_build_bayes_wire ----------------------------------------------------




# ---- header from test-persist.R ----
# ferx_save_fit() / ferx_load_fit() round-trip and CLI-flag tests.
#
# These exercise the .fitrx bundle format. The cross-tool test (Rust CLI
# writes, R reads) lives in ferx-core; here we only verify the R-side
# writer + reader are inverses, and that the `output` argument of
# ferx_fit() saves a valid bundle.

# Local aliases avoid the ::: operator (undesirable_operator_linter).
.fitrx_build_iov_wire   <- getFromNamespace(".fitrx_build_iov_wire",   "ferx")





















test_that("ferx_save_fit + ferx_load_fit round-trip a real fit", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_s3_class(loaded, "ferx_fit")
  expect_equal(loaded$method, fit$method)
  expect_equal(loaded$converged, fit$converged)
  expect_equal(loaded$ofv, fit$ofv, tolerance = 1e-12)
  expect_equal(loaded$aic, fit$aic, tolerance = 1e-12)
  expect_equal(loaded$bic, fit$bic, tolerance = 1e-12)
  expect_equal(unname(loaded$theta), unname(fit$theta), tolerance = 1e-12)
  expect_equal(names(loaded$theta), names(fit$theta))
  expect_equal(loaded$omega, fit$omega, tolerance = 1e-12)
  expect_equal(unname(loaded$sigma), unname(fit$sigma), tolerance = 1e-12)
  expect_equal(loaded$n_subjects, fit$n_subjects)
  expect_equal(loaded$n_obs, fit$n_obs)
  expect_equal(loaded$shrinkage_eta, fit$shrinkage_eta, tolerance = 1e-12)
  expect_equal(loaded$model_name, fit$model_name)
  expect_equal(loaded$ferx_version, fit$ferx_version)
})

test_that("bundle doubles are written and read back bit for bit", {
  # 15 significant digits (jsonlite's `digits = NA`, write.table()) cannot
  # identify a double, and R's own parser is not correctly rounded on every
  # platform, so both halves are checked across the whole double range.
  set.seed(20261002)
  x <- c(rnorm(5000) * 10^runif(5000, -300, 300), 0.1, 1 / 3, 5e-324,
         .Machine$double.xmax, -280.36396183775315, 1, 1e22)
  expect_identical(ferx:::.fitrx_parse_doubles(ferx:::.fitrx_double_text(x)), x)
  # The shortest text that round-trips, not always 17 digits.
  expect_identical(ferx:::.fitrx_double_text(c(0.1, 1, 1e22)), c("0.1", "1", "1e+22"))
  expect_identical(
    ferx:::.fitrx_parse_doubles(c("", "NA", "NaN", "Inf", "-Inf", "inf", "1.5")),
    c(NA, NA, NaN, Inf, -Inf, Inf, 1.5)
  )
  # A CSV entry: text stays text, missing stays missing, doubles are exact.
  df <- data.frame(ID = c("a,1", "2", "3"), x = c(0.1 + 0.2, NA, 1 / 3),
                   n = 1:3, stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_csv_exact(df, path)
  expect_identical(ferx:::.fitrx_read_csv_exact(path), df)

  # The random doubles above only catch an inexact parser on aarch64: x86_64
  # R parses through an 80-bit long double and reads every %.17g string
  # exactly, so on the ubuntu CI runner they pass under as.numeric() too. It
  # is still not correctly rounded on arbitrary text, and this literal is one
  # it misreads (as.numeric() gives ...4c5p+714 there, ...4cap+714 on
  # aarch64), so these two checks fail on both if the parse half is reverted.
  hard <- "1.6268116685806102e215"
  expect_identical(ferx:::.fitrx_parse_doubles(hard), 0x1.e33cb7bcbb4c6p+714)
  writeLines(c("ID,x", paste0("1,", hard)), path)
  expect_identical(ferx:::.fitrx_read_csv_exact(path)$x, 0x1.e33cb7bcbb4c6p+714)
})

test_that("a real fit round-trips bit for bit", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  # fit.json, names and dimnames included (#417)
  expect_identical(loaded$theta, fit$theta)
  expect_identical(loaded$omega, fit$omega)
  expect_identical(loaded$sigma, fit$sigma)
  expect_identical(loaded$cov_matrix, fit$cov_matrix)
  expect_identical(loaded$ofv, fit$ofv)
  expect_identical(loaded$se_theta, fit$se_theta)
  # the CSV entries, column types included (#417)
  expect_identical(loaded$sdtab, fit$sdtab)
  expect_identical(loaded$ebe_etas, fit$ebe_etas)
  # and what a reloaded fit drives
  ex <- ferx_example("warfarin")
  expect_identical(
    ferx_predict(ex$model, ex$data, fit = loaded),
    ferx_predict(ex$model, ex$data, fit = fit)
  )
  expect_identical(
    ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 7L, fit = loaded),
    ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 7L, fit = fit)
  )
})
test_that("a round trip keeps names, dimnames and column types (#417)", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "no covariance step in this fit")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  # The live fit names sigma like theta, and the bundle gives them back.
  expect_identical(names(fit$sigma), fit$sigma_names)
  expect_identical(loaded$sigma, fit$sigma)
  expect_identical(loaded$sigma_types, fit$sigma_types)
  # fit.json has no row names for the covariance matrix; the loader asks the
  # engine for the labels the live fit was named with.
  expect_false(is.null(dimnames(fit$cov_matrix)))
  expect_identical(dimnames(loaded$cov_matrix), dimnames(fit$cov_matrix))
  expect_identical(loaded$cov_matrix, fit$cov_matrix)
  expect_identical(loaded$cov_fixed, fit$cov_fixed)
  expect_identical(loaded$cor_matrix, fit$cor_matrix)
  # read.csv() types a column of whole numbers as integer; every sdtab column
  # and every eta is double on the live fit.
  expect_identical(loaded$sdtab, fit$sdtab)
  expect_identical(loaded$ebe_etas, fit$ebe_etas)
  expect_identical(vapply(loaded$sdtab, typeof, ""),
                   vapply(fit$sdtab, typeof, ""))
})

test_that("bundle tables keep double columns double when every value is whole (#417)", {
  df <- data.frame(ID = c("1", "2"), TIME = c(0, 24), WT = c(70, 100000),
                   EVID = c(1L, 0L), EMPTY = c(NA_real_, NA_real_),
                   stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_csv_exact(df, path)
  # Without the live types, read.csv() decides: WT integer, EMPTY logical.
  plain <- ferx:::.fitrx_read_csv_exact(path, id_cols = "ID")
  expect_type(plain$WT, "integer")
  expect_type(plain$EMPTY, "logical")
  expect_identical(
    ferx:::.fitrx_read_csv_exact(path, id_cols = "ID", integer_cols = "EVID"),
    df
  )
})

test_that("a cov_matrix is named on load for the engine-written bundle too (#417)", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "no covariance step in this fit")
  trip <- fitrx_engine_trip(fit)
  on.exit(unlink(c(trip$r_path, trip$core_path)), add = TRUE)
  # ferx-core's own writer stores no names either, so this is the path every
  # bundle written before #417 takes.
  loaded <- suppressWarnings(ferx_load_fit(trip$core_path))
  expect_identical(dimnames(loaded$cov_matrix), dimnames(fit$cov_matrix))
  expect_identical(dimnames(loaded$cor_matrix), dimnames(fit$cor_matrix))
})

test_that("a level-block fit's cov_matrix is named on load (#417)", {
  skip_on_cran()
  # The level thetas exist only once the model is laid out on the fit's
  # levels: labels from the unbound model are too few for the matrix.
  ex <- ferx_example("mbma_placebo")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE))
  skip_if(is.null(fit$cov_matrix), "no covariance step in this fit")
  expect_gt(nrow(fit$theta_levels), 0L)
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_false(is.null(dimnames(fit$cov_matrix)))
  expect_identical(loaded$cov_matrix, fit$cov_matrix)
  expect_identical(loaded$cov_fixed, fit$cov_fixed)
  expect_identical(loaded$cor_matrix, fit$cor_matrix)
})

test_that("bare (extension-less) output path defaults to .fitrx and round-trips (#268)", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  dir <- tempfile("fitrx_ext_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  bare <- file.path(dir, "run1_base_diag")
  ferx_save_fit(fit, bare)

  # The archive is written with the conventional `.fitrx` extension, not a
  # bare name and not Info-ZIP's `.zip`.
  expect_true(file.exists(paste0(bare, ".fitrx")))
  expect_false(file.exists(bare))
  expect_false(file.exists(paste0(bare, ".zip")))

  # ... and loads back from the same bare path the caller saved to.
  loaded <- ferx_load_fit(bare)
  expect_s3_class(loaded, "ferx_fit")
  expect_equal(unname(loaded$theta), unname(fit$theta), tolerance = 1e-12)

  # An explicit extension is honoured verbatim (not forced to .fitrx).
  custom <- file.path(dir, "run1.bundle")
  ferx_save_fit(fit, custom)
  expect_true(file.exists(custom))
  expect_s3_class(ferx_load_fit(custom), "ferx_fit")
})
test_that("covariance step results survive round-trip", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "no covariance step in this fit")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_true(is.matrix(loaded$cov_matrix))
  expect_equal(dim(loaded$cov_matrix), dim(fit$cov_matrix))
  expect_equal(loaded$cov_matrix, fit$cov_matrix, tolerance = 1e-12)
  expect_equal(unname(loaded$se_theta), unname(fit$se_theta), tolerance = 1e-12)
})
test_that("the model-selection surface survives round-trip", {
  skip_on_cran()
  fit <- warfarin_fit_cov()

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  # `bic_inputs` and `left_init` travel at wire level (they are ferx-core
  # FitWire fields); the three derived verdicts travel under r_extras.
  expect_equal(loaded$bic_inputs, fit$bic_inputs)
  expect_equal(loaded$left_init, fit$left_init)
  expect_identical(loaded$stationary_start, fit$stationary_start)
  expect_equal(loaded$stalled_at_init, fit$stalled_at_init)
  expect_equal(loaded$estimate_near_boundary, fit$estimate_near_boundary)
  expect_equal(loaded$max_abs_correlation, fit$max_abs_correlation,
               tolerance = 1e-12)

  for (ty in c("mixed", "fixed", "iiv", "random")) {
    expect_equal(ferx_bic(loaded, ty), ferx_bic(fit, ty), tolerance = 1e-12)
  }
  # The condition number and the packed-Omega layout travel at wire level too.
  # Reading only the wire spelling on save wrote a null for every bundle this
  # package writes, so the condition-number gate landed in `skipped` (not
  # `failures`) on every reload - compare both halves of the verdict.
  expect_equal(
    loaded$cov_condition_number %||% loaded$condition_number,
    fit$condition_number %||% fit$cov_condition_number,
    tolerance = 1e-8
  )
  expect_equal(loaded$omega_is_diagonal, fit$omega_is_diagonal)
  expect_equal(loaded$kappa_is_diagonal, fit$kappa_is_diagonal)

  # The engine's stalled_at_init() checks init/estimate shapes before reading
  # `left_init`, so the inits have to survive for the verdict to mean anything.
  expect_equal(loaded$theta_init, fit$theta_init, tolerance = 1e-12)
  expect_equal(loaded$sigma_init, fit$sigma_init, tolerance = 1e-12)
  expect_equal(loaded$omega_init_dim, as.integer(fit$omega_init_dim))
  expect_equal(loaded$omega_init, as.numeric(fit$omega_init), tolerance = 1e-12)

  strict_loaded <- check_strictness(loaded, max_condition_number = 1e10)
  strict_fresh  <- check_strictness(fit, max_condition_number = 1e10)
  expect_equal(strict_loaded$failures, strict_fresh$failures)
  expect_equal(strict_loaded$skipped, strict_fresh$skipped)
  expect_equal(strict_loaded$passed, strict_fresh$passed)
})

test_that("a warm start at an optimum reloads as one, not as a stall", {
  skip_on_cran()
  # Refit from the fitted estimates: the engine's restart probe confirms the
  # start (ferx-core #2506). The bundle has to carry `stationary_start`: the
  # engine's stalled_at_init() reads it beside `left_init`, and without it
  # calls the reloaded fit a stall. R's own `stalled_at_init` travels in
  # r_extras and is not recomputed on load, so it is no oracle here.
  ex <- ferx_example("warfarin")
  model <- tempfile(fileext = ".ferx")
  on.exit(unlink(model), add = TRUE)
  ferx_update_model(warfarin_fit(), model)
  fit <- suppressWarnings(ferx_fit(model, ex$data, method = "focei",
                                   verbose = FALSE,
                                   settings = list(maxiter = 30L)))
  expect_true(fit$stationary_start)
  expect_false(fit$left_init)
  expect_false(fit$stalled_at_init)

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  expect_true(ferx_load_fit(path)$stationary_start)

  # The engine's reader takes the key R wrote: its writer emits
  # `stationary_start` only when the FitResult it read holds TRUE.
  trip <- fitrx_engine_trip(fit)
  on.exit(unlink(c(trip$r_path, trip$core_path)), add = TRUE)
  expect_true(isTRUE(fitrx_fit_json(trip$core_path)$stationary_start))
  expect_true(ferx_load_fit(trip$core_path)$stationary_start)
})

test_that("bayes posterior summary survives round-trip", {
  skip_on_cran()
  fit <- warfarin_bayes_fit()
  skip_if(is.null(fit$bayes), "no bayes result in this fit")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$bayes))
  b0 <- fit$bayes
  b1 <- loaded$bayes
  expect_equal(b1$n_chains, b0$n_chains)
  expect_equal(b1$n_warmup, b0$n_warmup)
  expect_equal(b1$n_draws_per_chain, b0$n_draws_per_chain)
  expect_equal(b1$max_rhat, b0$max_rhat, tolerance = 1e-12)
  expect_equal(b1$param_names, b0$param_names)
  expect_equal(b1$mean, b0$mean, tolerance = 1e-12)
  expect_equal(b1$q975, b0$q975, tolerance = 1e-12)
  expect_equal(b1$ess_bulk, b0$ess_bulk, tolerance = 1e-9)
})
test_that("SAEM conditional-distribution results survive round-trip (real fit)", {
  skip_on_cran()
  fit <- warfarin_saem_conddist_fit()
  skip_if(is.null(fit$cond_dist), "no cond_dist result in this fit")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$cond_dist))
  expect_equal(loaded$cond_dist$data$ID, fit$cond_dist$data$ID)
  expect_equal(loaded$cond_dist$data$ETA, fit$cond_dist$data$ETA)
  expect_equal(loaded$cond_dist$data$COND_MEAN, fit$cond_dist$data$COND_MEAN,
               tolerance = 1e-9)
  expect_equal(loaded$cond_dist$data$COND_SD, fit$cond_dist$data$COND_SD,
               tolerance = 1e-9)
  expect_equal(loaded$cond_dist$data$COND_MODE, fit$cond_dist$data$COND_MODE,
               tolerance = 1e-9)
})
test_that(".fitrx_has_cond_dist / .fitrx_write_conddist_csv (synthetic)", {
  no_cd <- structure(list(cond_dist = NULL), class = "ferx_fit")
  expect_false(ferx:::.fitrx_has_cond_dist(no_cd))

  cd_data <- data.frame(
    ID = c("1", "1", "2", "2"),
    ETA = c("ETA_CL", "ETA_V", "ETA_CL", "ETA_V"),
    COND_MEAN = c(0.05, -0.02, -0.03, 0.01),
    COND_SD = c(0.12, 0.15, 0.11, 0.14),
    COND_MODE = c(0.04, -0.01, -0.02, 0.02),
    stringsAsFactors = FALSE
  )
  with_cd <- structure(
    list(cond_dist = list(data = cd_data, shrinkage = c(0.1, 0.05),
                           nsamp = 20L, burnin = 5L)),
    class = "ferx_fit"
  )
  expect_true(ferx:::.fitrx_has_cond_dist(with_cd))

  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_conddist_csv(with_cd, path)
  written <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(written$ID, c(1L, 1L, 2L, 2L))
  expect_equal(written$ETA, cd_data$ETA)
  expect_equal(written$COND_MEAN, cd_data$COND_MEAN, tolerance = 1e-12)
})
test_that("predictions survive round-trip (per-row equality)", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_true(is.data.frame(loaded$sdtab))
  expect_equal(nrow(loaded$sdtab), nrow(fit$sdtab))
  numeric_cols <- intersect(
    c("TIME", "DV", "PRED", "IPRED", "CWRES", "IWRES", "EBE_OFV"),
    intersect(names(loaded$sdtab), names(fit$sdtab))
  )
  for (col in numeric_cols) {
    expect_equal(
      as.numeric(loaded$sdtab[[col]]),
      as.numeric(fit$sdtab[[col]]),
      tolerance = 1e-9,
      info = sprintf("column %s", col)
    )
  }
})
test_that("ferx_fit(output = ...) writes a valid bundle", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  fit <- ferx_fit(
    ex$model, ex$data,
    method = "focei", verbose = FALSE, covariance = FALSE,
    settings = list(maxiter = 10L),
    output = path
  )
  expect_true(file.exists(path))

  loaded <- ferx_load_fit(path)
  expect_s3_class(loaded, "ferx_fit")
  expect_equal(loaded$ofv, fit$ofv, tolerance = 1e-12)
})
test_that("ferx_fit(output =) relays save-side conditions under its own name (#537)", {
  relay <- ferx:::.ferx_relay_save
  # A real save-side warning: a fit that carries no FIX flags. (The normal
  # SIR proposal this used to take is recorded since ferx-core#1819.)
  w <- tryCatch(
    relay(ferx:::.fitrx_warn_missing_fixed(list(theta = c(TVCL = 1)))),
    warning = function(w) w
  )
  expect_match(conditionMessage(w), "^ferx_fit: the FIX flags for theta are unknown")
  expect_null(conditionCall(w))
  # The real save-side refusal, prefixed by the engine's decoder.
  sc <- ferx:::.ferx_default_settings("scoring_settings")
  sc$inner_maxit <- 5L
  expect_error(
    relay(ferx:::.fitrx_settings_to_wire(list(scoring_settings = sc), "scoring_settings")),
    'ferx_fit: `fit$scoring_settings`: unknown key "inner_maxit"', fixed = TRUE
  )
  # The `ferx_save_fit(): ` spelling, with its class kept for handlers.
  e <- tryCatch(relay(stop(errorCondition(
    "ferx_save_fit(): the data file x has changed", class = "ferx_data_edited_error"
  ))), error = function(e) e)
  expect_s3_class(e, "ferx_data_edited_error")
  expect_identical(conditionMessage(e), "ferx_fit: the data file x has changed")
  w <- tryCatch(relay(warning(warningCondition(
    "ferx_save_fit(): the model file x has changed", class = "ferx_model_edited_warning"
  ))), warning = function(w) w)
  expect_s3_class(w, "ferx_model_edited_warning")
  expect_identical(conditionMessage(w), "ferx_fit: the model file x has changed")
  # Only the leading prefix is rewritten; anything else passes through as is,
  # and the value comes back.
  expect_warning(relay(warning("see ferx_save_fit: below")),
                 "^see ferx_save_fit: below$")
  expect_error(relay(stop("unrelated")), "^unrelated$")
  expect_identical(relay(42L), 42L)
})
test_that("ferx_fit(output =) never names ferx_save_fit in its warnings (#537)", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  testthat::local_mocked_bindings(
    ferx_save_fit = function(fit, path, ...) {
      warning("ferx_save_fit: not recorded", call. = FALSE)
      invisible(path)
    }
  )
  expect_warning(
    ferx_fit(
      ex$model, ex$data,
      method = "focei", verbose = FALSE, covariance = FALSE,
      settings = list(maxiter = 5L),
      output = path
    ),
    "^ferx_fit: not recorded$"
  )
})
test_that("include_data = TRUE bundles the input CSV", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  fit <- ferx_fit(
    ex$model, ex$data,
    method = "focei", verbose = FALSE, covariance = FALSE,
    settings = list(maxiter = 5L),
    output = path, include_data = TRUE
  )
  expect_true(file.exists(path))

  entries <- utils::unzip(path, list = TRUE)$Name
  expect_true("data.csv" %in% entries)

  loaded <- ferx_load_fit(path)
  expect_true(!is.na(loaded$data_path))
  expect_true(file.exists(loaded$data_path))
})
test_that("manifest reports format_version 2", {
  skip_on_cran()
  fit <- warfarin_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)

  staging <- tempfile("fitrx_inspect_")
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  dir.create(staging)
  utils::unzip(path, exdir = staging)
  manifest <- jsonlite::read_json(file.path(staging, "manifest.json"),
                                  simplifyVector = TRUE)
  expect_equal(manifest$format_version, "2")
  expect_true("fit.json" %in% manifest$entries)
})
test_that("ferx_save_fit rejects non-ferx_fit inputs", {
  expect_error(ferx_save_fit(list(), tempfile()), "ferx_fit")
})
test_that("save -> load -> save again works (warfarin)", {
  # End-to-end check: load a fit from a bundle and save it again. Warfarin's
  # numeric-string IDs come back from read.csv as integers, so this path
  # exercises only the numeric branch of the predictions writer.
  skip_on_cran()
  fit <- warfarin_fit()
  first <- tempfile(fileext = ".fitrx")
  second <- tempfile(fileext = ".fitrx")
  on.exit(unlink(c(first, second)), add = TRUE)

  ferx_save_fit(fit, first)
  reloaded <- ferx_load_fit(first)
  ferx_save_fit(reloaded, second)
  twice <- ferx_load_fit(second)
  expect_equal(twice$ofv, fit$ofv, tolerance = 1e-12)
})
test_that("init_as_sd flags survive ferx_save_fit / ferx_load_fit round-trip", {
  skip_if_not(nzchar(Sys.getenv("FERX_RUN_REAL_FIT", "")),
              "Set FERX_RUN_REAL_FIT=1 to run real-fit persist tests")
  fit <- warfarin_fit()
  # Inject flags so we can round-trip them without a real (sd)-annotated fit
  fit <- editable_fit(fit)
  fit$omega_init_as_sd <- c(TRUE, FALSE)
  fit$sigma_init_as_sd <- TRUE
  fit$kappa_init_as_sd <- logical(0L)

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(resealed_fit(fit), path)
  loaded <- ferx_load_fit(path)

  expect_equal(loaded$omega_init_as_sd, c(TRUE, FALSE))
  expect_equal(loaded$sigma_init_as_sd, TRUE)
  expect_equal(loaded$kappa_init_as_sd, logical(0L))
})
test_that("LTBS residual label survives a .fitrx round-trip", {
  skip_on_cran()
  # A fit of a log-transform-both-sides model must report its residual as
  # "additive (log-transformed)" (set by the Rust glue's residual_label from
  # CompiledModel.log_transform), and that label must survive save/load via the
  # r_extras-persisted model_structure.
  md <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.2, 0.001, 10.0)",
    "  theta TVV(10.0, 0.1, 500.0)",
    "  theta TVKA(1.5, 0.01, 50.0)",
    "  omega ETA_CL ~ 0.09",
    "  sigma ADD_LOG ~ 0.15 (sd)",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "  KA = TVKA",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  log(DV) ~ additive(ADD_LOG)"
  ), md)
  dd <- tempfile(fileext = ".csv")
  writeLines(c(
    "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV",
    "1,0,.,1,100,1,0,1",
    "1,1,8.5,0,.,1,0,0",
    "1,4,5.1,0,.,1,0,0",
    "1,12,1.9,0,.,1,0,0",
    "2,0,.,1,120,1,0,1",
    "2,1,9.8,0,.,1,0,0",
    "2,4,6.0,0,.,1,0,0",
    "2,12,2.2,0,.,1,0,0"
  ), dd)
  on.exit(unlink(c(md, dd)), add = TRUE)

  fit <- ferx_fit(md, data = dd, verbose = FALSE,
                  settings = list(maxiter = 30L))
  expect_equal(fit$model_structure$residual, "additive (log-transformed)")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_equal(loaded$model_structure$residual, "additive (log-transformed)")
})
test_that("FIX flags of an in-memory fit survive a .fitrx round-trip (#436)", {
  skip_on_cran()
  # Regression: the glue shipped no theta/omega/sigma/kappa_fixed, so
  # ferx_save_fit() fell back to rep(FALSE, ...) and the bundle claimed every
  # FIX parameter was estimated. One FIX and one free entry per class, so a
  # flag forced to either value is caught.
  ex <- ferx_example("warfarin_iov")
  md <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  theta TVV(8.1, 0.1, 500.0)",
    "  theta TVKA(1.0, FIX)",
    "  omega ETA_CL ~ 0.07",
    "  omega ETA_V  ~ 0.02 FIX",
    "  kappa KAPPA_CL ~ 0.04 FIX",
    "  kappa KAPPA_V  ~ 0.01",
    "  sigma PROP_ERR ~ 0.01 (sd)",
    "  sigma ADD_ERR  ~ 0.1 (sd) FIX",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL + KAPPA_CL)",
    "  V  = TVV  * exp(ETA_V + KAPPA_V)",
    "  KA = TVKA",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ combined(PROP_ERR, ADD_ERR)",
    "[fit_options]",
    "  method     = foce",
    "  iov_column = OCC",
    "  covariance = false"
  ), md)
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(c(md, path)), add = TRUE)

  fit <- ferx_fit(md, data = ex$data, verbose = FALSE,
                  settings = list(maxiter = 5L))
  expect_identical(fit$theta_fixed, c(FALSE, FALSE, TRUE))
  expect_identical(fit$omega_fixed, c(FALSE, TRUE))
  expect_identical(fit$sigma_fixed, c(FALSE, TRUE))
  expect_identical(fit$kappa_fixed, c(TRUE, FALSE))

  # ferx_runlog() reads the same flags: the held value, marked (FIX), in the
  # INITIAL and FINAL columns of the fixed rows only (#450). No covariance
  # step, so there is no SE column to say FIXED.
  runlog_rows <- function(f) {
    log_lines <- strsplit(ferx_runlog(f, verbose = FALSE), "\n", fixed = TRUE)[[1]]
    function(nm) grep(paste0("^\\s*", nm, "\\s"), log_lines, value = TRUE)[1]
  }
  check_runlog <- function(f, what) {
    row_of <- runlog_rows(f)
    held <- c(TVKA = f$theta[["TVKA"]], ETA_V = f$omega[2, 2], ADD_ERR = f$sigma[[2]])
    esc <- function(s) gsub("([.()+])", "\\\\\\1", s)
    for (nm in names(held)) {
      # INITIAL and FINAL both read "<value> (FIX)": the value was held.
      cell <- esc(paste(sprintf("%.6g", held[[nm]]), "(FIX)"))
      expect_match(row_of(nm), sprintf("^\\s*%s\\s+%s\\s+%s\\s*$", nm, cell, cell),
                   label = paste(what, nm))
      expect_no_match(row_of(nm), "FIXED", fixed = TRUE, label = paste(what, nm))
    }
    # expect_no_match(NA, ...) passes, so a missing row must fail on its own.
    for (nm in c("TVCL", "ETA_CL", "PROP_ERR")) {
      expect_false(is.na(row_of(nm)), label = paste(what, nm, "row missing"))
      expect_no_match(row_of(nm), "FIX", fixed = TRUE, label = paste(what, nm))
    }
  }
  expect_identical(fit$theta[["TVKA"]], 1)
  check_runlog(fit, "live")

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  for (f in c("theta_fixed", "omega_fixed", "sigma_fixed", "kappa_fixed")) {
    expect_identical(loaded[[f]], fit[[f]], label = paste0("loaded$", f))
  }
  check_runlog(loaded, "loaded")
})
test_that(".fitrx_build_iov_wire uses empty list when shrinkage_kappa_by_occ is absent", {
  fake_iov_fit <- list(
    omega_iov              = matrix(0.09, 1, 1),
    kappa_names            = "KAPPA_CL",
    kappa_fixed            = FALSE,
    se_kappa               = NULL,
    shrinkage_kappa        = 0.15,
    shrinkage_kappa_by_occ = NULL,
    omega_iov_param_corr   = NULL
  )
  wire <- .fitrx_build_iov_wire(fake_iov_fit)
  expect_equal(wire$shrinkage_kappa_by_occ, list())
})
test_that(".fitrx_build_iov_wire writes the required weight keys empty when unweighted", {
  # Regression: an explicit null fails core with "invalid type: null, expected
  # a sequence", and v2 requires the keys, so they are arrays, empty here.
  fake_iov_fit <- list(
    omega_iov       = matrix(0.09, 1, 1),
    kappa_names     = "KAPPA_CL",
    kappa_fixed     = FALSE,
    shrinkage_kappa = 0.15,
    kappa_weights   = NULL
  )
  json <- as.character(jsonlite::toJSON(.fitrx_build_iov_wire(fake_iov_fit),
                                        auto_unbox = TRUE, null = "null", na = "null"))
  expect_match(json, '"kappa_weights":[]', fixed = TRUE)
  expect_match(json, '"kappa_weight_typical":[]', fixed = TRUE)
})
test_that(".fitrx_build_iov_wire writes the #1031 weight keys as null-holed arrays", {
  fake_iov_fit <- list(
    omega_iov            = matrix(0.09, 2, 2),
    kappa_names          = c("KAPPA_CL", "KAPPA_V"),
    kappa_fixed          = c(FALSE, FALSE),
    shrinkage_kappa      = c(0.15, 0.2),
    kappa_weights        = c(NA_character_, "NARM"),
    kappa_weight_typical = c(NA_real_, 400)
  )
  wire <- .fitrx_build_iov_wire(fake_iov_fit)
  expect_equal(wire$kappa_weights, list(NA_character_, "NARM"))
  expect_equal(wire$kappa_weight_typical, list(NA_real_, 400))
  # `na = "null"` is what the writer uses; the holes must survive as JSON null
  # so `Vec<Option<..>>` keeps every entry lined up with its kappa.
  expect_equal(as.character(jsonlite::toJSON(wire$kappa_weights,
                                             auto_unbox = TRUE, na = "null")),
               '[null,"NARM"]')
})
test_that("shrinkage_kappa_by_occ survives a ferx_save/ferx_load round-trip", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_iov")
  fit <- ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE)
  skip_if(is.null(fit$shrinkage_kappa_by_occ),
          "IOV fit did not produce shrinkage_kappa_by_occ - skip round-trip")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_equal(colnames(loaded$shrinkage_kappa_by_occ),
               colnames(fit$shrinkage_kappa_by_occ))
  expect_equal(nrow(loaded$shrinkage_kappa_by_occ),
               nrow(fit$shrinkage_kappa_by_occ))
  expect_equal(loaded$shrinkage_kappa_by_occ$occ,
               fit$shrinkage_kappa_by_occ$occ)
})
# ---- header from test-persist-helpers.R ----
# Unit tests for the pure-R serialization helpers in persist.R. These convert
# between the in-memory ferx_fit fields and the cross-language .fitrx wire
# schema, and need no model fit to exercise — so they live in the fast tier.

# ---------------------------------------------------------------------------
# Method name <-> token round-trip
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Covariance status <-> token
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Error model inference from sigma types
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Matrix <-> wire (row-major)
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# Confidence-interval <-> wire
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Optional scalar wrappers (NULL/NA/empty collapse to NULL)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Named SE unwrap + subject string IDs
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Per-subject OFV / N_OBS extraction from sdtab
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# EBEs CSV writer — header-only path when there are no random effects
# ---------------------------------------------------------------------------


test_that(".fitrx_method_to_token normalises labels and defaults to focei", {
  expect_identical(ferx:::.fitrx_method_to_token(NULL), "focei")
  expect_identical(ferx:::.fitrx_method_to_token(""), "focei")
  expect_identical(ferx:::.fitrx_method_to_token("FOCEI"), "focei")   # case-insensitive
  expect_identical(ferx:::.fitrx_method_to_token("FOCE"), "foce")
  expect_identical(ferx:::.fitrx_method_to_token("foce-gn"), "foce_gn")
  expect_identical(ferx:::.fitrx_method_to_token("foce_gn"), "foce_gn")
  expect_identical(ferx:::.fitrx_method_to_token("foce-gn-hybrid"), "foce_gn_hybrid")
  expect_identical(ferx:::.fitrx_method_to_token("SAEM"), "saem")
  expect_identical(ferx:::.fitrx_method_to_token("VI"), "vi")
  expect_identical(ferx:::.fitrx_method_to_token("LAPLACE"), "laplace")
  expect_identical(ferx:::.fitrx_method_to_token("Mystery"), "mystery") # unknown -> tolower
})
test_that(".fitrx_covariance_status_to_token normalises and defaults", {
  expect_identical(ferx:::.fitrx_covariance_status_to_token(NULL), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token(""), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("Computed"), "computed")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("FAILED"), "failed")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("NotRequested"), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("sir_fallback"), "sir_fallback")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("SirFallback"), "sir_fallback")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("weird"), "weird")
})
test_that(".fitrx_error_model_from_sigma_types classifies sigma vectors", {
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(NULL), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(character(0)), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("proportional", "Proportional")), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("additive")), "additive")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("additive", "proportional")), "combined")
})
test_that(".fitrx_matrix_to_wire emits row-major data and dims", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE) # [[1,2,3],[4,5,6]]
  w <- ferx:::.fitrx_matrix_to_wire(m)
  expect_identical(w$rows, 2L)
  expect_identical(w$cols, 3L)
  expect_identical(w$data, c(1, 2, 3, 4, 5, 6)) # row-major
})
test_that(".fitrx_matrix_to_wire returns NULL for NULL / non-matrix input", {
  expect_null(ferx:::.fitrx_matrix_to_wire(NULL))
  expect_null(ferx:::.fitrx_matrix_to_wire(c(1, 2, 3)))
})
test_that(".fitrx_ci_to_wire accepts a 2-column matrix", {
  ci <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE) # rows (1,2) (3,4)
  w  <- ferx:::.fitrx_ci_to_wire(ci)
  expect_length(w, 2L)
  expect_identical(w[[1L]], c(1, 2))
  expect_identical(w[[2L]], c(3, 4))
})
test_that(".fitrx_ci_to_wire accepts the legacy flat layout and rejects odd length", {
  expect_equal(ferx:::.fitrx_ci_to_wire(c(1, 2, 3, 4)), list(c(1, 2), c(3, 4)))
  expect_null(ferx:::.fitrx_ci_to_wire(c(1, 2, 3)))   # odd length
  expect_null(ferx:::.fitrx_ci_to_wire(NULL))
})
test_that("CI round-trips matrix -> wire -> matrix", {
  ci <- matrix(c(0.5, 1.5, 2.5, 3.5), ncol = 2, byrow = TRUE,
               dimnames = list(c("TVCL", "TVV"), c("lower", "upper")))
  wire <- ferx:::.fitrx_ci_to_wire(ci)
  back <- ferx:::.ferx_sir_ci_matrices(wire, NULL, NULL, rownames(ci), NULL, NULL)
  expect_identical(back$sir_ci_theta, ci)
  expect_null(back$sir_ci_omega)
  expect_null(back$sir_ci_sigma)
})
test_that("SIR CI rows fall back to THETA<i> / OMEGA(i,i) / SIGMA(i) (#482)", {
  ci <- ferx:::.ferx_sir_ci_matrices(
    c(1, 2, 3, 4), list(c(1, 2), c(3, 4)), c(5, 6),
    theta_names = "only_one", eta_names = NULL, sigma_names = character()
  )
  expect_identical(rownames(ci$sir_ci_theta), c("THETA1", "THETA2"))
  expect_identical(rownames(ci$sir_ci_omega), c("OMEGA(1,1)", "OMEGA(2,2)"))
  expect_identical(rownames(ci$sir_ci_sigma), "SIGMA(1)")
  expect_identical(colnames(ci$sir_ci_sigma), c("lower", "upper"))
})
test_that("SIR CI rows fall back per element on an empty or NA name (#557)", {
  ci <- ferx:::.ferx_sir_ci_matrices(
    c(1, 2, 3, 4), c(1, 2, 3, 4), c(1, 2, 3, 4),
    theta_names = c("TVCL", NA), eta_names = c("ETA_CL", ""),
    sigma_names = c(NA, "EPS_ADD")
  )
  expect_identical(rownames(ci$sir_ci_theta), c("TVCL", "THETA2"))
  expect_identical(rownames(ci$sir_ci_omega), c("ETA_CL", "OMEGA(2,2)"))
  expect_identical(rownames(ci$sir_ci_sigma), c("SIGMA(1)", "EPS_ADD"))
  # The kappa intervals and the kappa print labels apply the same rule.
  expect_identical(
    rownames(ferx:::.ferx_sir_ci_kappa(c(1, 2, 3, 4), c("", "KAPPA_V"))),
    c("KAPPA1", "KAPPA_V")
  )
  expect_identical(ferx:::.ferx_kappa_labels(c("", "KAPPA_V"), 2L),
                   c("KAPPA1", "KAPPA_V"))
})
test_that("a malformed SIR CI vector is refused, not dropped or shifted (#557)", {
  # Odd length on the live (FFI) path: no row structure.
  expect_error(
    ferx:::.ferx_sir_ci_matrices(c(1, 2, 3), NULL, NULL, NULL, NULL, NULL),
    "`sir_ci_theta` is not a \\[lower, upper\\] number pair per row \\(got 3 values\\)"
  )
  expect_error(ferx:::.ferx_sir_ci_kappa(c(1, 2, 3), "KAPPA_CL"), "`sir_ci_kappa`")
  # A pair list from a bundle: a `null` cell, a one-cell pair, a text cell.
  bad <- list(
    list(list(1, NULL), list(3, 4)),
    list(list(1), list(3, 4)),
    list(list(1, "x"), list(3, 4)),
    list(1, 2, 3, 4)
  )
  for (b in bad) {
    expect_error(
      ferx:::.ferx_sir_ci_matrices(NULL, b, NULL, NULL, c("A", "B"), NULL, bundle = TRUE),
      "fit.json: `sir.ci_omega` is not a list of \\[lower, upper\\] number pairs; the bundle is corrupt",
      info = deparse(b)
    )
  }
  expect_error(ferx:::.ferx_sir_ci_kappa(bad[[1L]], c("A", "B"), bundle = TRUE),
               "`sir.ci_kappa`.*the bundle is corrupt")
  # A well-formed pair list still loads, in row order.
  ok <- ferx:::.ferx_sir_ci_matrices(NULL, list(list(1, 2), list(3, 4)), NULL,
                                      NULL, c("A", "B"), NULL, bundle = TRUE)
  expect_identical(ok$sir_ci_omega,
                   matrix(c(1, 3, 2, 4), 2L, dimnames = list(c("A", "B"), c("lower", "upper"))))
})
test_that(".fitrx_subject_string_ids extracts character IDs or NULL", {
  expect_null(ferx:::.fitrx_subject_string_ids(list()))
  expect_null(ferx:::.fitrx_subject_string_ids(list(ebe_etas = data.frame(x = 1))))
  fit <- list(ebe_etas = data.frame(ID = c(1, 2, 10), eta = c(0, 0, 0)))
  expect_identical(ferx:::.fitrx_subject_string_ids(fit), c("1", "2", "10"))
})
test_that(".fitrx_sdtab_subject assigns rows from the per-subject counts (#540)", {
  # Gaussian rows subject by subject, then discrete rows subject by subject;
  # a subject with no rows of a kind takes none.
  expect_identical(
    ferx:::.fitrx_sdtab_subject(c(2L, 0L, 1L), c(1L, 1L, 0L), 5L),
    list(gaussian = c(1L, 1L, 3L), discrete = c(1L, 2L), tte = integer())
  )
  # Time-to-event rows come last, subject by subject (ferx-core #901).
  expect_identical(
    ferx:::.fitrx_sdtab_subject(c(2L, 0L, 1L), c(1L, 1L, 0L), 7L, c(1L, 0L, 1L)),
    list(gaussian = c(1L, 1L, 3L), discrete = c(1L, 2L), tte = c(1L, 3L))
  )
  expect_null(ferx:::.fitrx_sdtab_subject(c(2L, 0L, 1L), c(1L, 1L, 0L), 5L, c(1L, 0L, 1L)))
  # Not one TTE count per subject: assigns nothing, even when the sum fits.
  expect_null(ferx:::.fitrx_sdtab_subject(c(2L, 0L, 1L), c(1L, 1L, 0L), 7L, c(1L, 1L)))
  # Counts that do not add up to the rows, or are missing, assign nothing.
  expect_null(ferx:::.fitrx_sdtab_subject(c(2L, 1L), c(0L, 0L), 4L))
  expect_null(ferx:::.fitrx_sdtab_subject(c(2L, NA), c(0L, 0L), 2L))
})
test_that(".fitrx_write_ebes_csv writes a stable header when ebes are absent", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_ebes_csv(list(ebe_etas = NULL), path)
  expect_identical(readLines(path), "ID,ofv_contribution,n_obs")
})

# ---- header from test-sde.R ----
# Tests for SDE (diffusion / EKF) surface on the R side.
# All tests use make_fake_fit() from helper-trace.R — no real fit needed.

# uses_sde flag in print.ferx_fit ----------------------------------------





# uses_sde flag in summary / print.ferx_summary ---------------------------






# DIFF_* theta parameters flow through naturally --------------------------


# ferx_estimates() with DIFF_* theta -------------------------------------


# ferx_save_fit / ferx_load_fit round-trip --------------------------------




# ---- header from test-covtab.R ----
# Covariate table (`fit$covtab`) + types (`fit$covariate_types`), surfaced from
# a `[covariates]` block. See ferx-core #182/#184 for the underlying DSL + table.

# Cached fit on the covariate example so the multiple tests below share one
# estimation (mirrors helper-warfarin-fit.R's memoised pattern).
cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("two_cpt_oral_cov")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method = "focei", verbose = FALSE,
        covariance = FALSE, settings = list(maxiter = 5L)
      )
    }
    fit
  }
})








test_that("covtab and covariate_types survive a .fitrx save/load round-trip", {
  fit <- cov_fit()
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)

  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)

  expect_s3_class(fit2$covtab, "data.frame")
  expect_identical(names(fit2$covtab), names(fit$covtab))
  expect_equal(nrow(fit2$covtab), nrow(fit$covtab))
  expect_type(fit2$covtab$ID, "character")
  expect_identical(fit2$covariate_types, fit$covariate_types)
})

# ---- header from test-data-selection.R ----
# Tests for [data_selection] / ferx_apply_selection() / ferx_fit(ignore=) support.

# Local alias avoids the ::: operator.
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# .sel_attr()/warfarin_sel_fit() come from helper-selection.R

# ---------------------------------------------------------------------------
# 1. Plain warfarin fit has NULL exclusions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 2. [data_selection] model populates exclusions
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 3. ferx_fit(data, ignore = ...) produces exclusions matching model-file path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4. ferx_apply_selection() pure-R preview
# ---------------------------------------------------------------------------







# ---------------------------------------------------------------------------
# 5. ferx_data passed directly to ferx_fit()
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 6. print.ferx_fit() shows DATA SELECTION block
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 7. ferx_runlog() shows exclusion line
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 8. persist round-trip preserves exclusions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 9. print.ferx_data
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 10. ferx_apply_selection(excluded = TRUE) on a ferx_data
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 11. accept= and ignore_ids= coverage for ferx_apply_selection()
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 12. ferx_apply_selection() data.frame input
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 13. ferx_apply_selection(fit, excluded = TRUE) replays fired conditions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 14. ferx_apply_selection(fit, excluded = TRUE) NULL exclusions path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 15. ferx_fit() type validation stops for ignore/accept/ignore_ids
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 16. ferx_fit() accept= and ignore_ids= settings_parts code paths
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 17. print.ferx_fit fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 18. ferx_runlog fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 19. .cmp() operator coverage (pure R, via ferx_apply_selection() with inline df)
# ---------------------------------------------------------------------------

# .sel_test_df() comes from helper-selection.R





# ---------------------------------------------------------------------------
# 20. .parse_filter_expr() coverage: character rhs, &&, !found
# ---------------------------------------------------------------------------





test_that("ferx_save/ferx_load round-trip preserves exclusions", {
  skip_on_cran()
  fit  <- warfarin_sel_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$exclusions))
  expect_equal(loaded$exclusions$n_obs_excluded,   fit$exclusions$n_obs_excluded)
  expect_equal(loaded$exclusions$n_records_total,  fit$exclusions$n_records_total)
  expect_equal(loaded$exclusions$fired_ignore,     fit$exclusions$fired_ignore)
  expect_equal(loaded$exclusions$excluded_subject_ids, fit$exclusions$excluded_subject_ids)
})

# A synthetic fit whose four subject tables (`ebe_etas`, `covtab`,
# `ebe_kappas`, `cond_dist$data`) carry the two subject IDs `ids`.
id_trip_fake <- function(ids) {
  cd_data <- data.frame(
    ID = rep(ids, each = 2L),
    ETA = c("ETA_CL", "ETA_V", "ETA_CL", "ETA_V"),
    COND_MEAN = c(0.05, -0.02, -0.03, 0.01),
    COND_SD = c(0.12, 0.15, 0.11, 0.14),
    COND_MODE = c(0.04, -0.01, -0.02, 0.02),
    stringsAsFactors = FALSE
  )
  fake <- structure(
    list(
      theta = c(TVCL = 1.0, TVV = 10.0),
      omega = matrix(c(0.04, 0, 0, 0.09), 2L, 2L,
                     dimnames = list(c("ETA_CL", "ETA_V"), c("ETA_CL", "ETA_V"))),
      eta_names = c("ETA_CL", "ETA_V"),
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      cond_dist = list(data = cd_data, shrinkage = c(0.1, 0.05),
                        nsamp = 20L, burnin = 5L),
      ebe_etas = data.frame(ID = ids, ETA_CL = c(0.04, -0.02),
                             ETA_V = c(-0.01, 0.02), stringsAsFactors = FALSE),
      covtab = data.frame(ID = ids, TIME = c(0, 0), WT = c(70.5, 80.25),
                          stringsAsFactors = FALSE),
      ebe_kappas = data.frame(ID = rep(ids, each = 2L), OCC = c(1L, 2L, 1L, 2L),
                              KAPPA_CL = c(0.01, -0.01, 0.02, -0.02),
                              stringsAsFactors = FALSE),
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 2L, n_iterations = 1L,
      method = "SAEM", method_chain = "SAEM",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = c(0, 0), shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "Enzyme AD",
      gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_
    ),
    class = "ferx_fit"
  )
}

# --- #462: the fit's reader settings and population fingerprint ------------

test_that("reader settings and fingerprint survive a .fitrx round trip verbatim (#462)", {
  fit <- rs_fit(ignore = rs_ignore)
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge - skipping")
  for (include_data in c(FALSE, TRUE)) {
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path, include_data = include_data)
    loaded <- ferx_load_fit(path)
    info <- paste("include_data:", include_data)
    expect_identical(loaded$reader_settings, fit$reader_settings, info = info)
    expect_identical(loaded$population_fingerprint, fit$population_fingerprint,
                     info = info)
    # Written as the objects they are, under the engine's own keys.
    staging <- withr::local_tempdir()
    utils::unzip(path, exdir = staging)
    wire <- jsonlite::read_json(file.path(staging, "fit.json"),
                                simplifyVector = FALSE)
    expect_type(wire$reader_settings, "list")
    expect_type(wire$population_fingerprint, "list")
    out <- suppressWarnings(ferx_covariance(loaded))
    expect_lt(rs_rel(out$se_theta, fit$se_theta), 1e-8)
  }
})

test_that("reader settings with categorical levels survive a .fitrx round trip verbatim (#526)", {
  # serde writes the declared levels as `[0.0, 1.0]`; jsonlite reads them back
  # as integers and would write `[0, 1]`. The loader hands the field back to the
  # engine to spell, so the string is the fit's own.
  ex <- ferx_example("warfarin")
  rows <- utils::read.csv(ex$data)
  rows$SEX <- rows$ID %% 2
  data <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
  model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(c(readLines(ex$model), "", "[covariates]",
               "  SEX categorical(levels = [0, 1])"), model)
  fit <- suppressWarnings(ferx_fit(model, data, method = "focei", verbose = FALSE,
                                   covariance = FALSE, settings = list(maxiter = 3L)))
  expect_match(fit$reader_settings, "[0.0,1.0]", fixed = TRUE)
  for (include_data in c(FALSE, TRUE)) {
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path, include_data = include_data)
    loaded <- ferx_load_fit(path)
    info <- paste("include_data:", include_data)
    expect_identical(loaded$reader_settings, fit$reader_settings, info = info)
    expect_identical(loaded$population_fingerprint, fit$population_fingerprint,
                     info = info)
  }
})

test_that("ferx_save_fit refuses reader settings that are not JSON, by name (#526)", {
  fit <- rs_fit(ignore = rs_ignore)
  for (bad in list("{not json", NA_character_, c("{}", "{}"))) {
    edited <- fit
    edited <- editable_fit(edited)
    edited$reader_settings <- bad
    path <- withr::local_tempfile(fileext = ".fitrx")
    expect_error(ferx_save_fit(resealed_fit(edited), path), "`fit$reader_settings` is not the engine's JSON",
                 fixed = TRUE)
  }
  edited <- fit
  edited <- editable_fit(edited)
  edited$population_fingerprint <- "{not json"
  path <- withr::local_tempfile(fileext = ".fitrx")
  expect_error(ferx_save_fit(resealed_fit(edited), path),
               "`fit$population_fingerprint` is not the engine's JSON", fixed = TRUE)
})

test_that("ferx_save_fit() warns when the model file changed since the fit (#531 review)", {
  skip_on_cran()
  skip_if_not(exists("sha256sum", envir = asNamespace("tools")),
              "tools::sha256sum() needs R >= 4.5")
  ex <- ferx_example("warfarin")
  model <- withr::local_tempfile(fileext = ".ferx")
  file.copy(ex$model, model)
  fit <- suppressWarnings(ferx_fit(model, ex$data, method = "focei",
                                   covariance = FALSE, verbose = FALSE,
                                   settings = list(maxiter = 5L)))
  path <- withr::local_tempfile(fileext = ".fitrx")
  # Unchanged: quiet, and a reloaded fit (whose model is the bundle's copy)
  # re-saves quietly too.
  expect_no_warning(ferx_save_fit(fit, path), class = "ferx_model_edited_warning")
  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_no_warning(ferx_save_fit(loaded, path), class = "ferx_model_edited_warning")
  # One comment line appended after the fit. The fit carries the fitted
  # bytes (#577), so the bundle is written from them and there is nothing to
  # warn about.
  cat("# edited after the fit\n", file = model, append = TRUE)
  expect_no_warning(ferx_save_fit(fit, path), class = "ferx_model_edited_warning")
  # A fit from before the fitted bytes were carried falls back to the file:
  # the bundle is still written, with the file as it is now, and the warning
  # names it.
  fit <- editable_fit(fit)
  fit$model_text <- NULL
  fit <- resealed_fit(fit)
  expect_warning(ferx_save_fit(fit, path), class = "ferx_model_edited_warning")
  w <- tryCatch(ferx_save_fit(fit, path), ferx_model_edited_warning = conditionMessage)
  expect_match(w, basename(model), fixed = TRUE)
  expect_match(w, "no longer matches fit$model_hash", fixed = TRUE)
  expect_true(file.exists(path))
  expect_identical(fitrx_fit_json(path)$model_hash, fit$model_hash)
  # No recorded hash: nothing to compare against.
  nohash <- editable_fit(fit)
  nohash$model_hash <- NA_character_
  expect_no_warning(ferx_save_fit(resealed_fit(nohash), path),
                    class = "ferx_model_edited_warning")
})

test_that("ferx_save_fit needs no external zip program (#662)", {
  # A Linux box without Info-ZIP leaves R_ZIPCMD empty, and utils::zip()
  # refused with "'zip' must be a non-empty character string" - after the
  # whole fit, when it came through ferx_fit(output =). The bundle is now
  # written by the Rust glue, so an empty R_ZIPCMD changes nothing.
  withr::local_envvar(R_ZIPCMD = "")
  fit <- warfarin_fit_cov()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path, include_data = TRUE)

  names <- utils::unzip(path, list = TRUE)$Name
  staging <- withr::local_tempdir()
  utils::unzip(path, files = "manifest.json", exdir = staging)
  manifest <- jsonlite::read_json(file.path(staging, "manifest.json"),
                                  simplifyVector = TRUE)
  expect_setequal(names, manifest$entries)
  expect_false(any(grepl("/", names, fixed = TRUE)))

  loaded <- ferx_load_fit(path)
  expect_equal(loaded$ofv, fit$ofv)
  expect_equal(loaded$theta, fit$theta)
  # The engine's own reader opens it too, not just R's.
  core_path <- withr::local_tempfile(fileext = ".fitrx")
  expect_no_error(ferx:::ferx_rust_fitrx_engine_resave(path, core_path))

  # A bare path still lands at <path>.fitrx, not at Info-ZIP's <path>.zip.
  bare <- file.path(withr::local_tempdir(), "run1")
  ferx_save_fit(fit, bare)
  expect_true(file.exists(paste0(bare, ".fitrx")))
  expect_false(file.exists(paste0(bare, ".zip")))
})
