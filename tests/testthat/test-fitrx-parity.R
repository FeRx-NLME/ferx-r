# The parity gate of ferx-core #2693 Phases 3 and 4 (#711, #712): the list of a fit
# equals the list from saving and loading that fit, field by field, because the
# bytes are what is saved and both lists come out of one converter. The fits cover
# what a bundle has a separate path for. Paths and the bundle bytes differ by
# design: a load stages files, and a save adds the data.
fitrx_parity_skip <- c("bundle", "model_path", "data_path", "trace_path")

# The fields of `fit` and `loaded` that are not identical().
fitrx_parity_diff <- function(fit, loaded, skip = fitrx_parity_skip) {
  a <- unclass(fit)
  b <- unclass(loaded)
  keys <- setdiff(union(names(a), names(b)), skip)
  keys[!vapply(keys, function(k) identical(a[[k]], b[[k]]), logical(1L))]
}

# The entries of a bundle, from its path or its bytes.
fitrx_entries <- function(bundle) {
  if (is.raw(bundle)) {
    path <- tempfile(fileext = ".fitrx")
    writeBin(bundle, path)
    bundle <- path
  }
  utils::unzip(bundle, list = TRUE)$Name
}

fitrx_round_trip <- function(fit, include_data = TRUE) {
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path, include_data = include_data)
  ferx_load_fit(path)
}

fitrx_parity_fit <- function(model, data, ...) {
  suppressMessages(suppressWarnings(ferx_fit(model, data, verbose = FALSE, ...)))
}

fitrx_parity_cases <- list(
  plain = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, ex$data, method = "foce")
  },
  iov = function() {
    ex <- ferx_example("warfarin_iov")
    fitrx_parity_fit(ex$model, ex$data, method = "focei", covariance = FALSE)
  },
  nested_levels = function() {
    ex <- ferx_example("nested_levels")
    fitrx_parity_fit(ex$model, ex$data, method = "focei", covariance = FALSE,
                     settings = list(maxiter = 5))
  },
  sir = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, ex$data, method = "focei", sir = TRUE,
                     settings = list(sir_samples = 200, sir_resamples = 50,
                                     sir_keep_samples = TRUE))
  },
  bayes = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, ex$data, method = "bayes", covariance = FALSE,
                     settings = list(bayes_warmup = 20L, bayes_iters = 40L,
                                     bayes_chains = 2L, bayes_seed = 1L))
  },
  vi = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, ex$data, method = "vi", covariance = FALSE,
                     settings = list(vi_iters = 30, vi_seed = 42))
  },
  mixture = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(test_path("fixtures", "mixture_override.ferx"), ex$data,
                     covariance = FALSE, settings = list(maxiter = 5))
  },
  dcm = function() {
    ex <- ferx_example("warfarin_dcm")
    fitrx_parity_fit(ex$model, ex$data, covariance = FALSE,
                     settings = list(maxiter = 2))
  },
  binary = function() {
    ex <- ferx_example("binary_logistic")
    fitrx_parity_fit(ex$model, ex$data, covariance = FALSE)
  },
  tte = function() {
    ex <- ferx_example("tte_weibull")
    fitrx_parity_fit(ex$model, ex$data, covariance = FALSE)
  },
  pk_tte = function() {
    ex <- ferx_example("pktte_joint")
    fitrx_parity_fit(ex$model, ex$data, covariance = FALSE)
  },
  saem_conddist = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, ex$data, method = "saem", covariance = FALSE,
                     settings = list(conddist = TRUE))
  },
  trace = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, ex$data, method = "foce", covariance = FALSE,
                     settings = list(optimizer_trace = TRUE))
  },
  data_selection = function() {
    ex <- ferx_example("warfarin_data_selection")
    fitrx_parity_fit(ex$model, ex$data, covariance = FALSE)
  },
  weighted_iov = function() {
    ex <- ferx_example("mbma_placebo")
    fitrx_parity_fit(ex$model, ex$data, covariance = FALSE,
                     settings = list(maxiter = 5))
  },
  text_ids = function() {
    ex <- ferx_example("warfarin")
    d <- read.csv(ex$data, colClasses = c(ID = "character"))
    ids <- unique(d$ID)
    labels <- c("NA", "007", "Smith, 2019", "say \"hi\"", "1e3", "two\nlines")
    d$ID <- c(labels, ids)[match(d$ID, ids)]
    fitrx_parity_fit(ex$model, d, method = "foce", covariance = FALSE)
  },
  data_frame = function() {
    ex <- ferx_example("warfarin")
    fitrx_parity_fit(ex$model, read.csv(ex$data), method = "foce",
                     covariance = FALSE)
  }
)

for (case in names(fitrx_parity_cases)) {
  test_that(sprintf("a %s fit's list survives save and load (#711)", case), {
    skip_on_cran()
    fit <- fitrx_parity_cases[[case]]()
    expect_true(is.raw(fit$bundle))
    expect_identical(fitrx_parity_diff(fit, fitrx_round_trip(fit)), character())
  })
}

test_that("a bundle without data.csv loads, without the data-dependent field (#711)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$iov()
  loaded <- fitrx_round_trip(fit, include_data = FALSE)
  expect_null(loaded$individual_estimates)
  expect_false(is.null(fit$individual_estimates))
  expect_identical(
    fitrx_parity_diff(fit, loaded, c(fitrx_parity_skip, "individual_estimates")),
    character()
  )
})

test_that("what the post-hoc steps add is saved, and a re-save is stable (#711)", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  fit <- fitrx_parity_fit(ex$model, ex$data, method = "foce", covariance = FALSE)
  with_cov <- suppressWarnings(ferx_covariance(fit))
  loaded <- fitrx_round_trip(with_cov)
  expect_identical(fitrx_parity_diff(with_cov, loaded), character())
  expect_identical(fitrx_parity_diff(loaded, fitrx_round_trip(loaded)), character())
})

test_that("the bytes on a fresh fit are a bundle core loads into the same fit (#711)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  raw <- ferx:::ferx_rust_load_fit(fit$bundle)$fit
  expect_identical(raw$ofv, fit$ofv)
  expect_identical(raw$theta, unname(fit$theta))
  expect_identical(raw$subject_ids, fit$subject_ids)
})

test_that("a bundle of the old format fails with a version error (#711)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  staging <- tempfile("fitrx_v1_")
  dir.create(staging)
  utils::unzip(path, exdir = staging)
  manifest <- file.path(staging, "manifest.json")
  writeLines(sub('"format_version": "2"', '"format_version": "1"', readLines(manifest)),
             manifest)
  old <- tempfile(fileext = ".fitrx")
  utils::zip(old, list.files(staging, full.names = TRUE), flags = "-j -q")
  expect_error(ferx_load_fit(old), "format_version")
})

test_that("a model the engine cannot read leaves cov_matrix unnamed (#417)", {
  bad <- withr::local_tempfile(fileext = ".ferx")
  writeLines("[not_a_block]", bad)
  raw <- list(cov_matrix_dim = 2L, model_path = bad, model_hash = "")
  expect_identical(ferx:::.fitrx_cov_labels(raw, NULL), character())
})

test_that("a subject ID holding a newline reads back from both writers (#711)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$text_ids()
  expect_true("two\nlines" %in% fit$subject_ids)
  raw <- ferx:::ferx_rust_load_fit(fit$bundle)$fit
  expect_identical(raw$subject_ids, fit$subject_ids)
  expect_identical(raw$ebe_etas$ID, fit$ebe_etas$ID)
})

# -- Phase 4 (#712): the bytes are the fit --

# Fails with the fields where `x`'s list and its bytes disagree.
expect_bytes_agree <- function(x, label) {
  expect_identical(fitrx_parity_diff(x, fitrx_round_trip(x)), character(), info = label)
}

fitrx_steps <- list(
  cov = function(f) suppressWarnings(ferx_covariance(f)),
  sir = function(f) suppressWarnings(ferx_sir(
    f, sir_samples = 50L, sir_resamples = 20L, sir_seed = 1L, sir_keep_samples = TRUE
  )),
  npde = function(f) suppressWarnings(ferx_calc_npde(f, nsim = 20L, seed = 1L)),
  trip = function(f) suppressWarnings(fitrx_round_trip(f))
)

test_that("save, load and post-hoc steps in every order keep bytes and list in agreement (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  orders <- list(
    c("cov", "sir", "npde"), c("cov", "npde", "sir"), c("npde", "cov", "sir"),
    c("sir", "cov", "npde"), c("sir", "npde", "cov"), c("npde", "sir", "cov"),
    c("trip", "cov", "trip", "sir", "trip", "npde", "trip"),
    c("npde", "trip", "sir", "trip", "cov", "trip")
  )
  for (order in orders) {
    x <- fit
    for (i in seq_along(order)) {
      x <- fitrx_steps[[order[[i]]]](x)
      expect_bytes_agree(x, paste(order[seq_len(i)], collapse = " > "))
    }
    # What each step added is there at the end, whatever came after it.
    expect_false(is.null(x$cov_matrix))
    expect_false(is.null(x$sir_ess))
    expect_true(all(c("NPDE", "NPD") %in% names(x$sdtab)))
    expect_false(is.null(x$sir_resamples))
  }
})

test_that("a step on a fit from the middle of a chain does not see what came after (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  with_npde <- fitrx_steps$npde(fit)
  expect_false("NPDE" %in% names(fit$sdtab))
  expect_false("NPDE" %in% names(fitrx_round_trip(fit)$sdtab))
  expect_true("NPDE" %in% names(fitrx_round_trip(with_npde)$sdtab))
  # ... and the columns reach core's own reader, not only the glue's.
  raw <- ferx:::ferx_rust_load_fit(with_npde$bundle)$fit
  expect_true("NPDE" %in% names(raw$sdtab))
})

test_that("the post-hoc steps leave their input fit as it was (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  bytes <- as.raw(as.integer(fit$bundle))
  snapshot <- unserialize(serialize(fit, NULL))
  for (step in c("cov", "sir", "npde")) {
    out <- fitrx_steps[[step]](fit)
    expect_identical(fit, snapshot, info = step)
    expect_identical(fit$bundle, bytes, info = step)
    expect_false(identical(out$bundle, bytes), info = step)
  }
})

test_that("ferx_covariance() on a fit without one adds it to the bytes, not to the input (#712)", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  fit <- fitrx_parity_fit(ex$model, ex$data, method = "foce", covariance = FALSE)
  out <- suppressWarnings(ferx_covariance(fit))
  expect_null(fit$cov_matrix)
  expect_null(fitrx_round_trip(fit)$cov_matrix)
  expect_false(is.null(out$cov_matrix))
  expect_identical(fitrx_round_trip(out)$cov_matrix, out$cov_matrix)
})

for (case in c("iov", "binary", "pk_tte", "data_frame", "text_ids")) {
  test_that(sprintf("a %s fit's list and bytes agree after ferx_covariance() (#712)", case), {
    skip_on_cran()
    fit <- fitrx_parity_cases[[case]]()
    out <- suppressWarnings(ferx_covariance(fit))
    expect_identical(out$sdtab, fit$sdtab)
    expect_identical(out$subject_ids, fit$subject_ids)
    expect_identical(out$data_source, fit$data_source)
    expect_bytes_agree(out, case)
  })
}

test_that("a step does not lose what the bundle cannot rebuild (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  out <- suppressWarnings(ferx_covariance(fit))
  expect_identical(out$individual_estimates, fit$individual_estimates)
  expect_identical(out$model_path, fit$model_path)
  expect_identical(out$data_path, fit$data_path)
  expect_identical(out$call_settings, fit$call_settings)
  expect_identical(out$model_name, fit$model_name)
})

test_that("a fit with no bytes is refused everywhere, saying to refit, with no fallback (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  fit$bundle <- NULL
  # The way such a fit arrives: an old `.rds`.
  rds <- withr::local_tempfile(fileext = ".rds")
  saveRDS(fit, rds)
  old <- readRDS(rds)
  out <- withr::local_tempfile(fileext = ".fitrx")
  for (call in list(
    function() ferx_save_fit(old, out),
    function() ferx_covariance(old),
    function() ferx_sir(old),
    function() ferx_calc_npde(old)
  )) {
    err <- expect_error(call(), "Refit", class = "ferx_no_bundle_error")
    expect_match(conditionMessage(err), "Old fits are not supported", fixed = TRUE)
  }
  expect_false(file.exists(out))
})

test_that("a hand edit of the list is not saved: the bytes are the fit (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  edited <- fit
  edited$theta[[1L]] <- 99
  expect_identical(fitrx_round_trip(edited)$theta, fit$theta)
})

test_that("a hand edit of the list reaches no step: each reads the bytes (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()

  # SIR: an edited record is not run; an explicit argument still is.
  sir_args <- list(sir_samples = 50L, sir_resamples = 20L, sir_seed = 1L)
  clean <- suppressWarnings(do.call(ferx_sir, c(list(fit), sir_args)))
  edited <- clean
  edited$sir_settings$df <- 1
  edited$sir_settings$samples <- 60L
  again <- suppressWarnings(ferx_sir(edited))
  expect_identical(again$sir_settings, clean$sir_settings)
  expect_identical(again$sir_ci_theta, clean$sir_ci_theta)

  # Covariance: the call's step settings and the scoring record are the bytes'.
  cov <- suppressWarnings(ferx_covariance(fit))
  edited <- fit
  edited$call_settings$covariance_method <- "s"
  out <- suppressWarnings(ferx_covariance(edited))
  expect_identical(out$covariance_method, cov$covariance_method)
  expect_identical(out$cov_matrix, cov$cov_matrix)
  edited <- fit
  edited$scoring_settings$inner_maxiter <- 1L
  expect_identical(
    suppressWarnings(ferx_covariance(edited, mu_referencing = TRUE))$cov_matrix,
    suppressWarnings(ferx_covariance(fit, mu_referencing = TRUE))$cov_matrix
  )

  # NPDE: computed at the fit's estimates, on the fit's table.
  clean <- fitrx_steps$npde(fit)
  edited <- fit
  edited$theta[[1L]] <- edited$theta[[1L]] * 3
  edited$sdtab$PRED <- 0
  out <- fitrx_steps$npde(edited)
  expect_identical(out$theta, fit$theta)
  expect_identical(out$sdtab, clean$sdtab)
  expect_bytes_agree(out, "npde on an edited list")
})

test_that("write_model renders the bytes, not an edited list (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  dir <- withr::local_tempdir()
  ferx_save_fit(fit, file.path(dir, "clean.fitrx"), write_model = TRUE)
  edited <- fit
  edited$theta[[1L]] <- 99
  edited$omega[1L, 1L] <- 9
  edited$sigma[[1L]] <- 9
  edited$model_text <- "not the model"
  ferx_save_fit(edited, file.path(dir, "edited.fitrx"), write_model = TRUE)
  expect_identical(readLines(file.path(dir, "edited_final.ferx")),
                   readLines(file.path(dir, "clean_final.ferx")))
})

test_that("a save whose model file cannot be written leaves the bundle (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  dir <- withr::local_tempdir()
  out <- file.path(dir, "run.fitrx")
  ferx_save_fit(fit, out)
  before <- readBin(out, "raw", n = file.size(out))
  dir.create(file.path(dir, "run_final.ferx"))
  later <- suppressWarnings(ferx_covariance(fit))
  expect_error(ferx_save_fit(later, out, write_model = TRUE), "is a directory")
  expect_identical(readBin(out, "raw", n = file.size(out)), before)
  expect_setequal(list.files(dir), c("run.fitrx", "run_final.ferx"))
})

test_that("a save replaces an existing bundle and model file, keeping their mode (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  later <- suppressWarnings(ferx_covariance(fit))
  dir <- withr::local_tempdir()
  out <- file.path(dir, "run.fitrx")
  model <- file.path(dir, "run_final.ferx")
  ferx_save_fit(fit, out, write_model = TRUE)
  cat("# stale\n", file = model)
  Sys.chmod(c(out, model), "600", use_umask = FALSE)

  ferx_save_fit(later, out, write_model = TRUE)

  expect_identical(ferx_load_fit(out)$cov_matrix, later$cov_matrix)
  expect_false(any(grepl("# stale", readLines(model), fixed = TRUE)))
  expect_setequal(list.files(dir, all.files = TRUE, no.. = TRUE),
                  c("run.fitrx", "run_final.ferx"))
  skip_on_os("windows")
  expect_identical(as.character(file.mode(c(out, model))), c("600", "600"))
})

test_that("a failed save leaves what was at the path, and no stray file (#712)", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  dir <- withr::local_tempdir()
  data <- file.path(dir, "data.csv")
  file.copy(ex$data, data)
  fit <- fitrx_parity_fit(ex$model, data, method = "foce", covariance = FALSE)
  out <- file.path(dir, "run.fitrx")
  ferx_save_fit(fit, out)
  before <- readBin(out, "raw", n = file.size(out))
  cat("1,0,0\n", file = data, append = TRUE)
  expect_error(ferx_save_fit(fit, out, include_data = TRUE),
               class = "ferx_data_edited_error")
  expect_identical(readBin(out, "raw", n = file.size(out)), before)
  expect_setequal(list.files(dir), c("data.csv", "run.fitrx"))
})

test_that("a re-save leaves a bundled copy of the data out unless asked (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  loaded <- fitrx_round_trip(fit, include_data = TRUE)
  expect_true("data.csv" %in% fitrx_entries(loaded$bundle))
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(loaded, path)
  expect_false("data.csv" %in% fitrx_entries(path))
  manifest <- jsonlite::fromJSON(unz(path, "manifest.json"))
  expect_setequal(manifest$entries, fitrx_entries(path))
})

test_that("an untraced job leaves no trace in the bytes (#712)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$trace()
  expect_false(is.null(fit$trace))
  bytes <- ferx:::ferx_rust_bundle_without_trace(fit$bundle)
  expect_false("trace.csv" %in% fitrx_entries(bytes))
  out <- ferx:::.ferx_fit_with_bundle(fit, bytes)
  expect_null(out$trace)
  expect_null(out$trace_path)
  expect_identical(out$theta, fit$theta)
  expect_identical(out$sdtab, fit$sdtab)
  expect_bytes_agree(out, "untraced")
})

test_that("a rewrite does not write through a planted temp path (#712)", {
  skip_on_cran()
  skip_on_os("windows")
  fit <- fitrx_parity_cases$trace()
  dir <- withr::local_tempdir()
  victim <- withr::local_tempfile(lines = "untouched")
  # Every name of the predictable scheme this process could have reached.
  planted <- file.path(dir, sprintf("ferx_bundle_%d_%d.fitrx", Sys.getpid(), 0:4999))
  expect_true(all(file.symlink(victim, planted)))
  withr::local_envvar(TMPDIR = dir)
  bytes <- ferx:::ferx_rust_bundle_without_trace(fit$bundle)
  expect_false("trace.csv" %in% fitrx_entries(bytes))
  expect_identical(readLines(victim), "untouched")
  # The private directory goes with the call.
  expect_setequal(list.files(dir, all.files = TRUE, no.. = TRUE), basename(planted))
})

# The CLI is ferx-core's, so CI has it only when pointed at one: FERX_CLI, else
# `ferx` on the PATH.
fitrx_cli <- function() {
  cli <- Sys.getenv("FERX_CLI", unname(Sys.which("ferx")))
  if (!nzchar(cli) || !file.exists(cli)) skip("no ferx CLI (set FERX_CLI)")
  cli
}

test_that("a bundle saved from R loads in the CLI (#712)", {
  skip_on_cran()
  cli <- fitrx_cli()
  for (case in c("plain", "iov", "binary", "pk_tte", "trace")) {
    fit <- fitrx_parity_cases[[case]]()
    path <- tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    out <- suppressWarnings(system2(cli, c("summary", shQuote(path)),
                                    stdout = TRUE, stderr = TRUE))
    expect_null(attr(out, "status"), info = case)
    expect_true(any(grepl(names(fit$theta)[[1L]], out, fixed = TRUE)), info = case)
  }
})

test_that("a bundle a post-hoc step wrote loads in the CLI (#712)", {
  skip_on_cran()
  cli <- fitrx_cli()
  fit <- fitrx_parity_cases$plain()
  x <- fitrx_steps$npde(fitrx_steps$sir(fitrx_steps$cov(fit)))
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(x, path, include_data = TRUE)
  out <- suppressWarnings(system2(cli, c("summary", shQuote(path)),
                                  stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"))
})

test_that("a bundle the CLI wrote loads in R, and saves back (#712)", {
  skip_on_cran()
  cli <- fitrx_cli()
  ex <- ferx_example("warfarin")
  dir <- withr::local_tempdir()
  run <- file.path(dir, "run.fitrx")
  withr::with_dir(dir, {
    log <- suppressWarnings(system2(
      cli, c(shQuote(ex$model), "--data", shQuote(ex$data), "--output", shQuote(run),
             "--include-data"),
      stdout = TRUE, stderr = TRUE
    ))
  })
  expect_null(attr(log, "status"))
  fit <- ferx_load_fit(run)
  expect_true(is.finite(fit$ofv))
  expect_identical(fit$bundle, readBin(run, "raw", n = file.size(run)))
  # The steps run on it, and the result goes back to the CLI.
  out <- suppressWarnings(ferx_covariance(fit))
  expect_false(is.null(out$cov_matrix))
  again <- file.path(dir, "again.fitrx")
  ferx_save_fit(out, again, include_data = TRUE)
  summary <- suppressWarnings(system2(cli, c("summary", shQuote(again)),
                                      stdout = TRUE, stderr = TRUE))
  expect_null(attr(summary, "status"))
})
