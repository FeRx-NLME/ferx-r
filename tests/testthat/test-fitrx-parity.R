# The parity gate of ferx-core #2693 Phase 3 (#711): the list of a fit equals the
# list from saving and loading that fit, field by field, because both come out of
# one converter. The fits cover what a bundle has a separate path for. Paths and
# the bundle bytes differ by design: a load stages files, and the R writer's
# archive is not the engine's until #712.
fitrx_parity_skip <- c("bundle", "model_path", "data_path", "trace_path")

# The fields of `fit` and `loaded` that are not identical().
fitrx_parity_diff <- function(fit, loaded, skip = fitrx_parity_skip) {
  a <- unclass(fit)
  b <- unclass(loaded)
  keys <- setdiff(union(names(a), names(b)), skip)
  keys[!vapply(keys, function(k) identical(a[[k]], b[[k]]), logical(1L))]
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
  ferx:::ferx_rust_zip_files(staging, list.files(staging), old)
  expect_error(ferx_load_fit(old), "format_version")
})

test_that("a fit that carries no bundle is refused by ferx_save_fit() (#711)", {
  skip_on_cran()
  fit <- fitrx_parity_cases$plain()
  fit$bundle <- NULL
  expect_error(ferx_save_fit(fit, tempfile(fileext = ".fitrx")), "Refit")
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
