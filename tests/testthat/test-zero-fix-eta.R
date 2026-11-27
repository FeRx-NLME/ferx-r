# ferx-core #3150: an `omega NAME ~ 0 FIX` placeholder eta is removed from the
# model, as NONMEM holds `$OMEGA 0 FIX`. The fit is the model written without
# it, bit for bit; the eta has no row in fit$omega, fit$eta_names, the EBEs, the
# shrinkage tables or the sdtab; and anything supplied for the declared shape is
# refused.

zf_scratch <- function(model_text) {
  path <- tempfile(fileext = ".ferx")
  writeLines(model_text, path)
  path
}

# The warfarin example with a placeholder eta on KA.
zf_models <- function() {
  ex <- ferx_example("warfarin")
  base <- readLines(ex$model)
  with_zero <- sub("omega ETA_KA ~ 0.40", "omega ETA_KA ~ 0.40\n  omega ETA_Z ~ 0 FIX",
                   base, fixed = TRUE)
  with_zero <- sub("KA = TVKA * exp(ETA_KA)", "KA = TVKA * exp(ETA_KA + ETA_Z)",
                   with_zero, fixed = TRUE)
  stopifnot(!identical(with_zero, base))
  list(data = ex$data, base = zf_scratch(base), zero = zf_scratch(with_zero))
}

test_that("a 0 FIX eta is removed from the fit, which equals the model without it", {
  skip_on_cran()
  m <- zf_models()
  fit_r <- function(model) {
    suppressWarnings(ferx_fit(model, m$data, covariance = FALSE, verbose = FALSE))
  }
  zero <- fit_r(m$zero)
  base <- fit_r(m$base)

  expect_identical(zero$eta_names, c("ETA_CL", "ETA_V", "ETA_KA"))
  expect_identical(dimnames(zero$omega), list(zero$eta_names, zero$eta_names))
  expect_identical(names(zero$shrinkage_eta), zero$eta_names)
  expect_false("ETA_Z" %in% names(zero$sdtab))
  expect_false("ETA_Z" %in% names(zero$ebe_etas))
  # Bit for bit the model written without the placeholder.
  expect_identical(zero$ofv, base$ofv)
  expect_identical(zero$theta, base$theta)
  expect_identical(zero$omega, base$omega)

  # The engine's note, with ferx_get_warnings() guidance under it.
  w <- ferx_get_warnings(zero, as_df = TRUE)
  row <- w[grepl("^W_ZERO_FIX_ETA_DROPPED", w$message), , drop = FALSE]
  expect_identical(nrow(row), 1L)
  expect_match(row$message, "omega ETA_Z ~ 0 FIX", fixed = TRUE)
  g <- ferx:::.ferx_warning_guidance(row$category, message = row$message)
  expect_match(g, "Informational", fixed = TRUE)
  expect_match(g, "must leave it out", fixed = TRUE)
  expect_false(any(grepl("ZERO_FIX", ferx_get_warnings(base, as_df = TRUE)$message)))

  # The removed eta stays removed through a .fitrx round trip.
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(zero, path)
  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_identical(loaded$eta_names, zero$eta_names)
  expect_identical(dim(loaded$omega), c(3L, 3L))
})

test_that("etas supplied for the declared shape are refused with E_SUPPLIED_ETAS", {
  skip_on_cran()
  ex <- ferx_example("adaptive_tdm")
  text <- readLines(ex$model)
  text <- sub("omega ETA_CL ~ 0.09", "omega ETA_CL ~ 0.09\n  omega ETA_Z ~ 0 FIX",
              text, fixed = TRUE)
  text <- sub("V  = TVV", "V  = TVV * exp(ETA_Z)", text, fixed = TRUE)
  model <- zf_scratch(text)
  ids <- unique(utils::read.csv(ex$data)$ID)

  # Sized for the declared shape (ETA_CL, ETA_Z): refused, quoting the note.
  declared <- data.frame(ID = ids, ETA_CL = 0.1, ETA_Z = 0)
  err <- tryCatch(
    suppressWarnings(ferx_simulate_adaptive(model, ex$data, n_sim = 1L, seed = 1L,
                                            etas = declared)),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_identical(err$code, "E_SUPPLIED_ETAS")
  expect_match(conditionMessage(err), "W_ZERO_FIX_ETA_DROPPED", fixed = TRUE)

  # Sized for the model as fitted (ETA_CL alone): accepted.
  fitted_shape <- data.frame(ID = ids, ETA_CL = 0.1)
  res <- suppressWarnings(ferx_simulate_adaptive(model, ex$data, n_sim = 1L, seed = 1L,
                                                 etas = fitted_shape))
  expect_true(all(res$metrics$ETA_SUPPLIED))
})
