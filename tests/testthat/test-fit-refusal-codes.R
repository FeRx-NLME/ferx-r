# Refusals that only the engine's fit() used to raise, without a code, carry
# one since ferx-core #1788 (#2321), and `ferx_model_validate()` reports the
# model-level ones, which it used to call valid. They reach R as `e$code`
# through the glue's engine_refusal().

rc_write <- function(txt) {
  p <- withr::local_tempfile(fileext = ".ferx", .local_envir = parent.frame())
  writeLines(txt, p)
  p
}

# warfarin with a dedicated residual-error eta (`iiv_on_ruv = <eta>`) under
# `method`.
rc_iiv_on_ruv <- function(eta, method) {
  t <- readLines(ferx_example("warfarin")$model)
  omegas <- grep("^  omega ", t)
  t <- append(t, "  omega ETA_RUV ~ 0.05", omegas[length(omegas)])
  t <- append(t, paste0("  iiv_on_ruv = ", eta), grep("^\\[error_model\\]", t) + 1L)
  sub("^  method *=.*", paste0("  method = ", method), t)
}

rc_code <- function(expr) {
  e <- tryCatch({
    suppressWarnings(suppressMessages(expr))
    NULL
  }, error = function(e) e)
  expect_false(is.null(e), label = "the refusal")
  e$code
}

test_that("iiv_on_ruv refusals carry their codes (ferx-core #1788)", {
  skip_on_cran()
  data <- ferx_example("warfarin")$data
  foce <- rc_write(rc_iiv_on_ruv("ETA_RUV", "foce"))
  expect_identical(rc_code(ferx_fit(foce, data, verbose = FALSE)), "E_IIV_ON_RUV_METHOD")
  shared <- rc_write(rc_iiv_on_ruv("ETA_CL", "focei"))
  expect_identical(rc_code(ferx_fit(shared, data, verbose = FALSE)),
                   "E_IIV_ON_RUV_SHARED_ETA")
  # The validator reports the method refusal, which it used to call valid.
  out <- paste(capture.output(suppressWarnings(ferx_model_validate(foce))), collapse = "\n")
  expect_match(out, "ERROR E_IIV_ON_RUV_METHOD", fixed = TRUE)
})

test_that("a derived iov_occasion with one occasion per subject is E_IOV_SINGLE_OCCASION", {
  skip_on_cran()
  t <- readLines(ferx_example("warfarin")$model)
  omegas <- grep("^  omega ", t)
  t <- append(t, "  kappa KAPPA_CL ~ 0.04", omegas[length(omegas)])
  t <- sub("exp(ETA_CL)", "exp(ETA_CL + KAPPA_CL)", t, fixed = TRUE)
  model <- rc_write(t)
  expect_identical(
    rc_code(ferx_fit(model, ferx_example("warfarin")$data, verbose = FALSE,
                     settings = list(iov_occasion = "time(1000)"))),
    "E_IOV_SINGLE_OCCASION"
  )
})
