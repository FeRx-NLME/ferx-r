# print.ferx_fit()'s OMEGA_IOV rows are ferx-core's own console rows (#470):
# format_kappa_rows_from over the R fit's fields, through
# ferx_rust_kappa_rows(). The oracle never touches that glue: the fit is saved,
# ferx-core's .fitrx loader rebuilds a FitResult, and the rows come from
# KappaRowsInput::from_result - the path the engine's console takes
# (ferx_rust_fitrx_kappa_rows). Any field the glue drops, misreads or
# reorders shows up as a difference between the two.

kappa_rows <- getFromNamespace(".ferx_kappa_rows", "ferx")
fitrx_kappa_rows <- getFromNamespace("ferx_rust_fitrx_kappa_rows", "ferx")

# The engine console's rows for `fit`, via a saved bundle.
core_rows <- function(fit) {
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  strsplit(fitrx_kappa_rows(f), "\n", fixed = TRUE)[[1L]]
}

# The OMEGA_IOV block of print(fit): the lines between its rule and the next
# blank line.
printed_kappa_rows <- function(fit) {
  out <- utils::capture.output(print(fit))
  start <- grep("OMEGA_IOV", out, fixed = TRUE)
  expect_length(start, 1L)
  body <- out[-seq_len(start + 1L)]
  body[seq_len(match(TRUE, body == "" | startsWith(body, "  ---"), nomatch = length(body) + 1L) - 1L)]
}

# Both sides, UTF-8 glyphs, plus print() agreeing with the helper it calls.
expect_core_rows <- function(fit, info) {
  core <- core_rows(fit)
  expect_gt(length(core), 0L)
  expect_identical(kappa_rows(fit, ascii = FALSE), core, info = info)
  expect_identical(printed_kappa_rows(fit), kappa_rows(fit), info = info)
  invisible(core)
}

mbma_variant <- local({
  cache <- list()
  function(name, edit = identity) {
    if (is.null(cache[[name]])) {
      ex <- ferx_example("mbma_placebo")
      txt <- readLines(ex$model)
      out <- edit(txt)
      mod <- tempfile(fileext = ".ferx")
      writeLines(out, mod)
      cache[[name]] <<- suppressWarnings(ferx_fit(mod, ex$data, verbose = FALSE))
    }
    cache[[name]]
  }
})

kappa_line <- "  kappa KAPPA_ARM ~ 4.0 (sd) weight = NARM"
edit_kappa <- function(new) function(txt) {
  stopifnot(sum(txt == kappa_line) == 1L)
  txt[txt == kappa_line] <- new
  txt
}

test_that("mbma_placebo (weighted kappa): print shows ferx-core's rows", {
  skip_on_cran()
  fit <- mbma_variant("weighted")
  skip_if(!identical(fit$covariance_status, "computed"), "covariance step did not run")
  core <- expect_core_rows(fit, "weighted")
  var <- fit$omega_iov[1L, 1L]
  n <- fit$kappa_weight_typical[["KAPPA_ARM"]]
  # Not a self-comparison of two empty strings: the row and the weight line,
  # with the SE and the typical-arm SD.
  expect_length(core, 2L)
  expect_match(core[1L], sprintf("(SD = %.4f at weight 1)  SE = %.6f",
                                 sqrt(var), fit$se_kappa[[1L]]), fixed = TRUE)
  expect_match(core[2L], sprintf("weight = NARM  →  SD = %.4f at NARM = %.4f",
                                 sqrt(var / n), n), fixed = TRUE)
})

test_that("mbma_placebo, kappa unweighted: print shows ferx-core's rows", {
  skip_on_cran()
  fit <- mbma_variant("unweighted", edit_kappa("  kappa KAPPA_ARM ~ 4.0 (sd)"))
  expect_null(fit$kappa_weights)
  core <- expect_core_rows(fit, "unweighted")
  expect_length(core, 1L)
  expect_match(core, sprintf("(SD = %.4f)", sqrt(fit$omega_iov[1L, 1L])), fixed = TRUE)
})

test_that("mbma_placebo, kappa FIX: print shows ferx-core's rows", {
  skip_on_cran()
  fit <- mbma_variant("fix", edit_kappa("  kappa KAPPA_ARM ~ 4.0 (sd) FIX weight = NARM"))
  expect_identical(unname(fit$kappa_fixed), TRUE)
  core <- expect_core_rows(fit, "fix")
  expect_match(core[1L], "^  KAPPA_ARM \\[FIX\\] +=")
  expect_match(core[1L], "SE = ---$")
})

test_that("a failed or SIR-fallback covariance step hides the SD on both lines", {
  skip_on_cran()
  fit <- mbma_variant("weighted")
  for (status in c("failed", "sir_fallback", "computed", "not_requested")) {
    f <- fit
    f$covariance_status <- status
    core <- expect_core_rows(f, status)
    shown <- !status %in% c("failed", "sir_fallback")
    # Both sides of the gate, on the row and on the weight line.
    expect_identical(grepl("(SD = ", core[1L], fixed = TRUE), shown, info = status)
    expect_identical(grepl("SD = ", core[2L], fixed = TRUE), shown, info = status)
  }
})

test_that("a reloaded .fitrx prints the same ferx-core rows", {
  skip_on_cran()
  fit <- mbma_variant("weighted")
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  loaded <- suppressWarnings(ferx_load_fit(f))
  expect_identical(kappa_rows(loaded, ascii = FALSE),
                   strsplit(fitrx_kappa_rows(f), "\n", fixed = TRUE)[[1L]])
  expect_identical(printed_kappa_rows(loaded), printed_kappa_rows(fit))
})

test_that("block_kappa: one diagonal SE per kappa, ferx-core's rows", {
  skip_on_cran()
  ex <- ferx_example("warfarin_iov")
  txt <- readLines(ex$model)
  txt <- sub("^\\s*kappa KAPPA_CL ~ .*$",
             "  block_kappa (KAPPA_CL, KAPPA_V) = [0.04, 0.005, 0.02]", txt)
  txt <- sub("TVV  \\* exp\\(ETA_V\\)", "TVV  * exp(ETA_V + KAPPA_V)", txt)
  stopifnot(sum(grepl("block_kappa", txt)) == 1L)
  mod <- withr::local_tempfile(fileext = ".ferx")
  writeLines(txt, mod)
  # warfarin_iov's [fit_options] says covariance = false.
  fit <- suppressWarnings(ferx_fit(mod, ex$data, method = "foce", verbose = FALSE,
                                   covariance = TRUE, settings = list(maxiter = 30L)))
  skip_if(is.null(fit$se_kappa), "covariance step did not run")
  expect_identical(names(fit$se_kappa), c("KAPPA_CL", "KAPPA_V"))
  expect_identical(unname(fit$estimates[c("KAPPA_CL", "KAPPA_V"), "se"]),
                   unname(fit$se_kappa))
  core <- expect_core_rows(fit, "block")
  expect_length(core, 2L)
  for (i in 1:2) {
    expect_match(core[i], sprintf("SE = %.6f$", fit$se_kappa[[i]]))
  }
})

test_that("ascii spells only the arrow and the kappa glyph differently", {
  skip_on_cran()
  fit <- mbma_variant("weighted")
  utf8 <- kappa_rows(fit, ascii = FALSE)
  ascii <- kappa_rows(fit, ascii = TRUE)
  expect_false(identical(ascii, utf8))
  expect_identical(ascii, gsub("κ", "kappa", gsub("→", "->", utf8, fixed = TRUE),
                               fixed = TRUE))
  expect_true(all(!grepl("[^ -~]", ascii)))
})

test_that("the glue refuses what it cannot read rather than guessing", {
  glue <- getFromNamespace("ferx_rust_kappa_rows", "ferx")
  call <- function(...) {
    args <- utils::modifyList(list(
      omega_iov = 0.04, n_kappa = 1L, kappa_names = "K", kappa_fixed = NULL,
      se_kappa = NULL, kappa_param_types = "additive", kappa_weights = NULL,
      kappa_weight_typical = NULL, covariance_status = "computed", ascii = TRUE
    ), list(...))
    do.call(glue, args)
  }
  expect_identical(call(), "  K                    = 0.040000  (SD = 0.2000)  SE = N/A\n")
  expect_error(call(covariance_status = "done"), "covariance_status")
  expect_error(call(kappa_param_types = "lognormal"), "random-effect scale")
  expect_error(call(omega_iov = c(0.04, 0)), "1 x 1")
  # A loaded .fitrx spells the status in CamelCase; same rows.
  expect_identical(call(covariance_status = "SirFallback"),
                   call(covariance_status = "sir_fallback"))
})
