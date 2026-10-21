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

test_that("unnamed kappas print the KAPPA<i> fallback on every line that names them", {
  # Review row 2: core's own fallback is a bare `KAPPA` for every row, which
  # AGENTS.md's label table does not allow and the correlations below did not use.
  fit <- make_fake_fit(
    omega = matrix(0.10, 1, 1),
    omega_iov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    se_kappa = c(0.01, 0.005), shrinkage_kappa = c(0.1, 0.2),
    kappa_param_types = c("log_normal", "log_normal")
  )
  rows <- printed_kappa_rows(fit)
  expect_length(rows, 2L)
  expect_match(rows[1L], "^  KAPPA1 +=")
  expect_match(rows[2L], "^  KAPPA2 +=")
  out <- utils::capture.output(print(fit))
  expect_true(any(startsWith(out, "  KAPPA2 ~ KAPPA1 : cov = 0.010000")))
  expect_match(grep("KAPPA1: ", out, value = TRUE), "KAPPA1: 10.0%   KAPPA2: 20.0%", fixed = TRUE)
  # An empty name in an otherwise named fit falls back too, at its own index.
  fit$kappa_names <- c("KAPPA_CL", "")
  rows <- printed_kappa_rows(fit)
  expect_match(rows[1L], "^  KAPPA_CL +=")
  expect_match(rows[2L], "^  KAPPA2 +=")
  out <- utils::capture.output(print(fit))
  expect_true(any(startsWith(out, "  KAPPA2 ~ KAPPA_CL : cov = 0.010000")))
})

test_that("an NA se_kappa entry prints SE = N/A, not NaN (#546)", {
  # A .fitrx whose se_kappa has a null entry reloads as NA; the glue hands it
  # to core as NaN, and core prints a non-finite SE as a missing one
  # (ferx-core #2225).
  fit <- make_fake_fit(
    omega = matrix(0.10, 1, 1),
    omega_iov = diag(c(0.04, 0.09)), kappa_names = c("KA", "KB"),
    se_kappa = c(0.01, NA), kappa_param_types = c("log_normal", "log_normal"),
    covariance_status = "computed"
  )
  rows <- printed_kappa_rows(fit)
  expect_length(rows, 2L)
  expect_match(rows[1L], "SE = 0.010000$")
  expect_match(rows[2L], "^  KB +=")
  expect_match(rows[2L], "SE = N/A$")
  expect_false(any(grepl("NaN", rows, fixed = TRUE)))
})

test_that("print refuses an unreadable kappa scale or covariance status by name (#545)", {
  base <- make_fake_fit(
    omega = matrix(0.10, 1, 1),
    omega_iov = diag(c(0.04, 0.02)), kappa_names = c("KAPPA_CL", "KAPPA_V"),
    kappa_param_types = c("additive", "log_normal"), covariance_status = "computed"
  )
  expect_length(printed_kappa_rows(base), 2L)
  refuses <- function(fit, field, pattern) {
    err <- expect_error(print(fit), class = "ferx_fit_field_error")
    expect_identical(err$field, field)
    expect_match(conditionMessage(err), pattern, fixed = TRUE)
  }
  f <- base; f$kappa_param_types <- c("additive", NA)
  refuses(f, "kappa_param_types", "has NA at position 2")
  f <- base; f$kappa_param_types <- c("lognormal", "additive")
  refuses(f, "kappa_param_types", "has \"lognormal\" at position 1")
  # A short vector would print KAPPA_V log-normal by default (#562).
  f <- base; f$kappa_param_types <- "additive"
  refuses(f, "kappa_param_types",
          "fit$kappa_param_types has 1 entry; expected 2, one per row of fit$omega_iov, or none.")
  f <- base; f$covariance_status <- NA_character_
  refuses(f, "covariance_status", "fit$covariance_status is NA_character_")
  f <- base; f$covariance_status <- "done"
  refuses(f, "covariance_status", "fit$covariance_status is \"done\"")
  f <- base; f$covariance_status <- c("computed", "failed")
  refuses(f, "covariance_status", "fit$covariance_status is c(\"computed\", \"failed\")")
  # Both sides of each gate: an empty scale vector (a fit from before
  # ferx-core #1643), a missing status, and a reloaded fit's CamelCase status
  # are all readable.
  for (ok in list(list(kappa_param_types = character()),
                  list(covariance_status = NULL),
                  list(covariance_status = "SirFallback"),
                  list(covariance_status = "NotRequested"))) {
    f <- utils::modifyList(base, ok)
    class(f) <- class(base)
    expect_length(printed_kappa_rows(f), 2L)
  }
})

test_that("print's field pre-check and the glue agree on every scale and status (#551)", {
  # The pre-check takes its lists from the glue, so every spelling the glue
  # lists passes it and reaches the formatter, and an unknown one is refused
  # naming the glue's list.
  types <- ferx:::ferx_rust_known_eta_param_types()
  statuses <- ferx:::ferx_rust_known_covariance_statuses()
  expect_true(all(c("log_normal", "additive", "logit", "logit_probability", "custom") %in% types))
  expect_true(all(c("computed", "failed", "not_requested", "sir_fallback") %in% statuses))
  for (t in types) {
    for (s in statuses) {
      fit <- make_fake_fit(
        omega = matrix(0.10, 1, 1),
        omega_iov = matrix(0.04, 1, 1), kappa_names = "KAPPA_CL",
        kappa_param_types = t, covariance_status = s
      )
      expect_length(printed_kappa_rows(fit), 1L)
    }
  }
  check <- getFromNamespace(".ferx_check_row_fields", "ferx")
  err <- expect_error(check("lognormal", "kappa_param_types", 1L, "omega_iov", "computed"),
                      class = "ferx_fit_field_error")
  expect_match(conditionMessage(err), paste(types, collapse = ", "), fixed = TRUE)
  err <- expect_error(check(types[1L], "kappa_param_types", 1L, "omega_iov", "done"),
                      class = "ferx_fit_field_error")
  expect_match(conditionMessage(err), paste(statuses, collapse = ", "), fixed = TRUE)
})

test_that("kappa shrinkage alone opens SHRINKAGE, and over 30% is flagged", {
  # Review row 5: the kappa term of has_shrinkage and the [!] flag.
  fit <- make_fake_fit(
    omega = matrix(0.10, 1, 1),
    omega_iov = diag(c(0.04, 0.02)), kappa_names = c("KAPPA_CL", "KAPPA_V"),
    shrinkage_kappa = c(0.4, 0.1)
  )
  expect_null(fit$shrinkage_eta)
  expect_null(fit$shrinkage_eps)
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("SHRINKAGE", out, fixed = TRUE)))
  line <- grep("KAPPA_CL: ", out, value = TRUE)
  expect_length(line, 1L)
  expect_match(line, "KAPPA_CL: 40.0% [!]   KAPPA_V: 10.0%", fixed = TRUE)
  expect_false(grepl("KAPPA_V: 10.0% [!]", line, fixed = TRUE))
  # No kappa shrinkage either: no SHRINKAGE section at all.
  fit$shrinkage_kappa <- NULL
  expect_false(any(grepl("SHRINKAGE", utils::capture.output(print(fit)), fixed = TRUE)))
})

test_that("the glue refuses what it cannot read rather than guessing", {
  glue <- getFromNamespace("ferx_rust_kappa_rows", "ferx")
  call <- function(...) {
    args <- utils::modifyList(list(
      omega_iov = 0.04, n_kappa = 1L, kappa_names = "K", kappa_fixed = NULL,
      se_kappa = NULL, kappa_param_types = "additive", kappa_weights = NULL,
      kappa_weight_typical = NULL, kappa_levels = character(),
      covariance_status = "computed", ascii = TRUE
    ), list(...))
    do.call(glue, args)
  }
  expect_identical(call(), "  K                    = 0.040000  (SD = 0.2000)  SE = N/A\n")
  # Nested levels (ferx-core #1082): the kappa's level follows its name.
  expect_match(call(kappa_levels = "STUDY"), "^  K at STUDY +=")
  expect_error(call(covariance_status = "done"), "covariance_status")
  expect_error(call(kappa_param_types = "lognormal"), "random-effect scale")
  expect_error(call(omega_iov = c(0.04, 0)), "1 x 1")
  # One scale per kappa or none (#562): a short vector would read the missing
  # kappa as log-normal, a long one is misaligned. Empty is a pre-#1643 fit.
  two <- list(omega_iov = c(0.04, 0, 0, 0.09), n_kappa = 2L, kappa_names = c("K1", "K2"))
  expect_error(do.call(call, c(two, list(kappa_param_types = "additive"))),
               "kappa_param_types has length 1, not 2 (one per row) or 0", fixed = TRUE)
  expect_error(call(kappa_param_types = c("additive", "additive")),
               "kappa_param_types has length 2, not 1 (one per row) or 0", fixed = TRUE)
  expect_match(call(kappa_param_types = character()), "(CV% = ", fixed = TRUE)
  # A loaded .fitrx spells the status in CamelCase; same rows.
  expect_identical(call(covariance_status = "SirFallback"),
                   call(covariance_status = "sir_fallback"))
})
