# fit$held_at_bound (ferx-core #2797): every coordinate the covariance step
# reports no standard error for. A real guard collapse is hard to provoke on a
# bundled example, so these tests put one in the fit's bytes with
# fitrx_edit_fit() and let the engine's loader and verdict read it back.

quiet <- function(x) suppressWarnings(suppressMessages(x))

# `fit` with ETA `eta` of warfarin_fit_cov() listed as collapsed at the guard
# (`reason`, `kind` as the wire spells them), held by the step when `held`.
with_collapse <- function(fit, eta, held, reason = "collapsed_at_guard",
                          kind = "Omega") {
  i <- match(eta, fit$eta_names)
  packed <- match(eta, rownames(fit$cov_matrix))
  stopifnot(!is.na(i), !is.na(packed))
  fitrx_edit_fit(fit, function(w) {
    w$held_at_bound <- list(list(
      name = eta, kind = kind, index = i - 1L, packed_index = packed - 1L,
      reason = reason
    ))
    if (held) {
      mask <- as.logical(unlist(w$cov_fixed))
      if (length(mask) == 0L) mask <- rep(FALSE, nrow(fit$cov_matrix))
      mask[packed] <- TRUE
      w$cov_fixed <- as.list(mask)
    }
    w
  })
}

se_cell <- function(lines, label) {
  row <- grep(paste0("^  ", label, " "), lines, value = TRUE)
  # Not a block's covariance row (`ETA_V ~ ETA_CL : cov = ...`).
  row <- row[!grepl(" ~ ", row, fixed = TRUE)]
  # Nor a [mixture] override of it (`ETA_CL (class 2)`), unless asked for.
  if (!grepl("(class ", label, fixed = TRUE)) {
    row <- row[!grepl(paste0("^  ", label, " \\(class "), row)]
  }
  expect_length(row, 1L)
  sub(".*SE = ([^ ]+).*", "\\1", row)
}

test_that("a real guard collapse prints N/A and is held quietly (ferx-core #2797)", {
  skip_on_cran()
  # ferx-core's tests/guard_hold.rs "stationary" case: ETA_KA sits on the floor
  # of its log packing at the start, `maxiter = 0` keeps it there, and the
  # covariance step holds it.
  # Not `fixtures/guard_hold`: R CMD build drops every directory whose name
  # ends in "old" (`([Oo]ld|\.Rcheck)$`), so the test failed only on CI.
  dir <- test_path("fixtures", "guard_collapse")
  fit <- expect_no_warning(ferx_fit(file.path(dir, "guard_collapse.ferx"),
                                    file.path(dir, "guard_collapse.csv"),
                                    verbose = FALSE))
  expect_identical(
    fit$held_at_bound,
    data.frame(parameter = "ETA_KA", kind = "omega", index = 3L,
               reason = "collapsed_at_guard", class = NA_integer_,
               stringsAsFactors = FALSE)
  )
  # Held by the step: a FIX-like row, so the engine's 0 stays in se_omega,
  # cov_fixed marks it, and cor_matrix does not warn about a zero variance
  # (expect_no_warning above).
  expect_identical(unname(fit$se_omega[3]), 0)
  expect_true(fit$cov_fixed[["ETA_KA"]])
  expect_false(any(startsWith(fit$status_reasons, "collapsed_at_guard")))
  # No surface reports that 0 as a standard error.
  out <- utils::capture.output(print(fit))
  expect_identical(se_cell(out, "ETA_KA"), "N/A")
  expect_true(is.na(fit$estimates$se[fit$estimates$param == "ETA_KA"]))
  rl <- utils::capture.output(ferx_runlog(fit))
  row <- strsplit(trimws(grep("^  ETA_KA ", rl, value = TRUE)), "[[:space:]]+")[[1]]
  expect_identical(tail(row, 2L), c("N/A", "N/A"))
})

test_that("a fit that holds nothing has a zero-row held_at_bound", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  h <- fit$held_at_bound
  expect_s3_class(h, "data.frame")
  expect_identical(names(h), c("parameter", "kind", "index", "reason", "class"))
  expect_identical(nrow(h), 0L)
})

test_that("a collapse the step held prints N/A and leaves the status clean", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  fit <- with_collapse(base, eta, held = TRUE)
  expect_identical(
    fit$held_at_bound,
    data.frame(parameter = eta, kind = "omega", index = 2L,
               reason = "collapsed_at_guard", class = NA_integer_,
               stringsAsFactors = FALSE)
  )
  # Held by the step: core lists no reason for it.
  expect_identical(fit$status, base$status)
  expect_false(any(startsWith(fit$status_reasons, "collapsed_at_guard")))
  # The step's mask reaches cov_fixed, so cor_matrix reads the row as held.
  expect_true(fit$cov_fixed[[eta]])
  # Every surface shows no SE for it; the other ETA keeps its own.
  est <- fit$estimates
  expect_true(all(is.na(est[est$param == eta, c("se", "rse_pct", "lower_95", "upper_95")])))
  expect_true(is.finite(est$se[est$param == base$eta_names[1]]))
  out <- utils::capture.output(quiet(print(fit)))
  expect_identical(se_cell(out, eta), "N/A")
  expect_false(identical(se_cell(out, base$eta_names[1]), "N/A"))
  rl <- utils::capture.output(quiet(ferx_runlog(fit)))
  row <- strsplit(trimws(grep(paste0("^  ", eta, " "), rl, value = TRUE)), "[[:space:]]+")[[1]]
  expect_identical(tail(row, 2L), c("N/A", "N/A"))
  # The fit's se_omega keeps the engine's number, as for a FIX parameter.
  expect_identical(fit$se_omega, base$se_omega)
})

test_that("a collapse the step did not hold is a status reason", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  fit <- with_collapse(base, eta, held = FALSE)
  expect_true(fit$converged)
  expect_identical(fit$status, "converged_with_issues")
  expect_identical(fit$status_label, "CONVERGED WITH ISSUES")
  expect_identical(tail(fit$status_reasons, 1L), paste0("collapsed_at_guard: ", eta))
  line <- grep("STATUS:", utils::capture.output(quiet(print(fit))), value = TRUE)
  expect_match(line, paste0("CONVERGED WITH ISSUES   (", paste(fit$status_reasons, collapse = ", "), ")"),
               fixed = TRUE)
  # Not held, still no SE.
  expect_true(is.na(fit$estimates$se[fit$estimates$param == eta]))
})

test_that("held_at_bound survives ferx_save_fit() / ferx_load_fit()", {
  skip_on_cran()
  fit <- with_collapse(warfarin_fit_cov(), warfarin_fit_cov()$eta_names[2], held = FALSE)
  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$held_at_bound, fit$held_at_bound)
  expect_identical(back$status_reasons, fit$status_reasons)
  expect_identical(back$estimates, fit$estimates)
})

test_that("a held reason or family a later ferx added still loads (ferx-core #2817)", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  later <- with_collapse(base, eta, held = FALSE, reason = "some_future_reason")
  expect_identical(later$held_at_bound$reason, "unknown")
  # Still held: no SE, and no guard-collapse reason it is not.
  expect_true(is.na(later$estimates$se[later$estimates$param == eta]))
  expect_false(any(startsWith(later$status_reasons, "collapsed_at_guard")))
  # A family this build does not know is dropped, not a load failure.
  alien <- with_collapse(base, eta, held = FALSE, kind = "Hyperprior")
  expect_identical(nrow(alien$held_at_bound), 0L)
})

# `fit` with parameter `name` (family `kind`, 1-based `index` within it) listed
# in its bytes as collapsed at the guard and not held by the step.
with_held <- function(fit, name, kind, index) {
  packed <- match(name, rownames(fit$cov_matrix))
  stopifnot(!is.na(packed))
  fitrx_edit_fit(fit, function(w) {
    w$held_at_bound <- list(list(
      name = name, kind = kind, index = index - 1L, packed_index = packed - 1L,
      reason = "collapsed_at_guard"
    ))
    w
  })
}

# The SE column of `ferx_runlog()`'s row for `name`: its last two cells.
runlog_se <- function(fit, name) {
  rl <- utils::capture.output(quiet(ferx_runlog(fit)))
  row <- rl[startsWith(rl, paste0("  ", name, " "))]
  # Not a [mixture] override of it (`ETA_CL (class 2)`).
  row <- row[!startsWith(row, paste0("  ", name, " (class "))]
  expect_length(row, 1L)
  tail(strsplit(trimws(row), "[[:space:]]+")[[1]], 2L)
}

block_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin_block_omega")
      fit <<- quiet(ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                             covariance = TRUE, settings = list(maxiter = 30L)))
    }
    fit
  }
})

mbma_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("mbma_placebo")
      fit <<- quiet(ferx_fit(ex$model, ex$data, verbose = FALSE))
    }
    fit
  }
})

test_that("a held block-omega variance masks its own se_omega slot (#737 review, finding 4)", {
  skip_on_cran()
  base <- block_fit()
  skip_if(is.null(base$cov_matrix), "no covariance on the block fit")
  # The 2x2 block plus ETA_KA: se_omega is the lower triangle, so ETA_V's
  # variance is slot 3, not 2 (slot 2 is the ETA_V ~ ETA_CL covariance).
  expect_length(base$se_omega, 6L)
  fit <- with_held(base, "ETA_V", "Omega", 2L)
  out <- utils::capture.output(quiet(print(fit)))
  expect_identical(se_cell(out, "ETA_V"), "N/A")
  expect_false(identical(se_cell(out, "ETA_CL"), "N/A"))
  expect_false(identical(se_cell(out, "ETA_KA"), "N/A"))
  est <- fit$estimates
  expect_true(is.na(est$se[est$param == "ETA_V"]))
  expect_true(is.finite(est$se[est$param == "ETA_CL"]))
  expect_identical(runlog_se(fit, "ETA_V"), c("N/A", "N/A"))
})

test_that("a held sigma prints N/A everywhere (#737 review, finding 4)", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  sg <- base$sigma_names[1]
  row_of <- function(f) {
    out <- utils::capture.output(quiet(print(f)))
    out[startsWith(out, paste0("  ", sg, " "))]
  }
  fit <- with_held(base, sg, "Sigma", 1L)
  expect_length(row_of(fit), 1L)
  expect_match(row_of(fit), "SE = N/A", fixed = TRUE)
  expect_false(grepl("SE = N/A", row_of(base), fixed = TRUE))
  expect_true(is.na(fit$estimates$se[fit$estimates$param == sg]))
  expect_identical(runlog_se(fit, sg), c("N/A", "N/A"))
})

test_that("a held kappa prints N/A in the KAPPA rows and estimates (#737 review, finding 4)", {
  skip_on_cran()
  base <- mbma_fit()
  skip_if(is.null(base$cov_matrix), "no covariance on the MBMA fit")
  kappa_row <- function(f) {
    rows <- ferx:::.ferx_kappa_rows(f)
    rows[grepl("KAPPA_ARM", rows, fixed = TRUE)][1]
  }
  fit <- with_held(base, "KAPPA_ARM", "Kappa", 1L)
  expect_match(kappa_row(fit), "SE = N/A", fixed = TRUE)
  expect_false(grepl("SE = N/A", kappa_row(base), fixed = TRUE))
  expect_true(is.na(fit$estimates$se[fit$estimates$param == "KAPPA_ARM"]))
})

test_that("a guard-collapsed level theta has no se in theta_levels either (#737 review, finding 2)", {
  skip_on_cran()
  base <- mbma_fit()
  skip_if(is.null(base$cov_matrix), "no covariance on the MBMA fit")
  lv <- "PLACEBO[STUDY=1,TIME=0]"
  i <- match(lv, names(base$theta))
  skip_if(is.na(i), "no such level")
  at <- function(f) which(f$theta_levels$theta_name == lv)
  expect_true(is.finite(base$theta_levels$se[at(base)]))
  fit <- with_held(base, lv, "Theta", i)
  expect_true(is.na(fit$theta_levels$se[at(fit)]))
  expect_true(is.na(fit$estimates$se[fit$estimates$param == lv]))
  # The other levels keep theirs.
  others <- fit$theta_levels$theta_name != lv & !is.na(fit$theta_levels$theta_name)
  expect_identical(fit$theta_levels$se[others], base$theta_levels$se[others])
})

# The real guard-collapse fit (ETA_KA held), with TVKA also held through `how`:
# its covariance row and column zeroed, as the step leaves a held coordinate,
# and listed in held_at_bound ("collapsed") or profiled_levels_at_bound
# ("profiled").
collapse_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      dir <- test_path("fixtures", "guard_collapse")
      fit <<- quiet(ferx_fit(file.path(dir, "guard_collapse.ferx"),
                             file.path(dir, "guard_collapse.csv"), verbose = FALSE))
    }
    fit
  }
})

with_tvka_held <- function(fit, how) {
  k <- match("TVKA", rownames(fit$cov_matrix))
  stopifnot(!is.na(k))
  fitrx_edit_fit(fit, function(w) {
    n <- w$covariance_matrix$rows
    cells <- expand.grid(i = seq_len(n), j = seq_len(n))
    w$covariance_matrix$data[cells$i == k | cells$j == k] <- list(0)
    if (how == "collapsed") {
      w$held_at_bound <- c(w$held_at_bound, list(list(
        name = "TVKA", kind = "Theta", index = k - 1L, packed_index = k - 1L,
        reason = "collapsed_at_guard"
      )))
    } else if (how == "profiled") {
      w$profiled_levels_at_bound <- list(k - 1L)
    }
    w
  })
}

# TVKA as each simulated subject was drawn with: KA = TVKA * exp(ETA_KA).
drawn_tvka <- function(fit) {
  dir <- test_path("fixtures", "guard_collapse")
  sims <- quiet(ferx_simulate_with_uncertainty(
    file.path(dir, "guard_collapse.ferx"), file.path(dir, "guard_collapse.csv"),
    fit, n_uncertainty_draws = 20L, seed = 3L
  ))
  p <- attr(sims, "parameters")
  p$KA / exp(p$ETA_KA)
}

test_that("uncertainty draws keep a held coordinate on its value (#737 review, finding 3)", {
  skip_on_cran()
  base <- collapse_fit()
  tvka <- unname(base$theta[["TVKA"]])
  # Control: TVKA has a standard error, so its draws spread. Only what the fit
  # holds is held (ETA_KA's hold landing on the wrong coordinate fails here).
  free <- drawn_tvka(base)
  expect_gt(max(abs(free / tvka - 1)), 1e-3)
  # Held, it is drawn at its estimate every time, whichever list holds it. A
  # zero covariance row alone lets it scatter about 1e-4 around it, since the
  # regularised Cholesky adds 1e-8 to the diagonal.
  for (how in c("collapsed", "profiled")) {
    held <- drawn_tvka(with_tvka_held(base, how))
    expect_lt(max(abs(held / tvka - 1)), 1e-12, label = how)
  }
  # The zero row alone, listed nowhere: the scatter the lists prevent. ETA_KA,
  # omega index 3, is listed and must be held as itself, not as TVKA (theta
  # index 3), whose row is now zero too.
  zero_only <- drawn_tvka(with_tvka_held(base, "zero_row_only"))
  expect_gt(max(abs(zero_only / tvka - 1)), 1e-6)
})

# ferx-core's tests/guard_hold.rs `mixture_source` case (#2809): two classes,
# class 2 overrides ETA_CL with a variance on the internal guard, and
# `maxiter = 0` keeps it there for the covariance step to hold.
mix_override_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      dir <- test_path("fixtures", "mix_override_collapse")
      fit <<- quiet(ferx_fit(file.path(dir, "mix_override_collapse.ferx"),
                             file.path(dir, "mix_override_collapse.csv"),
                             verbose = FALSE))
    }
    fit
  }
})

test_that("a held mixture override is its own row, not its base's (ferx-core #2809)", {
  skip_on_cran()
  fit <- mix_override_fit()
  expect_false(is.null(fit$cov_matrix))
  expect_identical(
    fit$held_at_bound,
    data.frame(parameter = "ETA_CL (class 2)", kind = "omega", index = 1L,
               reason = "collapsed_at_guard", class = 2L,
               stringsAsFactors = FALSE)
  )
  # One name for the parameter across the fit: the label cov_matrix and
  # cov_fixed give the override, and the one cov_fixed marks as held.
  expect_true(fit$cov_fixed[["ETA_CL (class 2)"]])
  expect_identical(names(fit$cov_fixed)[unlist(fit$cov_fixed)], "ETA_CL (class 2)")
  # The row names ETA_CL by kind + index, but ETA_CL itself is not held: it
  # keeps a standard error on every surface.
  se <- fit$estimates$se[fit$estimates$param == "ETA_CL"]
  expect_length(se, 1L)
  expect_true(is.finite(se) && se > 0)
  out <- utils::capture.output(quiet(print(fit)))
  expect_false(identical(se_cell(out, "ETA_CL"), "N/A"))
  expect_false(any(runlog_se(fit, "ETA_CL") == "N/A"))
  # The override has its own row (ferx-core #2833), free but with no standard
  # error: the engine's 0 for it is not one.
  expect_identical(fit$se_mixture_overrides[["ETA_CL (class 2)"]], 0)
  row <- fit$estimates["ETA_CL (class 2)", ]
  expect_identical(row$class, 2L)
  expect_false(row$fixed)
  expect_true(is.na(row$se))
  expect_identical(se_cell(out, "ETA_CL \\(class 2\\)"), "N/A")
  expect_identical(runlog_se(fit, "ETA_CL (class 2)"), c("N/A", "N/A"))
})

test_that("a held mixture override keeps its class through save/load (ferx-core #2809)", {
  skip_on_cran()
  fit <- mix_override_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$held_at_bound, fit$held_at_bound)
  expect_identical(back$held_at_bound$class, 2L)
  expect_identical(back$estimates, fit$estimates)
})

test_that("uncertainty draws hold a held mixture override, found by its class (ferx-core #2809)", {
  skip_on_cran()
  fit <- mix_override_fit()
  dir <- test_path("fixtures", "mix_override_collapse")
  draws <- function(f) {
    quiet(ferx_simulate_with_uncertainty(
      file.path(dir, "mix_override_collapse.ferx"),
      file.path(dir, "mix_override_collapse.csv"),
      f, n_uncertainty_draws = 20L, seed = 3L
    ))
  }
  held <- draws(fit)
  expect_gt(nrow(held), 0L)
  # The row reaches the engine: without it the override is drawn, not held.
  free <- fit
  free$held_at_bound <- free$held_at_bound[0, ]
  expect_false(identical(draws(free), held))
  # It is found by kind, index and class: the label alone names no coordinate,
  # and a class the model does not override is refused rather than dropped.
  for (cls in c(NA_integer_, 1L)) {
    wrong <- fit
    wrong$held_at_bound$class <- cls
    expect_error(draws(wrong), "names no single coordinate", label = format(cls))
  }
})

test_that("a fit object from before the class column reads every row as a base", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  fit <- quiet(with_collapse(base, eta, held = TRUE))
  fit$held_at_bound$class <- NULL
  expect_true(ferx:::.ferx_se_held(fit, "omega", 2L))
  expect_false(ferx:::.ferx_se_held(fit, "omega", 1L))
})
