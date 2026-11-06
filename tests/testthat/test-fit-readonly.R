# Read-only core fit fields (#713). The fixture is a small sealed list, so the
# assignment rules need no engine; what needs a real fit uses `warfarin_fit()`.

sealed <- function() {
  ferx:::.ferx_fit_seal(list(theta = c(TVCL = 0.1), omega = diag(2), ofv = 100))
}

test_that("each assignment form on a core field is an error", {
  fit <- sealed()
  expect_error(fit$theta <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit[["theta"]] <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit["theta"] <- list(1), class = "ferx_fit_readonly_error")
  expect_error(fit$theta <- NULL, class = "ferx_fit_readonly_error")
  expect_error(fit$theta[1] <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit$theta["TVCL"] <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit[["theta"]][1] <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit[[1]] <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit[1] <- list(1), class = "ferx_fit_readonly_error")
  expect_error(fit[-1] <- list(1, 2), class = "ferx_fit_readonly_error")
  expect_error(fit[c(TRUE, FALSE, FALSE)] <- list(1),
               class = "ferx_fit_readonly_error")
  expect_error(fit[] <- list(1, 2, 3), class = "ferx_fit_readonly_error")
  expect_error(fit[c("mine", "theta")] <- list(1, 2),
               class = "ferx_fit_readonly_error")
  expect_error(fit[[c("theta", "TVCL")]] <- 1, class = "ferx_fit_readonly_error")
  expect_identical(fit, sealed())
})

test_that("the error says fits are read-only and names the field", {
  fit <- sealed()
  err <- tryCatch(fit$theta <- 1, ferx_fit_readonly_error = identity)
  expect_match(conditionMessage(err), "fit$theta is read-only", fixed = TRUE)
  expect_identical(err$field, "theta")
})

test_that("a core field the fit lacks is still core", {
  fit <- sealed()
  expect_false("se_theta" %in% names(fit))
  expect_error(fit$se_theta <- 0.1, class = "ferx_fit_readonly_error")
  expect_error(fit[["se_theta"]] <- 0.1, class = "ferx_fit_readonly_error")
  expect_error(fit["se_theta"] <- list(0.1), class = "ferx_fit_readonly_error")
})

test_that("a name `$` would read in place of an absent core field is refused", {
  fit <- sealed()
  expect_error(fit$se_theta_note <- "x", "partially",
               class = "ferx_fit_readonly_error")
  expect_null(fit$se_theta)
})

test_that("each assignment form works for a new name", {
  fit <- sealed()
  fit$a <- 1
  fit[["b"]] <- 2
  fit["c"] <- list(3)
  fit[c("d", "e")] <- list(4, 5)
  expect_identical(unlist(fit[c("a", "b", "c", "d", "e")]),
                   c(a = 1, b = 2, c = 3, d = 4, e = 5))
  expect_s3_class(fit, "ferx_fit")
  expect_true(ferx:::.ferx_fit_intact(fit))
})

test_that("an addition is the user's to change and remove", {
  fit <- sealed()
  fit$note <- "a"
  fit$note <- "b"
  fit$note[2] <- "c"
  fit[["note"]][3] <- "d"
  expect_identical(fit$note, c("b", "c", "d"))
  fit[[match("note", names(fit))]] <- "e"
  expect_identical(fit$note, "e")
  fit$note <- NULL
  expect_false("note" %in% names(fit))
  expect_true(ferx:::.ferx_fit_intact(fit))
})

test_that("a new field needs a name", {
  fit <- sealed()
  expect_error(fit[[length(fit) + 1L]] <- 1, "needs a name",
               class = "ferx_fit_readonly_error")
  expect_error(fit[length(fit) + 1L] <- list(1), "needs a name",
               class = "ferx_fit_readonly_error")
  expect_error(fit[[""]] <- 1, "needs a name", class = "ferx_fit_readonly_error")
})

test_that("a fit with no seal still guards the fields it holds", {
  # An `.rds` from before seals, or a seal stripped or overwritten by hand.
  for (fit in list(structure(list(theta = 1), class = "ferx_fit"),
                   structure(list(theta = 1), class = "ferx_fit",
                             ferx_seal = "nonsense"))) {
    expect_error(fit$theta <- 2, class = "ferx_fit_readonly_error")
    expect_error(fit[["theta"]] <- 2, class = "ferx_fit_readonly_error")
    expect_error(fit["theta"] <- list(2), class = "ferx_fit_readonly_error")
    expect_error(fit$se_theta <- 2, class = "ferx_fit_readonly_error")
    fit$note <- "new"
    expect_identical(fit$note, "new")
  }
})

test_that("the engine's wire names cannot be added", {
  fit <- sealed()
  for (name in ferx:::.FERX_FIT_WIRE_NAMES) {
    expect_false(name %in% names(fit))
    expect_error(fit[[name]] <- 1, class = "ferx_fit_readonly_error", info = name)
  }
  expect_error(fit$cov_condition_number <- 1, class = "ferx_fit_readonly_error")
  expect_error(fit["cov_eigenvalues"] <- list(1),
               class = "ferx_fit_readonly_error")
})

test_that("every interface field the writer keeps is reserved", {
  fit <- sealed()
  for (name in c(ferx:::.FITRX_INTERFACE_KEYS, "model_structure", "cond_dist")) {
    expect_error(fit[[name]] <- 1, class = "ferx_fit_readonly_error", info = name)
  }
})

test_that("the digest follows the values, not how R holds them", {
  fit <- ferx:::.ferx_fit_seal(list(a = 1:10, b = data.frame(x = 1:3)))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(fit, path)
  expect_true(ferx:::.ferx_fit_intact(readRDS(path)))
})

test_that("every way past the assignment methods breaks the digest", {
  fit <- sealed()
  edited <- unclass(fit)
  edited$theta[[1]] <- 0.2
  class(edited) <- "ferx_fit"
  expect_false(ferx:::.ferx_fit_intact(edited))

  dropped <- unclass(fit)
  dropped$ofv <- NULL
  class(dropped) <- "ferx_fit"
  expect_false(ferx:::.ferx_fit_intact(dropped))

  renamed <- fit
  names(renamed)[[1]] <- "theta2"
  expect_false(ferx:::.ferx_fit_intact(renamed))

  forged <- fit
  attr(forged, "ferx_seal") <- "nonsense"
  expect_false(ferx:::.ferx_fit_intact(forged))
})

test_that("a name the assignment methods refuse breaks the digest when forced in", {
  force_in <- function(name) {
    x <- unclass(sealed())
    x[[name]] <- 1
    class(x) <- "ferx_fit"
    x
  }
  # Absent core field, wire name, and a name `$` would read as an absent field.
  for (name in c("se_theta", "cov_condition_number", "se_theta_note")) {
    expect_false(ferx:::.ferx_fit_intact(force_in(name)), info = name)
  }
  # An ordinary name is an addition however it arrived.
  expect_true(ferx:::.ferx_fit_intact(force_in("note")))
})

test_that("internal code building a fit is not tripped by the guards", {
  skip_on_cran()
  fit <- warfarin_fit()
  expect_s3_class(fit, "ferx_fit")
  expect_true(ferx:::.ferx_fit_intact(fit))
  expect_length(ferx:::.ferx_fit_user_fields(fit), 0L)

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_true(ferx:::.ferx_fit_intact(loaded))
  expect_length(ferx:::.ferx_fit_user_fields(loaded), 0L)
})

test_that("the post-hoc steps return an intact fit and keep additions", {
  skip_on_cran()
  fit <- warfarin_fit()
  fit$note <- "kept"

  cov <- suppressWarnings(ferx_covariance(fit))
  expect_true(ferx:::.ferx_fit_intact(cov))
  expect_identical(cov$note, "kept")
  expect_identical(ferx:::.ferx_fit_user_fields(cov), "note")
  # What the step filled in is core, though the input fit lacked it.
  expect_error(cov$se_theta <- 1, class = "ferx_fit_readonly_error")

  npde <- ferx_calc_npde(cov, nsim = 5L, seed = 1L)
  expect_true(ferx:::.ferx_fit_intact(npde))
  expect_identical(npde$note, "kept")
  expect_true(ferx:::.ferx_fit_intact(cov))

  sir <- suppressWarnings(ferx_sir(cov, sir_samples = 50L, sir_resamples = 10L,
                                   sir_seed = 1L))
  expect_true(ferx:::.ferx_fit_intact(sir))
  expect_identical(sir$note, "kept")
})

test_that("a post-hoc step cannot launder an edit", {
  skip_on_cran()
  edited <- unclass(warfarin_fit())
  edited$theta[[1]] <- edited$theta[[1]] * 2
  class(edited) <- "ferx_fit"
  expect_error(ferx_covariance(edited), class = "ferx_fit_modified_error")
  expect_error(ferx_sir(edited), class = "ferx_fit_modified_error")
  expect_error(ferx_calc_npde(edited), class = "ferx_fit_modified_error")
})

test_that("a post-hoc step does not make a fit with no seal saveable", {
  raw <- unclass(sealed())
  attr(raw, "ferx_seal") <- NULL
  expect_s3_class(ferx:::.ferx_fit_reseal(raw), "ferx_fit")
  expect_null(attr(ferx:::.ferx_fit_reseal(raw), "ferx_seal"))
})

test_that("unclass(), edit, reclass, save is refused", {
  skip_on_cran()
  path <- tempfile(fileext = ".fitrx")
  edited <- unclass(warfarin_fit())
  edited$theta[[1]] <- edited$theta[[1]] * 2
  class(edited) <- "ferx_fit"
  expect_error(ferx_save_fit(edited, path), "core fields",
               class = "ferx_fit_modified_error")
  expect_false(file.exists(path))
})

test_that("nothing added past the methods reaches the bundle", {
  skip_on_cran()
  path <- tempfile(fileext = ".fitrx")
  fit <- warfarin_fit()
  expect_null(fit$se_theta)
  for (name in c("se_theta", "se_theta_note", "cov_condition_number",
                 "cov_eigenvalues")) {
    edited <- unclass(fit)
    edited[[name]] <- rep(0.1, length(fit$theta))
    class(edited) <- "ferx_fit"
    expect_error(ferx_save_fit(edited, path), class = "ferx_fit_modified_error",
                 info = name)
    expect_false(file.exists(path), info = name)
  }
})

test_that("a fit with its digest removed is refused", {
  skip_on_cran()
  path <- tempfile(fileext = ".fitrx")
  expect_error(ferx_save_fit(editable_fit(warfarin_fit()), path), "no digest",
               class = "ferx_fit_modified_error")
  expect_false(file.exists(path))
})

test_that("a user addition survives a round trip, exactly", {
  skip_on_cran()
  fit <- warfarin_fit()
  additions <- list(
    note = "final run",
    number = 1 / 3,
    special = c(a = NA, b = Inf, c = -Inf, d = NaN, e = 1e-300),
    int = 1:3,
    flag = TRUE,
    nothing = NULL,
    frame = data.frame(id = 1:2, grp = factor(c("a", "b")), wt = c(70.1, 80.2)),
    mat = matrix(c(0.1, 0.2, 0.3, 0.4), 2, dimnames = list(NULL, c("p", "q"))),
    nested = list(when = as.Date("2026-11-05"), inner = list(x = 2 / 3))
  )
  fit[names(additions)] <- additions

  path <- tempfile(fileext = ".fitrx")
  resaved <- tempfile(fileext = ".fitrx")
  on.exit(unlink(c(path, resaved)))
  expect_no_warning(ferx_save_fit(fit, path))
  loaded <- ferx_load_fit(path)
  expect_setequal(ferx:::.ferx_fit_user_fields(loaded), names(additions))
  for (key in names(additions)) {
    expect_identical(loaded[[key]], additions[[key]], info = key)
  }
  expect_true("nothing" %in% names(loaded))
  expect_true(ferx:::.ferx_fit_intact(loaded))

  # Still an addition after loading: editable, and saved again.
  loaded$note <- "reviewed"
  ferx_save_fit(loaded, resaved)
  again <- ferx_load_fit(resaved)
  expect_identical(again$note, "reviewed")
  expect_identical(again$frame, additions$frame)
})

test_that("an unserializable addition warns and is left out", {
  skip_on_cran()
  fit <- warfarin_fit()
  fit$keep <- 1
  fit$fun <- function(x) x
  fit$env <- new.env()
  fit$model_lm <- lm(dist ~ speed, cars)
  fit$deep <- list(ok = 1, bad = list(sum))

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  w <- expect_warning(ferx_save_fit(fit, path),
                      class = "ferx_fit_addition_dropped_warning")
  expect_setequal(w$fields, c("fun", "env", "model_lm", "deep"))
  expect_match(conditionMessage(w), "fit$fun", fixed = TRUE)

  loaded <- ferx_load_fit(path)
  expect_identical(loaded$keep, 1)
  expect_false(any(c("fun", "env", "model_lm", "deep") %in% names(loaded)))
})

test_that("the reader unpacks data and nothing else", {
  ok <- as.character(jsonlite::serializeJSON(list(a = 1, b = "x"), digits = I(17)))
  expect_identical(ferx:::.fitrx_addition_from_wire(ok), list(a = 1, b = "x"))

  # `unserializeJSON()` would rebuild each of these; the first two call
  # `unserialize()` and `getNamespace()` on what the file says.
  hostile <- c(
    as.character(jsonlite::serializeJSON(sum)),
    '{"type":"namespace","attributes":{},"value":{"name":"stats"}}',
    as.character(jsonlite::serializeJSON(function(x) x)),
    as.character(jsonlite::serializeJSON(quote(a + b))),
    '{"type":"S4","attributes":{},"value":{"class":"x","package":"y"}}',
    # Reached only through a list, or an attribute.
    sprintf('{"type":"list","attributes":{},"value":[%s]}',
            as.character(jsonlite::serializeJSON(sum))),
    sprintf('{"type":"double","attributes":{"x":%s},"value":[1]}',
            as.character(jsonlite::serializeJSON(sum))),
    '{"type":"double"',
    "[1, 2]"
  )
  for (text in hostile) {
    expect_error(ferx:::.fitrx_addition_from_wire(text), info = text)
  }
})

test_that("an addition that cannot come back costs that addition, with a warning", {
  fit <- sealed()
  good <- as.character(jsonlite::serializeJSON(1, digits = I(17)))
  bad <- as.character(jsonlite::serializeJSON(sum))
  wire <- list(good = good, bad = bad, theta = good)
  w <- expect_warning(out <- ferx:::.fitrx_restore_user_additions(fit, wire),
                      class = "ferx_fit_addition_dropped_warning")
  expect_setequal(w$fields, c("bad", "theta"))
  expect_identical(out$good, 1)
  expect_identical(out$theta, fit$theta)
  expect_false("bad" %in% names(out))
})
