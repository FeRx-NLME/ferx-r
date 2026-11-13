# Fits are read-only (#713): core fields refuse assignment, new names are the
# user's own and survive a save, and a save refuses a fit edited around the
# guards.

test_that("each assignment form on a core field errors", {
  fit <- warfarin_fit()
  expect_error(fit$theta <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit[["theta"]] <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit["theta"] <- list(1), class = "ferx_fit_read_only_error")
  # Deleting is changing.
  expect_error(fit$ofv <- NULL, class = "ferx_fit_read_only_error")
  # By position, by a nested path, and by every field at once.
  i <- match("bundle", names(unclass(fit)))
  expect_error(fit[[i]] <- raw(0), class = "ferx_fit_read_only_error")
  expect_error(fit[i] <- list(raw(0)), class = "ferx_fit_read_only_error")
  expect_error(fit[[c("theta", "TVCL")]] <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit[c("my_note", "omega")] <- list(1, 2),
               class = "ferx_fit_read_only_error")
  expect_error(fit$sdtab$DV <- 0, class = "ferx_fit_read_only_error")
  # A factor subscript replaces by its integer code, as a position.
  f <- factor("theta", levels = names(unclass(fit)))
  expect_error(fit[[f]] <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit[f] <- list(1), class = "ferx_fit_read_only_error")
  # A factor whose label is a free name but whose code is a core field's
  # position: the position is what is replaced.
  pos <- match("ofv", names(unclass(fit)))
  g <- factor("my_note", levels = c(paste0("l", seq_len(pos - 1L)), "my_note"))
  expect_identical(as.integer(g), pos)
  expect_error(fit[[g]] <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit[g] <- list(1), class = "ferx_fit_read_only_error")
  # Classed numbers (Date, POSIXct, difftime) replace by position too.
  d <- as.Date(pos, origin = "1970-01-01")
  expect_identical(as.numeric(d), as.numeric(pos))
  expect_error(fit[[d]] <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit[d] <- list(1), class = "ferx_fit_read_only_error")
  expect_error(fit[[as.difftime(pos, units = "secs")]] <- 1,
               class = "ferx_fit_read_only_error")
  expect_error(fit[as.POSIXct(pos, origin = "1970-01-01", tz = "UTC")] <- list(1),
               class = "ferx_fit_read_only_error")
  expect_error(fit[[TRUE]] <- 1, class = "ferx_fit_read_only_error")
  expect_error(fit[-1L] <- list(NULL), class = "ferx_fit_read_only_error")
  err <- tryCatch(fit$theta <- 1, ferx_fit_read_only_error = function(e) e)
  expect_identical(err$fields, "theta")
  expect_match(conditionMessage(err), "read-only")
  # Nothing changed.
  expect_identical(fit, warfarin_fit())
})

test_that("each assignment form works for a new name", {
  fit <- warfarin_fit()
  fit$note1 <- "a"
  fit[["note2"]] <- 1:3
  fit["note3"] <- list(list(x = 1))
  expect_s3_class(fit, "ferx_fit")
  expect_identical(fit$note1, "a")
  expect_identical(fit$note2, 1:3)
  expect_identical(fit$note3, list(x = 1))
  # A user addition is the user's: it can be changed and removed.
  fit$note1 <- "b"
  fit$note2 <- NULL
  expect_identical(fit$note1, "b")
  expect_false("note2" %in% names(fit))
  expect_setequal(ferx:::.ferx_fit_user_fields(fit), c("note1", "note3"))
  # `[[TRUE]] <-` replaces the first element only: allowed when that is the
  # user's own, whatever core fields follow it.
  first <- ferx:::.ferx_fit_seal(c(list(mine = 0), unclass(warfarin_fit())),
                                 user = "mine")
  first[[TRUE]] <- 1
  expect_identical(first$mine, 1)
  # Putting back the value a field already has changes nothing.
  same <- warfarin_fit()
  same$theta <- same$theta
  expect_identical(same, warfarin_fit())
  # Where the fit's files are now is not a core field.
  fit$model_path <- "elsewhere.ferx"
  expect_identical(fit$model_path, "elsewhere.ferx")
})

test_that("unclass(), edit, reclass, save is refused", {
  fit <- unclass(warfarin_fit())
  fit$theta[[1L]] <- fit$theta[[1L]] * 2
  class(fit) <- "ferx_fit"
  out <- tempfile(fileext = ".fitrx")
  err <- tryCatch(ferx_save_fit(fit, out), ferx_fit_edited_error = function(e) e)
  expect_s3_class(err, "ferx_fit_edited_error")
  expect_identical(err$fields, "theta")
  expect_false(file.exists(out))

  # A removed field, and a stripped digest, too.
  gone <- unclass(warfarin_fit())
  gone$sdtab <- NULL
  class(gone) <- "ferx_fit"
  expect_error(ferx_save_fit(gone, out), class = "ferx_fit_edited_error")
  bare <- warfarin_fit()
  attr(bare, "ferx_core_digest") <- NULL
  expect_error(ferx_save_fit(bare, out), class = "ferx_fit_edited_error")
  expect_false(file.exists(out))
})

test_that("a user addition survives a save and load", {
  fit <- warfarin_fit()
  fit$notes <- "run 1, reviewed"
  fit$table <- data.frame(x = c(0.1, 1 / 3), y = c("a", "b"))
  fit$cfg <- list(seed = 42L, tags = c("a", "b"))
  out <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, out)
  back <- ferx_load_fit(out)
  expect_identical(back$notes, fit$notes)
  expect_identical(back$table, fit$table)
  expect_identical(back$cfg, fit$cfg)
  expect_setequal(ferx:::.ferx_fit_user_fields(back), c("notes", "table", "cfg"))
  # Still the user's after a load: changeable, and core fields still are not.
  back$notes <- "edited"
  expect_identical(back$notes, "edited")
  expect_error(back$theta <- 1, class = "ferx_fit_read_only_error")

  # Removed after the load, it does not come back on the next save.
  back$cfg <- NULL
  out2 <- tempfile(fileext = ".fitrx")
  ferx_save_fit(back, out2)
  again <- ferx_load_fit(out2)
  expect_false("cfg" %in% names(again))
  expect_identical(again$notes, "edited")

  # A step that returns a new fit keeps the input fit's additions.
  stepped <- ferx_covariance(back)
  expect_identical(stepped$notes, "edited")
  expect_false("cfg" %in% names(stepped))
})

test_that("an unserializable user addition warns and is left out", {
  fit <- warfarin_fit()
  fit$env <- new.env()
  fit$ok <- 1
  out <- tempfile(fileext = ".fitrx")
  expect_warning(ferx_save_fit(fit, out),
                 class = "ferx_fit_user_field_dropped_warning")
  back <- ferx_load_fit(out)
  expect_false("env" %in% names(back))
  expect_identical(back$ok, 1)

  # Functions, S4 objects, and classes or attributes outside the allowlist are
  # not plain data either.
  for (value in list(function(x) x, methods::getClass("numeric"),
                     structure(1, class = "my_class"), structure(1, note = "x"),
                     as.POSIXlt("2026-01-01", tz = "UTC"), 1i, as.raw(1))) {
    other <- warfarin_fit()
    other$bad <- value
    expect_warning(ferx_save_fit(other, tempfile(fileext = ".fitrx")),
                   class = "ferx_fit_user_field_dropped_warning")
  }
})

test_that("the allowlisted kinds of plain data round-trip exactly", {
  fit <- warfarin_fit()
  kinds <- list(
    dbl = c(a = 1.5, b = NaN, c = Inf, d = -Inf, e = NA, f = 1 / 3),
    int = c(1L, NA, .Machine$integer.max), lgl = c(TRUE, NA), chr = c("x", NA),
    empty = character(), nested = list(a = list(b = NULL, c = "q")),
    df = data.frame(x = c(0.1, 2), f = factor(c("u", "v")), s = c("a", "b")),
    ord = ordered(c("lo", "hi"), levels = c("lo", "hi")),
    mat = matrix(c(0.1, 2, 3, 4), 2, dimnames = list(c("r1", "r2"), c("a", "b"))),
    date = as.Date("2026-11-12"),
    time = as.POSIXct("2026-11-12 10:11:12.5", tz = "UTC"),
    dt = as.difftime(5, units = "mins"), tab = table(c(1, 1, 2))
  )
  for (k in names(kinds)) fit[[k]] <- kinds[[k]]
  out <- tempfile(fileext = ".fitrx")
  expect_no_warning(ferx_save_fit(fit, out))
  back <- ferx_load_fit(out)
  for (k in names(kinds)) expect_identical(back[[k]], kinds[[k]], info = k)
})

test_that("a crafted bundle's user objects are decoded as data, never deserialized", {
  fit <- warfarin_fit()
  # A jsonlite::serializeJSON() payload (which unserializeJSON() would hand to
  # base unserialize() or a namespace) is not a node of the encoding, and nor
  # are a node of a type, class or attribute outside the allowlist.
  payload <- as.character(jsonlite::serializeJSON(quote(stop("ran"))))
  # The helper rebuilds the fit through the same decoder, which warns there too.
  crafted <- suppressWarnings(fitrx_edit_extras(fit, function(r) {
    r$user <- list(
      serialized = payload,
      closure = list(type = "closure", values = list()),
      s4 = list(type = "S4", values = list()),
      classed = list(type = "double", values = list(1),
                     attributes = list(class = list(type = "character",
                                                    values = list("evil")))),
      attr = list(type = "double", values = list(1),
                  attributes = list(srcref = list(type = "integer",
                                                  values = list(1L)))),
      fine = list(type = "character", values = list("kept"))
    )
    r
  }))
  path <- tempfile(fileext = ".fitrx")
  writeBin(crafted$bundle, path)
  w <- NULL
  loaded <- withCallingHandlers(
    ferx_load_fit(path),
    ferx_fit_user_field_dropped_warning = function(cnd) {
      w <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_setequal(w$fields, c("serialized", "closure", "s4", "classed", "attr"))
  expect_identical(loaded$fine, "kept")
  expect_identical(ferx:::.ferx_fit_user_fields(loaded), "fine")
})

test_that("internal code building a fit is not tripped by the guards", {
  fit <- warfarin_fit()
  expect_false(is.null(ferx:::.ferx_fit_core_fields(fit)))
  out <- tempfile(fileext = ".fitrx")
  # Save (with data, so the canonical view is built), load, and the steps that
  # rebuild a fit from new bytes.
  expect_no_error(ferx_save_fit(fit, out, include_data = TRUE))
  loaded <- ferx_load_fit(out)
  expect_no_error(ferx_save_fit(loaded, tempfile(fileext = ".fitrx")))
  cov <- ferx_covariance(loaded)
  expect_no_error(ferx_save_fit(cov, tempfile(fileext = ".fitrx")))
  expect_identical(ferx:::.ferx_fit_user_fields(cov), character())
  # A fit built by hand (no digest) is a plain list: the guards stay out of it.
  stub <- structure(list(theta = 1), class = "ferx_fit")
  stub$theta <- 2
  expect_identical(stub$theta, 2)
})
