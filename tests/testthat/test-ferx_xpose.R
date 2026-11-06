
# ---- header from test-xpose.R ----
# Tests for the Xpose interop layer. These exercise the pure table/role
# builder (`.ferx_xpose_frame()` and helpers), which has no Xpose dependency.
# Backend object construction is only smoke-tested for its guard behaviour,
# since xpose / xpose4 are optional Suggests and may be absent.

# Minimal fake ferx_fit with covariates, etas, and individual parameters.
fake_fit <- function(with_covtab = TRUE, echo_wt = FALSE) {
  sdtab <- data.frame(
    ID    = c(1, 1, 2, 2),
    TIME  = c(1, 4, 1, 4),
    DV    = c(10, 8, 12, 9),
    PRED  = c(11, 7, 13, 8),
    IPRED = c(10.5, 7.5, 12.5, 8.5),
    CWRES = c(0.1, -0.2, 0.3, -0.1),
    IWRES = c(0.2, -0.1, 0.2, -0.2),
    NPDE  = c(0.3, -0.4, 0.1, -0.2),
    NPD   = c(0.25, -0.35, 0.15, -0.15),
    TAD   = c(1, 4, 1, 4),
    stringsAsFactors = FALSE
  )
  if (echo_wt) sdtab$WT <- c(70, 70, 80, 80)

  # Mirror production typing: sdtab IDs are numeric (Rust f64), but the
  # per-subject / covtab tables carry character IDs.
  fit <- list(
    sdtab = sdtab,
    individual_estimates = data.frame(ID = c("1", "2"), CL = c(1.1, 1.4),
                                      V = c(30, 35), stringsAsFactors = FALSE),
    ebe_etas = data.frame(ID = c("1", "2"), ETA_CL = c(0.05, -0.03),
                          ETA_V = c(-0.1, 0.2), stringsAsFactors = FALSE),
    covariate_types = c(WT = "continuous", SEX = "categorical"),
    ofv = 123.45,
    theta = c(CL = 1.0, V = 32),
    omega = matrix(c(0.04, 0, 0, 0.09), 2, 2)
  )
  if (with_covtab) {
    fit$covtab <- data.frame(
      ID   = c("1", "1", "2", "2"),
      TIME = c(0, 1, 0, 1),
      EVID = c(1, 0, 1, 0),
      WT   = c(70, 70, 80, 80),
      SEX  = c(0, 0, 1, 1),
      stringsAsFactors = FALSE
    )
  }
  class(fit) <- c("ferx_editable_fit", "ferx_fit")
  fit
}

# A widened optimizer trace (issue #168): the fixed scalar columns the Rust
# backend writes, plus one `val:<param>` and one `grad:<param>` column per
# optimized coordinate (declared names, as ferx-core emits them).
fake_widened_trace <- function(n = 6L, method = "gn", grad = TRUE,
                               params = c("TVCL", "TVV", "ETA_CL", "PROP_ERR")) {
  tr <- data.frame(
    iter      = seq_len(n),
    method    = method,
    phase     = "",
    ofv       = seq(100, by = -1, length.out = n),
    wall_ms   = seq(10, by = 10, length.out = n),
    grad_norm = if (grad) seq(1, by = -0.1, length.out = n) else NA_real_,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (i in seq_along(params)) {
    tr[[paste0("val:", params[i])]] <- seq(i, by = 0.1, length.out = n)
  }
  for (i in seq_along(params)) {
    tr[[paste0("grad:", params[i])]] <-
      if (grad) seq(0.5 * i, by = -0.01, length.out = n) else NA_real_
  }
  tr
}

# A fake fit carrying a widened trace (via fit$trace, which ferx_trace() reads).
fake_fit_trace <- function(...) {
  fit <- fake_fit()
  fit$converged <- TRUE
  fit$trace <- fake_widened_trace(...)
  fit
}




















# --- backend integration (only when the packages are installed) --------------






test_that("frame derives RES/IRES and NA WRES", {
  tbl <- .ferx_xpose_frame(fake_fit())
  df <- tbl$data
  expect_equal(df$RES, df$DV - df$PRED)
  expect_equal(df$IRES, df$DV - df$IPRED)
  expect_true(all(is.na(df$WRES)))
})
test_that("etas and individual parameters are joined per subject (repeated by row)", {
  tbl <- .ferx_xpose_frame(fake_fit())
  df <- tbl$data
  expect_true(all(c("CL", "V", "ETA_CL", "ETA_V") %in% names(df)))
  # subject 1 rows get subject 1's values
  expect_equal(df$CL[df$ID == 1], c(1.1, 1.1))
  expect_equal(df$ETA_V[df$ID == 2], c(0.2, 0.2))
  expect_setequal(tbl$params, c("CL", "V"))
  expect_setequal(tbl$etas, c("ETA_CL", "ETA_V"))
})
test_that("covariates split into continuous / categorical via covariate_types", {
  tbl <- .ferx_xpose_frame(fake_fit())
  expect_setequal(tbl$cont, "WT")
  expect_setequal(tbl$cat, "SEX")
  expect_true(all(c("WT", "SEX") %in% names(tbl$data)))
})
test_that("continuous/categorical args override classification", {
  tbl <- .ferx_xpose_frame(fake_fit(), categorical = "WT")
  expect_false("WT" %in% tbl$cont)
  expect_true("WT" %in% tbl$cat)
})
test_that("unknown covariate override warns", {
  expect_warning(
    .ferx_xpose_frame(fake_fit(), continuous = "NOPE"),
    "not declared"
  )
})
test_that("covariate values come from covtab via LOCF when not echoed in sdtab", {
  tbl <- .ferx_xpose_frame(fake_fit())
  df <- tbl$data
  # WT is time-constant per subject in the fixture
  expect_equal(df$WT[df$ID == 1], c(70, 70))
  expect_equal(df$WT[df$ID == 2], c(80, 80))
})
test_that("covariate already echoed in sdtab is used directly (no covtab needed)", {
  # SEX (categorical, declared) has neither an echoed column nor a covtab here,
  # so it is dropped with a warning; this test only asserts the echoed WT.
  tbl <- suppressWarnings(
    .ferx_xpose_frame(fake_fit(with_covtab = FALSE, echo_wt = TRUE)))
  expect_true("WT" %in% names(tbl$data))
  expect_equal(tbl$data$WT, c(70, 70, 80, 80))
})
test_that("LOCF join carries the latest record at or before each obs time", {
  id  <- c(1, 1)
  t   <- c(2, 6)
  sid <- c(1, 1, 1)
  st  <- c(0, 3, 5)
  sv  <- c(10, 20, 30)
  expect_equal(.ferx_locf_join(id, t, sid, st, sv), c(10, 30))
})
test_that("LOCF tie-break carries the LAST source record at the latest time", {
  # two source records tied at t=3 with different values; true LOCF -> 99
  expect_equal(.ferx_locf_join(1, 5, c(1, 1, 1), c(0, 3, 3), c(10, 20, 99)), 99)
})
test_that("LOCF join matches across numeric vs character ID types", {
  # target IDs numeric (sdtab), source IDs character (covtab)
  out <- .ferx_locf_join(c(1, 2), c(5, 5),
                         c("1", "2"), c(0, 0), c(70, 80))
  expect_equal(out, c(70, 80))
})
test_that("ID join works when sdtab ID is numeric and ebe/ind ID is character", {
  tbl <- .ferx_xpose_frame(fake_fit())  # fixture now uses mismatched ID types
  df <- tbl$data
  expect_false(anyNA(df$ETA_CL))
  expect_false(anyNA(df$CL))
  expect_equal(df$CL[df$ID == 2], c(1.4, 1.4))
})
test_that("override naming an undeclared covariate (no covariate_types) is dropped + warned", {
  ff <- fake_fit()
  ff$covariate_types <- NULL
  ff$covtab <- NULL
  expect_warning(
    tbl <- .ferx_xpose_frame(ff, continuous = "WT"),
    "not declared"
  )
  expect_false("WT" %in% tbl$cont)
})
test_that("unrecognized covariate_types value is omitted with a warning", {
  ff <- fake_fit()
  ff$covariate_types <- c(WT = "continuous", SEX = NA, AGE = "Continuous")
  expect_warning(split <- .ferx_xpose_covariates(ff, NULL, NULL),
                 "unrecognized type")
  expect_false("SEX" %in% c(split$continuous, split$categorical))  # NA omitted
  expect_setequal(split$continuous, c("WT", "AGE"))                # case-insensitive
})
test_that("missing GOF column in sdtab warns at build time", {
  ff <- fake_fit()
  ff$sdtab$IPRED <- NULL
  expect_warning(.ferx_xpose_frame(ff), "no IPRED")
})
test_that("param/eta echoed into sdtab prefers the per-subject value, with a warning", {
  ff <- fake_fit()
  ff$sdtab$CL <- c(9, 9, 9, 9)   # stale per-observation CL echoed into sdtab
  expect_warning(tbl <- .ferx_xpose_frame(ff), "per-subject value")
  # per-subject value from individual_estimates wins (1.1 for ID 1, 1.4 for ID 2)
  expect_equal(tbl$data$CL[tbl$data$ID == 1], c(1.1, 1.1))
  expect_true("CL" %in% tbl$params)
})
test_that("new-xpose summary uses the 'run' key xpose reads (not 'runno')", {
  skip_if_not_installed("xpose")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose", runno = 7)
  expect_true("run" %in% xpdb$summary$label)
  expect_false("runno" %in% xpdb$summary$label)
})
test_that("roles map points at the expected columns", {
  roles <- .ferx_xpose_frame(fake_fit())$roles
  expect_equal(roles$id, "ID")
  expect_equal(roles$idv, "TIME")
  expect_equal(roles$dv, "DV")
  expect_setequal(roles$pred, c("PRED", "IPRED"))
  expect_setequal(roles$res,
                  c("RES", "IRES", "WRES", "CWRES", "IWRES", "NPDE", "NPD"))
  expect_setequal(roles$contcov, "WT")
  expect_setequal(roles$catcov, "SEX")
})
test_that("ferx_xpose rejects non-fit input", {
  expect_error(ferx_xpose(list()), "ferx_fit")
})
test_that("ferx_xpose errors clearly when the backend package is absent", {
  skip_if(requireNamespace("xpose", quietly = TRUE),
          "xpose installed; guard path not exercised")
  expect_error(ferx_xpose(fake_fit(), backend = "xpose"), "xpose")
})
test_that("backend 'xpose' builds an xpose_data with correct column roles", {
  skip_if_not_installed("xpose")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose")
  expect_s3_class(xpdb, "xpose_data")
  idx <- xpdb$data$index[[1]]
  role <- function(col) idx$type[idx$col == col]
  expect_equal(role("ID"), "id")
  expect_equal(role("TIME"), "idv")
  expect_equal(role("DV"), "dv")
  expect_equal(role("PRED"), "pred")
  expect_equal(role("IPRED"), "ipred")   # distinct from pred in xpose
  expect_equal(role("CWRES"), "res")
  expect_equal(role("NPDE"), "res")
  expect_equal(role("NPD"), "res")
  expect_equal(role("CL"), "param")
  expect_equal(role("ETA_CL"), "eta")
  expect_equal(role("WT"), "contcov")
  expect_equal(role("SEX"), "catcov")
})
test_that("backend 'xpose' GOF templates render", {
  skip_if_not_installed("xpose")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(suppressWarnings(print(xpose::dv_vs_ipred(xpdb))))
  expect_no_error(suppressWarnings(print(xpose::res_vs_idv(xpdb, res = "CWRES"))))
})
test_that("backend 'xpose4' builds an xpose.data with covariate/param vardefs", {
  skip_if_not_installed("xpose4")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose4")
  expect_s4_class(xpdb, "xpose.data")
  expect_setequal(xpdb@Prefs@Xvardef$covariates, c("WT", "SEX"))
  expect_setequal(xpdb@Prefs@Xvardef$parms, c("CL", "V", "ETA_CL", "ETA_V"))
  # categorical covariate coerced to factor in the data slot
  expect_s3_class(xpdb@Data$SEX, "factor")
})
test_that("backend 'xpose4' default plot renders", {
  skip_if_not_installed("xpose4")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose4")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(suppressWarnings(print(
    xpose4::xpose.plot.default("PRED", "DV", xpdb)
  )))
})

# --- estimation-iteration trace -> .ext / .grd files (issue #168) ------------

test_that("iteration files: ext + grd tables built from a widened trace", {
  files <- .ferx_xpose_iteration_files(fake_fit_trace())
  expect_s3_class(files, "data.frame")
  expect_setequal(files$extension, c("ext", "grd"))
  expect_setequal(names(files),
                  c("name", "extension", "problem", "subprob", "method",
                    "data", "modified"))

  ext <- files$data[[which(files$extension == "ext")]]
  # ITERATION, one column per declared parameter name, then OBJ.
  expect_equal(names(ext),
               c("ITERATION", "TVCL", "TVV", "ETA_CL", "PROP_ERR", "OBJ"))
  expect_equal(ext$ITERATION, 1:6)
  expect_equal(ext$OBJ, seq(100, by = -1, length.out = 6))

  grd <- files$data[[which(files$extension == "grd")]]
  # ITERATION + one GRD(n) column per parameter (NONMEM-style numeric names,
  # which grd_vs_iteration() keys its facet on).
  expect_equal(names(grd),
               c("ITERATION", "GRD(1)", "GRD(2)", "GRD(3)", "GRD(4)"))
})

test_that("iteration files: grd omitted when no gradients were recorded", {
  files <- .ferx_xpose_iteration_files(fake_fit_trace(grad = FALSE))
  expect_equal(files$extension, "ext")   # ext only, no grd table
})

test_that("iteration files: method label is the last non-empty trace method", {
  fit <- fake_fit_trace()
  fit$trace$method <- rep(c("gn", "focei"), each = 3L)
  files <- .ferx_xpose_iteration_files(fit)
  expect_true(all(files$method == "focei"))
})

test_that("iteration files: NULL when the fit has no trace", {
  expect_null(.ferx_xpose_iteration_files(fake_fit()))
})

test_that("iteration files: NULL for a scalar-only (pre-widening) trace", {
  fit <- fake_fit()
  fit$trace <- data.frame(iter = 1:3, method = "gn", ofv = c(3, 2, 1),
                          grad_norm = c(1, 0.5, 0.1), stringsAsFactors = FALSE)
  expect_null(.ferx_xpose_iteration_files(fit))
})

test_that("ferx_xpose validates the iterations flag", {
  expect_error(ferx_xpose(fake_fit(), iterations = NA), "iterations")
  expect_error(ferx_xpose(fake_fit(), iterations = "yes"), "iterations")
})

test_that("ferx_xpose populates / skips the files slot per iterations + trace", {
  skip_if_not_installed("xpose")
  # trace present + iterations = TRUE -> ext + grd populated
  xpdb <- ferx_xpose(fake_fit_trace(), iterations = TRUE)
  expect_false(is.null(xpdb$files))
  expect_setequal(xpdb$files$extension, c("ext", "grd"))
  # iterations = FALSE -> slot left empty
  expect_null(ferx_xpose(fake_fit_trace(), iterations = FALSE)$files)
  # no trace -> slot left empty even when iterations = TRUE (graceful no-op)
  expect_null(ferx_xpose(fake_fit(), iterations = TRUE)$files)
})

test_that("prm_vs_iteration / grd_vs_iteration render on a ferx_xpose trace", {
  skip_if_not_installed("xpose")
  xpdb <- ferx_xpose(fake_fit_trace())
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(suppressWarnings(print(xpose::prm_vs_iteration(xpdb))))
  expect_no_error(suppressWarnings(print(xpose::grd_vs_iteration(xpdb))))
})

# --- discrete (#900) and time-to-event (#901) sdtab rows (#688) --------------

# ferx-core writes every subject's Gaussian rows, then every subject's
# discrete rows, then every subject's TTE rows. The Xpose data keeps the
# Gaussian ones only.

# Fit `example` with `edit` (a function of the data frame) applied to its data
# and `model_edit` (a function of the model's lines) to its model.
xpose_edited_fit <- local({
  cache <- list()
  function(example, key, edit = identity, model_edit = identity) {
    id <- paste(example, key, sep = "/")
    if (is.null(cache[[id]])) {
      ex <- ferx_example(example)
      model <- tempfile(paste0("xpose_", key, "_"), fileext = ".ferx")
      writeLines(model_edit(readLines(ex$model)), model)
      d <- edit(utils::read.csv(ex$data, stringsAsFactors = FALSE))
      path <- tempfile(paste0("xpose_", key, "_"), fileext = ".csv")
      utils::write.csv(d, path, row.names = FALSE, quote = FALSE)
      cache[[id]] <<- suppressWarnings(
        ferx_fit(model, path, verbose = FALSE, covariance = FALSE,
                 settings = list(maxiter = 30L)))
    }
    cache[[id]]
  }
})

# Give the third subject the first subject's ID: a non-contiguous reuse, so
# an ID join hands the third subject the first one's etas and parameters.
reuse_third_id <- function(d) {
  subj <- .ferx_subject_index(d$ID)
  d$ID[subj == 3L] <- d$ID[match(1L, subj)]
  d
}

# Warfarin plus a binary endpoint on CMT 3: two binary records per subject,
# placed in time order among its concentrations.
add_binary_rows <- function(d) {
  subj <- .ferx_subject_index(d$ID)
  n <- max(subj)
  bin <- d[match(seq_len(n), subj)[rep(seq_len(n), each = 2L)], ]
  bin$TIME <- rep(c(12, 36), n)
  bin$DV <- as.character(rep(seq_len(n), each = 2L) %% 2L)
  bin$EVID <- 0L
  bin$AMT <- "."
  bin$CMT <- 3L
  bin$RATE <- 0L
  bin$MDV <- 0L
  out <- rbind(d, bin)
  s <- c(subj, rep(seq_len(n), each = 2L))
  out[order(s, out$TIME, -out$EVID), ]
}
add_binary_model <- function(m) {
  at <- grep("^\\s*sigma ", m)[1L]
  c(append(m, "  theta TH0(0.0, -10.0, 10.0)", after = at - 1L),
    "", "[binary_model]", "  cmt   = 3", "  logit = TH0")
}

# The checks every case shares: the Xpose data is exactly the leading
# sum(n_obs) sdtab rows, no TTE residual survives, and each row carries its
# own subject's etas / parameters (by position, never by ID).
expect_gaussian_only <- function(fit, n_dropped) {
  summ <- fit$subject_summary
  n_gauss <- sum(summ$n_obs)
  expect_equal(nrow(fit$sdtab), n_gauss + n_dropped)
  expect_message(tbl <- .ferx_xpose_frame(fit),
                 class = "ferx_xpose_rows_dropped_message")
  df <- tbl$data
  expect_equal(nrow(df), n_gauss)
  expect_equal(df$TIME, fit$sdtab$TIME[seq_len(n_gauss)])
  expect_equal(df$DV, fit$sdtab$DV[seq_len(n_gauss)])
  expect_false(any(.FITRX_TTE_COLS %in% names(df)))
  subject <- rep.int(seq_len(nrow(summ)), summ$n_obs)
  for (eta in setdiff(names(fit$ebe_etas), "ID")) {
    expect_equal(df[[eta]], fit$ebe_etas[[eta]][subject], info = eta)
  }
  for (p in setdiff(names(fit$individual_estimates), "ID")) {
    expect_equal(df[[p]], fit$individual_estimates[[p]][subject], info = p)
  }
  df
}

test_that("ferx_xpose keeps only pktte_joint's Gaussian rows (#688)", {
  fit <- xpose_edited_fit("pktte_joint", "plain")
  expect_equal(sum(fit$subject_summary$n_obs_tte), 6L)
  df <- expect_gaussian_only(fit, n_dropped = 6L)
  expect_equal(nrow(df), 12L)
  expect_false(any(c("COXSNELL", "MARTINGALE", "DEVIANCE", "PIT") %in% names(df)))
  expect_message(.ferx_xpose_frame(fit), "6 time-to-event sdtab rows")
  # message(<condition>) prints the text as is: the line ends in the text.
  m <- tryCatch(.ferx_xpose_frame(fit), message = conditionMessage)
  expect_match(m, "\n$")
})

test_that("ferx_xpose reports pktte_joint's nobs as its 12 observations (#688)", {
  skip_if_not_installed("xpose")
  fit <- xpose_edited_fit("pktte_joint", "plain")
  xpdb <- suppressMessages(ferx_xpose(fit))
  s <- xpdb$summary
  expect_equal(s$value[s$label == "nobs"], "12")
  expect_equal(s$value[s$label == "nind"], "6")
})

test_that("ferx_xpose joins pktte_joint etas by subject, not ID (#688)", {
  fit <- xpose_edited_fit("pktte_joint", "reused_id", edit = reuse_third_id)
  # The reuse leaves the third subject with etas of its own to pick up.
  expect_false(isTRUE(all.equal(fit$ebe_etas$ETA_CL[1], fit$ebe_etas$ETA_CL[3])))
  expect_gaussian_only(fit, n_dropped = 6L)
  # One subject per subject_summary row, though the data has five distinct IDs.
  s <- .ferx_xpose_new_summary(fit, 1L)
  expect_equal(s$value[s$label == "nind"], "6")
})

test_that("ferx_xpose joins by subject position when IDs share a number (#688)", {
  # Subject 2 is "0" + subject 1's ID: a distinct subject in the data, but
  # sdtab writes both IDs as the same number, so the run-length index merges
  # their rows and an ID join gives subject 2 subject 1's etas.
  fit <- xpose_edited_fit("warfarin", "zero_padded_id", edit = function(d) {
    subj <- .ferx_subject_index(d$ID)
    d$ID <- as.character(d$ID)
    d$ID[subj == 2L] <- paste0("0", d$ID[match(1L, subj)])
    d
  })
  summ <- fit$subject_summary
  n_obs <- summ$n_obs
  expect_equal(fit$sdtab$ID[n_obs[1] + 1L], fit$sdtab$ID[1])
  expect_false(isTRUE(all.equal(fit$ebe_etas$ETA_CL[1], fit$ebe_etas$ETA_CL[2])))
  expect_no_message(df <- .ferx_xpose_frame(fit)$data)
  subject <- rep.int(seq_len(nrow(summ)), n_obs)
  for (eta in setdiff(names(fit$ebe_etas), "ID")) {
    expect_equal(df[[eta]], fit$ebe_etas[[eta]][subject], info = eta)
  }
  for (p in setdiff(names(fit$individual_estimates), "ID")) {
    expect_equal(df[[p]], fit$individual_estimates[[p]][subject], info = p)
  }
})

# --- the covariate LOCF join keys on the subject (#695) ----------------------

# Warfarin with WT declared, and a WT that differs on every record: an
# observation that picks up any record but its own shows in the value.
add_time_varying_wt <- function(d) {
  d$WT <- 60 + 5 * .ferx_subject_index(d$ID) + d$TIME / 10
  d
}
declare_wt <- function(m) c(m, "", "[covariates]", "  WT continuous")

# Every observation record of `d`, whose WT the Xpose data should hold in order.
expect_own_wt <- function(fit, d) {
  obs <- d[d$EVID == 0L, ]
  df <- .ferx_xpose_frame(fit)$data
  expect_equal(df$TIME, obs$TIME)
  expect_equal(df$WT, obs$WT)
}

test_that("ferx_xpose carries covariates by subject when IDs share a number (#695)", {
  # Subject 2 is "0" + subject 1's ID, so sdtab has one run fewer than covtab,
  # and "01" itself is no sdtab ID: an ID join leaves subject 2 with no WT.
  edit <- function(d) {
    subj <- .ferx_subject_index(d$ID)
    d <- add_time_varying_wt(d)
    d$ID <- as.character(d$ID)
    d$ID[subj == 2L] <- paste0("0", d$ID[match(1L, subj)])
    d
  }
  fit <- xpose_edited_fit("warfarin", "zero_padded_id_wt", edit = edit,
                          model_edit = declare_wt)
  expect_equal(.ferx_n_subjects(fit$sdtab$ID), nrow(fit$subject_summary) - 1L)
  expect_equal(.ferx_n_subjects(fit$covtab$ID), nrow(fit$subject_summary))
  ex <- ferx_example("warfarin")
  expect_own_wt(fit, edit(utils::read.csv(ex$data, stringsAsFactors = FALSE)))
})

test_that("ferx_xpose carries covariates past a subject with no observations (#695)", {
  # Subject 2 keeps only its dose, so it has covtab records and no sdtab rows.
  # Subject 3 reuses subject 1's ID: an ID join puts their records in one
  # group, and subject 1's later observations take subject 3's WT.
  edit <- function(d) {
    subj <- .ferx_subject_index(d$ID)
    d <- reuse_third_id(add_time_varying_wt(d))
    d[subj != 2L | d$EVID != 0L, ]
  }
  fit <- xpose_edited_fit("warfarin", "dose_only_subject_wt", edit = edit,
                          model_edit = declare_wt)
  expect_equal(fit$subject_summary$n_obs[2], 0L)
  ex <- ferx_example("warfarin")
  expect_own_wt(fit, edit(utils::read.csv(ex$data, stringsAsFactors = FALSE)))
})

test_that("covariate join: covtab runs that are not the subjects fall back to the ID (#695)", {
  # A subject [data_selection] removed whole keeps its covtab records.
  fit <- fake_fit()
  fit$subject_summary <- data.frame(ID = c("1", "2"), n_obs = c(2L, 2L))
  fit$covtab <- rbind(fit$covtab[1:2, ],
                      data.frame(ID = "9", TIME = 0, EVID = 1, WT = 99, SEX = 1),
                      fit$covtab[3:4, ])
  df <- .ferx_xpose_frame(fit)$data
  expect_equal(df$WT, c(70, 70, 80, 80))
})

test_that("covariate join: without counts, sdtab's runs key a reused ID (#695)", {
  fit <- fake_fit()
  fit$covtab <- rbind(fit$covtab, transform(fit$covtab[1:2, ], WT = 90))
  fit$sdtab <- rbind(fit$sdtab, fit$sdtab[1:2, ])
  fit$individual_estimates <- NULL
  fit$ebe_etas <- NULL
  df <- .ferx_xpose_frame(fit)$data
  expect_equal(df$WT, c(70, 70, 80, 80, 90, 90))
})

test_that("ferx_xpose refuses a fit with no Gaussian observations (#688)", {
  for (case in list(c("binary_logistic", "discrete-endpoint"),
                    c("tte_exponential", "time-to-event"))) {
    fit <- xpose_edited_fit(case[1], "plain")
    expect_equal(sum(fit$subject_summary$n_obs), 0L)
    expect_error(.ferx_xpose_frame(fit), class = "ferx_xpose_no_gaussian_rows_error")
    expect_error(.ferx_xpose_frame(fit),
                 sprintf("no Gaussian observation rows \\(%d %s sdtab rows only\\)",
                         nrow(fit$sdtab), case[2]))
  }
})

test_that("ferx_xpose keeps only a binary model's Gaussian rows (#688)", {
  fit <- xpose_edited_fit("warfarin", "binary_reused_id",
                          edit = function(d) add_binary_rows(reuse_third_id(d)),
                          model_edit = add_binary_model)
  n_subj <- nrow(fit$subject_summary)
  expect_equal(sum(fit$subject_summary$n_obs_discrete), 2L * n_subj)
  expect_false(isTRUE(all.equal(fit$ebe_etas$ETA_CL[1], fit$ebe_etas$ETA_CL[3])))
  df <- expect_gaussian_only(fit, n_dropped = 2L * n_subj)
  expect_false(any(df$CMT == 3, na.rm = TRUE))
  expect_message(.ferx_xpose_frame(fit),
                 sprintf("%d discrete-endpoint sdtab rows", 2L * n_subj))
})

test_that("xpose rows: missing n_obs_discrete / n_obs_tte count as 0 (#688)", {
  fit <- fake_fit()
  fit$subject_summary <- data.frame(ID = c("1", "2"), n_obs = c(2L, 2L))
  expect_no_message(tbl <- .ferx_xpose_frame(fit))
  expect_equal(nrow(tbl$data), 4L)
  expect_equal(tbl$data$ETA_CL, fit$ebe_etas$ETA_CL[c(1, 1, 2, 2)])
})

test_that("xpose rows: without counts, the COXSNELL rows are left out (#688)", {
  fit <- fake_fit()
  tte <- fit$sdtab[c(1, 3), ]
  tte$TIME <- c(9, 9)
  tte$DV <- c(1, 0)
  fit$sdtab$COXSNELL <- fit$sdtab$MARTINGALE <- NA_real_
  fit$sdtab$DEVIANCE <- fit$sdtab$PIT <- NA_real_
  tte$COXSNELL <- c(0.4, 0.2)
  tte$MARTINGALE <- tte$DEVIANCE <- tte$PIT <- c(0.1, 0.1)
  fit$sdtab <- rbind(fit$sdtab, tte)
  expect_message(tbl <- .ferx_xpose_frame(fit), "2 time-to-event sdtab rows")
  expect_equal(nrow(tbl$data), 4L)
  expect_false("COXSNELL" %in% names(tbl$data))
  expect_equal(tbl$data$ETA_CL, fit$ebe_etas$ETA_CL[c(1, 1, 2, 2)])
})

test_that("xpose rows: the rows-left-out note is an ordinary message (#702)", {
  fit <- fake_fit()
  tte <- fit$sdtab[c(1, 3), ]
  tte$TIME <- c(9, 9)
  for (col in .FITRX_TTE_COLS) {
    fit$sdtab[[col]] <- NA_real_
    tte[[col]] <- c(0.4, 0.2)
  }
  fit$sdtab <- rbind(fit$sdtab, tte)
  cond <- tryCatch(.ferx_xpose_frame(fit), message = identity)
  expect_s3_class(cond, c("ferx_xpose_rows_dropped_message", "message",
                          "condition"), exact = TRUE)
  expect_match(conditionMessage(cond), "2 time-to-event sdtab rows")
  expect_null(conditionCall(cond))
  expect_silent(suppressMessages(.ferx_xpose_frame(fit)))
})

test_that("xpose rows: without counts, rows with no N_OBS are left out (#688)", {
  fit <- fake_fit()
  fit$sdtab$N_OBS <- 2
  disc <- fit$sdtab[c(1, 3, 3), ]
  disc$TIME <- c(6, 6, 7)
  disc$DV <- c(1, 0, 1)
  disc$N_OBS <- NaN
  tte <- disc[1:2, ]
  tte$TIME <- c(9, 9)
  for (col in .FITRX_TTE_COLS) {
    fit$sdtab[[col]] <- NA_real_
    disc[[col]] <- NA_real_
    tte[[col]] <- c(0.4, 0.2)
  }
  fit$sdtab <- rbind(fit$sdtab, disc, tte)
  expect_message(tbl <- .ferx_xpose_frame(fit),
                 "3 discrete-endpoint and 2 time-to-event sdtab rows")
  expect_equal(nrow(tbl$data), 4L)
  expect_equal(tbl$data$TIME, c(1, 4, 1, 4))
  expect_false(any(.FITRX_TTE_COLS %in% names(tbl$data)))
})
