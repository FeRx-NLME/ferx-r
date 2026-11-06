# Which subject owns each sdtab row in a .fitrx bundle (#540), and the
# per-subject ofv_contribution / n_obs in ebes.csv (#554).
#
# Both used to be read back from sdtab: rows by the run length of the ID
# number (the ID parsed as one, else the subject's position), values from the
# first row of each run. Two adjacent subjects with the same number are one
# run, a subject without observations has no rows, and a discrete (binary)
# row has a blank N_OBS. They now come from `fit$subject_summary`, the
# engine's per-subject result. The oracle is the data the fit read (row counts
# per subject), and ferx-core's loader, which reads every bundle here.

# `ferx_fit()` without a covariance step, muffling only the notes that a
# call-time argument overrides the model file's [fit_options].
quiet_fit <- function(model, data, maxiter = NULL) {
  settings <- if (!is.null(maxiter)) list(maxiter = maxiter) else list()
  withCallingHandlers(
    ferx_fit(model, data, verbose = FALSE, covariance = FALSE,
             settings = settings),
    warning = function(w) {
      if (grepl("ferx_fit() argument overrides it", conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# A table with its all-empty columns as doubles. A column with no value in a
# bundle reads back as logical NA where the live one holds NaN, as it always
# has; the rows are what is under test here.
as_doubles <- function(df) {
  df[] <- lapply(df, function(v) if (is.logical(v)) as.numeric(v) else v)
  df
}

# Fit `example` on its data after `edit` (a function of the data frame), and
# its model after `model_edit` (a function of the model's lines).
edited_fit <- local({
  cache <- list()
  function(example, key, edit, model_edit = NULL) {
    id <- paste(example, key, sep = "/")
    if (is.null(cache[[id]])) {
      ex <- ferx_example(example)
      if (!is.null(model_edit)) {
        model <- tempfile(paste0("edited_", key, "_"), fileext = ".ferx")
        writeLines(model_edit(readLines(ex$model)), model)
        ex$model <- model
      }
      d <- edit(utils::read.csv(ex$data, stringsAsFactors = FALSE))
      path <- tempfile(paste0("edited_", key, "_"), fileext = ".csv")
      utils::write.csv(d, path, row.names = FALSE,
                       quote = if (is.character(d$ID)) match("ID", names(d)) else FALSE)
      fit <- quiet_fit(ex$model, path, maxiter = 30L)
      cache[[id]] <<- list(fit = fit, data = d)
    }
    cache[[id]]
  }
})

# Warfarin's subjects in data order, relabelled by position: `labels[k]`
# replaces the k-th subject's ID; `no_obs` positions lose their observations
# (the doses stay, so each is still a subject).
warfarin_relabelled <- function(key, labels, no_obs = integer()) {
  edited_fit("warfarin", key, function(d) {
    subj <- ferx:::.ferx_subject_index(d$ID)
    keep <- !(subj %in% no_obs & d$EVID == 0)
    d$ID <- as.character(d$ID)
    for (k in names(labels)) d$ID[subj == as.integer(k)] <- labels[[k]]
    d[keep, , drop = FALSE]
  })
}

# "a" (number 1, its position) then "1": adjacent subjects with the same
# sdtab number, so their rows are one run. A later "1" repeats an ID
# non-contiguously.
repeat_case <- function() {
  warfarin_relabelled("repeat", c("1" = "a", "2" = "1", "4" = "1"))
}
# "p", then "q" without observations, then "p" again: the two "p" blocks are
# adjacent in predictions.csv and the text cannot tell them apart, while
# their sdtab numbers are their positions, 1 and 3.
gap_case <- function() {
  warfarin_relabelled("gap", c("1" = "p", "2" = "q", "3" = "p"), no_obs = 2L)
}

# Warfarin plus a binary endpoint on CMT 3: three binary records per subject
# between its concentrations, subject 2's ID text. Each subject has rows of
# both kinds, so sdtab's order - every subject's continuous rows, then every
# subject's binary rows - differs from subject-by-subject order (#561).
joint_case <- function() {
  edited_fit("warfarin", "joint_binary",
    edit = function(d) {
      subj <- ferx:::.ferx_subject_index(d$ID)
      n <- max(subj)
      bin <- d[match(seq_len(n), subj)[rep(seq_len(n), each = 3L)], ]
      k <- rep(0:2, n)
      bin$TIME <- c(12, 36, 72)[k + 1L]
      bin$DV <- as.character((rep(seq_len(n), each = 3L) + k) %% 2L)
      bin$EVID <- 0L
      bin$AMT <- "."
      bin$CMT <- 3L
      bin$RATE <- 0L
      bin$MDV <- 0L
      out <- rbind(d, bin)
      s <- c(subj, rep(seq_len(n), each = 3L))
      out <- out[order(s, out$TIME, -out$EVID), ]
      out$ID <- ifelse(out$ID == 2L, "b", as.character(out$ID))
      out
    },
    model_edit = function(m) {
      at <- grep("^\\s*sigma ", m)[1L]
      c(append(m, "  theta TH0(0.0, -10.0, 10.0)", after = at - 1L),
        "", "[binary_model]", "  cmt   = 3", "  logit = TH0")
    }
  )
}

# Per-subject continuous observation count, from the data the fit read.
data_n_obs <- function(d) {
  subj <- ferx:::.ferx_subject_index(d$ID)
  obs <- d$EVID == 0 & d$MDV == 0
  as.integer(tabulate(subj[obs], nbins = max(subj)))
}

# A bundle entry, read as text.
bundle_csv <- function(path, entry) {
  staging <- tempfile("fitrx_entry_")
  dir.create(staging)
  utils::unzip(path, files = entry, exdir = staging)
  utils::read.csv(file.path(staging, entry), stringsAsFactors = FALSE,
                  check.names = FALSE, colClasses = "character",
                  na.strings = character())
}

test_that("fit$subject_summary holds the engine's per-subject values (#554)", {
  for (case in list(repeat_case(), gap_case())) {
    fit <- case$fit
    summ <- fit$subject_summary
    expect_identical(names(summ),
                     c("ID", "ofv_contribution", "n_obs", "n_obs_discrete", "n_obs_tte"))
    expect_identical(summ$ID, fit$subject_ids)
    expect_identical(summ$ID, fit$ebe_etas$ID)
    expect_identical(summ$n_obs, data_n_obs(case$data))
    expect_identical(summ$n_obs_discrete, integer(nrow(summ)))
    expect_identical(summ$n_obs_tte, integer(nrow(summ)))
    # Each subject's sdtab rows, counted off from the data, carry its value.
    owner <- rep.int(seq_len(nrow(summ)), summ$n_obs)
    expect_identical(fit$sdtab$EBE_OFV, summ$ofv_contribution[owner])
  }
  # The cases are what they claim: the repeat case's first two subjects share
  # an sdtab number, the gap case has a subject without rows.
  expect_identical(unique(repeat_case()$fit$sdtab$ID)[1:2], c(1, 3))
  gap <- gap_case()$fit
  expect_identical(gap$subject_summary$n_obs[2], 0L)
  expect_identical(unique(gap$sdtab$ID)[1:2], c(1, 3))
})

test_that("ebes.csv carries each subject's own ofv_contribution and n_obs (#540)", {
  for (case in list(repeat_case(), gap_case())) {
    fit <- case$fit
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    ebes <- bundle_csv(path, "ebes.csv")
    expect_identical(ebes$ID, fit$subject_ids)
    expect_identical(as.integer(ebes$n_obs), data_n_obs(case$data))
    expect_identical(ferx:::.fitrx_parse_doubles(ebes$ofv_contribution),
                     fit$subject_summary$ofv_contribution)
    # No discrete rows: the engine's own ebes.csv columns, nothing added.
    expect_false("n_obs_discrete" %in% names(ebes))
    preds <- bundle_csv(path, "predictions.csv")
    expect_identical(preds$ID, rep.int(fit$subject_ids, fit$subject_summary$n_obs))
  }
})

test_that("a bundle with repeated and gapped IDs reloads every row on its subject (#540)", {
  for (case in list(repeat_case(), gap_case())) {
    fit <- case$fit
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    loaded <- ferx_load_fit(path)
    expect_identical(loaded$sdtab$ID, fit$sdtab$ID)
    expect_identical(loaded$subject_summary, fit$subject_summary)
    expect_equal(loaded$sdtab, fit$sdtab)
    expect_identical(loaded$ebe_etas$ID, fit$ebe_etas$ID)
  }
})

# binary_logistic with its second subject's ID text: three binary records per
# subject, no continuous observations, no etas.
binary_case <- function() {
  edited_fit("binary_logistic", "text_b", function(d) {
    d$ID <- ifelse(d$ID == 2L, "b", as.character(d$ID))
    d
  })
}

test_that("a binary fit's bundle loads in ferx-core, discrete rows apart (#554)", {
  fit <- binary_case()$fit
  summ <- fit$subject_summary
  expect_identical(summ$ID[1:3], c("1", "b", "3"))
  # Binary only: no continuous observations, three binary records each.
  expect_identical(summ$n_obs, integer(nrow(summ)))
  expect_identical(summ$n_obs_discrete, rep(3L, nrow(summ)))
  expect_identical(nrow(fit$sdtab), 3L * nrow(summ))

  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  entries <- utils::unzip(path, list = TRUE)$Name
  expect_true("predictions_discrete.csv" %in% entries)
  ebes <- bundle_csv(path, "ebes.csv")
  expect_identical(ebes$n_obs, rep("0", nrow(summ)))
  expect_identical(ebes$n_obs_discrete, rep("3", nrow(summ)))
  expect_identical(ebes$ID, summ$ID)
  expect_identical(ferx:::.fitrx_parse_doubles(ebes$ofv_contribution),
                   summ$ofv_contribution)
  expect_identical(nrow(bundle_csv(path, "predictions.csv")), 0L)
  disc <- bundle_csv(path, "predictions_discrete.csv")
  expect_identical(disc$ID, rep(summ$ID, each = 3L))

  loaded <- ferx_load_fit(path)
  expect_identical(loaded$sdtab$ID, fit$sdtab$ID)
  expect_equal(as_doubles(loaded$sdtab), fit$sdtab)
  expect_identical(loaded$subject_summary, summ)
  # n_obs_discrete is a per-subject column, not an eta: still no etas.
  expect_null(fit$ebe_etas)
  expect_null(loaded$ebe_etas)
})

test_that("a fit with continuous and binary rows keeps both on their subjects (#561)", {
  case <- joint_case()
  fit <- case$fit
  summ <- fit$subject_summary
  d <- case$data
  pk <- d[d$CMT != 3L, ]
  expect_identical(summ$n_obs, data_n_obs(pk))
  expect_identical(summ$n_obs_discrete, rep(3L, nrow(summ)))
  expect_identical(summ$ID[1:2], c("1", "b"))
  ids <- summ$ID

  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  expect_identical(bundle_csv(path, "predictions.csv")$ID, rep.int(ids, summ$n_obs))
  expect_identical(bundle_csv(path, "predictions_discrete.csv")$ID,
                   rep(ids, each = 3L))

  loaded <- ferx_load_fit(path)
  expect_identical(loaded$sdtab$ID, fit$sdtab$ID)
  expect_equal(as_doubles(loaded$sdtab), fit$sdtab)
  expect_identical(loaded$subject_summary, summ)
  # The EBE table has the live shape: the per-subject columns are in
  # subject_summary, not counted as etas.
  expect_identical(names(loaded$ebe_etas), names(fit$ebe_etas))

  live <- ferx_covariance(fit)
  skip_if(is.null(live$cov_matrix), "covariance step did not converge - skipping")
  expect_identical(ferx_covariance(loaded)$se_theta, live$se_theta)
})

test_that("a bundle ferx-core wrote alone rebuilds its binary rows as core's sdtab fills them", {
  skip_if(!nzchar(Sys.which("zip")), "needs a zip binary")
  fit <- joint_case()$fit
  summ <- fit$subject_summary
  saved <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, saved)
  # No sdtab_endpoint_rows.csv: the binary rows are only in core's
  # predictions_discrete.csv, and the glue rebuilds them.
  path <- fitrx_core_only_bundle(saved)

  loaded <- ferx_load_fit(path)
  expect_identical(nrow(loaded$sdtab), nrow(fit$sdtab))
  expect_identical(loaded$subject_summary$n_obs_discrete, summ$n_obs_discrete)
  disc <- sum(summ$n_obs) + seq_len(sum(summ$n_obs_discrete))
  # Every column the bundle carries, EBE_OFV included, on the binary rows.
  for (col in names(loaded$sdtab)) {
    expect_identical(as_doubles(loaded$sdtab)[[col]][disc], fit$sdtab[[col]][disc],
                     info = col)
  }
})

# fixtures/mixture_override.ferx on mixture_iv.csv with three binary records
# per subject on CMT 3, and a per-subject data covariate called MIXNUM that the
# model's [output] block lists: a user column, blank on the binary rows as
# every [output] column is. The engine's own MIXNUM (the class a subject was
# simulated under, ferx-core #2372) is only on a fit of data drawn from a
# [simulation] block, which ferx_fit() never makes.
mixnum_output_case <- local({
  case <- NULL
  function() {
    if (is.null(case)) {
      d <- utils::read.csv(test_path("fixtures", "mixture_iv.csv"))
      ids <- unique(d$ID)
      n <- length(ids)
      bin <- d[match(ids, d$ID)[rep(seq_len(n), each = 3L)], ]
      k <- rep(0:2, n)
      bin$TIME <- c(2, 6, 12)[k + 1L]
      bin$DV <- (rep(seq_len(n), each = 3L) + k) %% 2L
      bin$AMT <- 0
      bin$EVID <- 0L
      bin$CMT <- 3L
      d <- rbind(d, bin)
      d <- d[order(match(d$ID, ids), d$TIME, -d$EVID), ]
      d$MIXNUM <- rep_len(1:2, n)[match(d$ID, ids)]
      data <- tempfile(fileext = ".csv")
      utils::write.csv(d, data, row.names = FALSE, quote = FALSE)

      m <- readLines(test_path("fixtures", "mixture_override.ferx"))
      at <- grep("^\\s*omega ", m)[1L]
      m <- c(append(m, "  theta TH0(0.0, -10.0, 10.0)", after = at - 1L),
             "", "[binary_model]", "  cmt   = 3", "  logit = TH0",
             "", "[output]", "  MIXNUM")
      model <- tempfile(fileext = ".ferx")
      writeLines(m, model)
      fit <- muffle_mixl_lint(ferx_fit(model, data, method = "focei",
                                       covariance = FALSE, verbose = FALSE))
      case <<- list(fit = fit, data = d)
    }
    case
  }
})

test_that("an [output] MIXNUM column stays blank on binary rows rebuilt from core's", {
  skip_if(!nzchar(Sys.which("zip")), "needs a zip binary")
  fit <- mixnum_output_case()$fit
  summ <- fit$subject_summary
  n_gauss <- sum(summ$n_obs)
  disc <- n_gauss + seq_len(sum(summ$n_obs_discrete))
  expect_identical(summ$n_obs_discrete, rep(3L, nrow(summ)))
  # The live fit: the covariate on the Gaussian rows, both values among them,
  # and blank on the binary rows.
  expect_setequal(unique(fit$sdtab$MIXNUM[seq_len(n_gauss)]), c(1, 2))
  expect_true(all(is.nan(fit$sdtab$MIXNUM[disc])))

  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  # The R package's wide predictions.csv and core's predictions_discrete.csv
  # without the glue's sdtab_endpoint_rows.csv: the shape of every bundle the
  # R writer of ferx 0.4.0 and before saved. The glue rebuilds the binary rows
  # from core's, and used to fill MIXNUM on them from the Gaussian rows.
  status <- utils::zip(path, "sdtab_endpoint_rows.csv", flags = "-qd")
  skip_if(status != 0L, "zip -d failed")
  expect_false("sdtab_endpoint_rows.csv" %in% utils::unzip(path, list = TRUE)$Name)
  expect_true("MIXNUM" %in% names(bundle_csv(path, "predictions.csv")))

  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_identical(nrow(loaded$sdtab), nrow(fit$sdtab))
  expect_identical(loaded$sdtab$MIXNUM, fit$sdtab$MIXNUM)
})

test_that("ebes.csv puts ofv_contribution and n_obs before a mixture's columns (#561)", {
  # ferx-core reads the two by position right after the etas, and MIXEST /
  # PMIX_* by name; a reloaded mixture bundle carries those in fit$ebe_etas.
  fit <- mixture_override_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, path))
  cols <- names(bundle_csv(path, "ebes.csv"))
  etas <- fit$eta_names
  expect_identical(cols[seq_len(length(etas) + 3L)],
                   c("ID", etas, "ofv_contribution", "n_obs"))
  expect_true("MIXEST" %in% cols)
  expect_gt(match("MIXEST", cols), match("n_obs", cols))
})

test_that("a TTE fit's bundle keeps its time-to-event rows (#554)", {
  ex <- ferx_example("tte_exponential")
  fit <- quiet_fit(ex$model, ex$data)
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$subject_summary, fit$subject_summary)
  expect_identical(loaded$sdtab, fit$sdtab)
})

test_that("a no-eta fit's ebes.csv has the engine's columns (#554)", {
  ex <- ferx_example("one_cpt_iv_pooled")
  fit <- quiet_fit(ex$model, ex$data)
  expect_null(fit$ebe_etas)
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  ebes <- bundle_csv(path, "ebes.csv")
  expect_identical(names(ebes), c("ID", "ofv_contribution", "n_obs"))
  expect_identical(ebes$ID, fit$subject_ids)
  expect_identical(as.integer(ebes$n_obs), fit$subject_summary$n_obs)
  expect_identical(ferx:::.fitrx_parse_doubles(ebes$ofv_contribution),
                   fit$subject_summary$ofv_contribution)
})
