# Which subject owns each sdtab row in a .fitrx bundle (#540), and the
# per-subject ofv_contribution / n_obs in ebes.csv (#554).
#
# Both used to be read back from sdtab: rows by the run length of the ID
# number (the ID parsed as one, else the subject's position), values from the
# first row of each run. Two adjacent subjects with the same number are one
# run, a subject without observations has no rows, and a discrete (binary)
# row has a blank N_OBS. They now come from `fit$subject_summary`, the
# engine's per-subject result. The oracles are the data the fit read (row
# counts per subject) and ferx-core's own loader (`fitrx_engine_trip()`),
# never R's loader agreeing with R's writer.

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

# ferx-core's .fitrx writer prints doubles with six decimals; R's is exact.
core_tol <- 1e-6

# A table with its all-empty columns as doubles. A column with no value in a
# bundle reads back as logical NA where the live one holds NaN, as it always
# has; the rows are what is under test here.
as_doubles <- function(df) {
  df[] <- lapply(df, function(v) if (is.logical(v)) as.numeric(v) else v)
  df
}

# A fit reloaded from the bundle ferx-core wrote after reading R's carries
# the same per-subject values: IDs and counts exactly, ofv to six decimals.
expect_core_summary <- function(from_core, fit) {
  expect_identical(from_core$subject_summary[c("ID", "n_obs")],
                   fit$subject_summary[c("ID", "n_obs")])
  expect_equal(from_core$subject_summary$ofv_contribution,
               fit$subject_summary$ofv_contribution, tolerance = core_tol)
}

# Fit `example` on its data after `edit` (a function of the data frame).
edited_fit <- local({
  cache <- list()
  function(example, key, edit) {
    id <- paste(example, key, sep = "/")
    if (is.null(cache[[id]])) {
      ex <- ferx_example(example)
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
                     c("ID", "ofv_contribution", "n_obs", "n_obs_discrete"))
    expect_identical(summ$ID, fit$subject_ids)
    expect_identical(summ$ID, fit$ebe_etas$ID)
    expect_identical(summ$n_obs, data_n_obs(case$data))
    expect_identical(summ$n_obs_discrete, integer(nrow(summ)))
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

    trip <- fitrx_engine_trip(fit)
    from_core <- ferx_load_fit(trip$core_path)
    expect_identical(from_core$sdtab$ID, fit$sdtab$ID)
    expect_identical(from_core$ebe_etas$ID, fit$ebe_etas$ID)
    expect_core_summary(from_core, fit)
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

  # The engine reads the R bundle: ebes.csv's n_obs is a count, and
  # predictions.csv has no row past it.
  trip <- fitrx_engine_trip(fit)
  core_ebes <- bundle_csv(trip$core_path, "ebes.csv")
  expect_identical(core_ebes$ID, summ$ID)
  expect_identical(core_ebes$n_obs, ebes$n_obs)
  expect_equal(as.numeric(core_ebes$ofv_contribution), summ$ofv_contribution,
               tolerance = core_tol)
})

test_that("a TTE fit's bundle loads in ferx-core (#554)", {
  ex <- ferx_example("tte_exponential")
  fit <- quiet_fit(ex$model, ex$data)
  trip <- fitrx_engine_trip(fit)
  from_core <- ferx_load_fit(trip$core_path)
  expect_core_summary(from_core, fit)
  expect_identical(from_core$sdtab$ID, fit$sdtab$ID)
  loaded <- ferx_load_fit(trip$r_path)
  expect_identical(loaded$sdtab$ID, fit$sdtab$ID)
})

test_that("a no-eta fit's ebes.csv is the engine's, up to its six decimals (#554)", {
  # Byte for byte only where ferx-core's writer is exact: it prints ofv to six
  # decimals, R to the last bit. Header, IDs and n_obs are the same text.
  ex <- ferx_example("one_cpt_iv_pooled")
  fit <- quiet_fit(ex$model, ex$data)
  expect_null(fit$ebe_etas)
  trip <- fitrx_engine_trip(fit)
  r_ebes <- bundle_csv(trip$r_path, "ebes.csv")
  core_ebes <- bundle_csv(trip$core_path, "ebes.csv")
  expect_identical(names(r_ebes), c("ID", "ofv_contribution", "n_obs"))
  expect_identical(names(r_ebes), names(core_ebes))
  expect_identical(r_ebes[c("ID", "n_obs")], core_ebes[c("ID", "n_obs")])
  expect_identical(r_ebes$ID, fit$subject_ids)
  expect_equal(as.numeric(core_ebes$ofv_contribution),
               as.numeric(r_ebes$ofv_contribution), tolerance = core_tol)
})

test_that("a bundle written before fit$subject_summary still loads (#540, #554)", {
  # The old writer's shape: no subject_summary, so ebes.csv's per-subject
  # columns are empty and every sdtab row, discrete ones too, is in
  # predictions.csv. The loader walks the ID text instead.
  for (fit in list(repeat_case()$fit, binary_case()$fit)) {
    old <- fit
    old$subject_summary <- NULL
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(old, path)
    expect_false("predictions_discrete.csv" %in% utils::unzip(path, list = TRUE)$Name)
    expect_true(all(bundle_csv(path, "ebes.csv")$n_obs == ""))
    loaded <- ferx_load_fit(path)
    expect_identical(loaded$sdtab$ID, fit$sdtab$ID)
    expect_equal(as_doubles(loaded$sdtab), fit$sdtab)
    expect_true(all(is.na(loaded$subject_summary$n_obs)))
  }
})
