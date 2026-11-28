# Every ID column ferx returns is the subject's ID text, verbatim (#786).
#
# `fit$sdtab$ID` used to be the engine's number for the subject: the ID parsed
# as a number, else its 1-based position. So "S1" came back as 1, "001" as 1,
# subjects "1" and "1.0" merged, and a join of `ferx_simulate(fit = )` (text
# IDs) onto `fit$sdtab` matched nothing or stopped on a type mismatch. Now
# every ID column is character and holds what `fit$subject_ids` holds.
#
# The relabelled fits are twins of the bundled warfarin fit: the same data
# under other labels, so their numbers agree exactly (helper-relabelled-fit.R).

id_fit <- function(key) {
  switch(key,
    numeric = warfarin_fit(),
    letters = relabelled_fit("warfarin", "id786_letters", function(id) paste0("S", id),
                             covariance = FALSE),
    padded  = relabelled_fit("warfarin", "id786_padded", function(id) sprintf("%03d", id),
                             covariance = FALSE)
  )
}
id_keys <- c("numeric", "letters", "padded")

# A cached warfarin fit on the bundled data after `edit` (IDs read as text).
edited_warfarin <- local({
  cache <- list()
  function(key, edit) {
    if (is.null(cache[[key]])) {
      ex <- ferx_example("warfarin")
      d <- utils::read.csv(ex$data, stringsAsFactors = FALSE, colClasses = c(ID = "character"))
      d <- edit(d)
      path <- tempfile(paste0("id786_", key, "_"), fileext = ".csv")
      utils::write.csv(d, path, row.names = FALSE, quote = match("ID", names(d)))
      cache[[key]] <<- ferx_fit(ex$model, path, method = "focei", verbose = FALSE,
                                covariance = FALSE, settings = list(maxiter = 30L))
    }
    cache[[key]]
  }
})

# The data file a fit was made from, as the fit records it.
id_data <- function(fit) fit$data_path

# Every data frame column named `ID` anywhere in `x`, as `path = typeof`: list
# elements and attributes alike, since some tables ride as an attribute
# (`ferx_simulate()`'s `parameters`, which `ferx_sim_parameters()` reads).
id_column_types <- function(x, path = "x") {
  out <- character()
  for (a in setdiff(names(attributes(x)), c("names", "row.names", "class", "dim", "dimnames"))) {
    value <- attr(x, a, exact = TRUE)
    if (is.list(value)) out <- c(out, id_column_types(value, paste0(path, "@", a)))
  }
  if (is.data.frame(x)) {
    if ("ID" %in% names(x)) out[path] <- typeof(x[["ID"]])
    x <- unclass(x)
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      out <- c(out, id_column_types(x[[i]], paste0(path, "$", nms[i])))
    }
  }
  out
}

expect_all_ids_character <- function(x, label) {
  types <- id_column_types(x, label)
  expect_gt(length(types), 0L)
  bad <- types[types != "character"]
  expect(length(bad) == 0L,
         sprintf("ID columns that are not character: %s",
                 paste(names(bad), bad, sep = " = ", collapse = ", ")))
}

test_that("ferx_simulate(fit = ) joins onto fit$sdtab by ID and TIME for every ID form (#786)", {
  skip_on_cran()
  for (key in id_keys) {
    fit <- id_fit(key)
    sim <- ferx_simulate(fit$model_path, id_data(fit), n_sim = 1L, seed = 1L, fit = fit)
    sim <- sim[sim$TIME %in% fit$sdtab$TIME, , drop = FALSE]
    expect_identical(typeof(sim$ID), "character", label = key)
    expect_identical(typeof(fit$sdtab$ID), "character", label = key)
    # What `dplyr::left_join()` does: it refuses mismatched key types, and keeps
    # every left row. Base merge() would coerce the types instead.
    joined <- merge(sim, fit$sdtab[, c("ID", "TIME", "PRED")],
                    by = c("ID", "TIME"), all.x = TRUE, sort = FALSE)
    expect_identical(nrow(joined), nrow(sim), label = key)
    expect_false(anyNA(joined$PRED), label = key)
  }
})

test_that("fit$sdtab$ID holds fit$subject_ids, never a position (#786)", {
  skip_on_cran()
  for (key in id_keys) {
    fit <- id_fit(key)
    expect_identical(unique(fit$sdtab$ID), fit$subject_ids, label = key)
    expect_identical(fit$ebe_etas$ID, fit$subject_ids, label = key)
  }
  expect_identical(id_fit("padded")$subject_ids[1:2], c("001", "002"))
  expect_identical(id_fit("letters")$subject_ids[1:2], c("S1", "S2"))
})

test_that("subjects 1 and 1.0 stay two subjects in fit$sdtab (#786)", {
  skip_on_cran()
  fit <- edited_warfarin("one_point_zero", function(d) {
    d$ID[d$ID == "2"] <- "1.0"
    d
  })
  expect_identical(fit$subject_ids[1:2], c("1", "1.0"))
  expect_identical(length(unique(fit$sdtab$ID)), length(fit$subject_ids))
  expect_identical(unique(fit$sdtab$ID), fit$subject_ids)
  expect_identical(sum(fit$sdtab$ID == "1.0"), sum(fit$sdtab$ID == "1"))
})

test_that("a dose-only subject is in ebe_etas, not sdtab, and IDs stay character (#786)", {
  skip_on_cran()
  fit <- edited_warfarin("dose_only", function(d) {
    dose <- d[d$ID == "1" & d$EVID == 1, , drop = FALSE]
    dose$ID <- "X"
    rbind(d, dose)
  })
  expect_true("X" %in% fit$subject_ids)
  expect_true("X" %in% fit$ebe_etas$ID)
  expect_false("X" %in% fit$sdtab$ID)
  expect_identical(unique(fit$sdtab$ID), setdiff(fit$subject_ids, "X"))
  expect_all_ids_character(unclass(fit), "fit")
})

test_that("a numeric-ID fit's sdtab is its relabelled twin's, ID aside (#786)", {
  skip_on_cran()
  num <- id_fit("numeric")
  # The bundled IDs are 1..10, so the text reads back as the old number.
  n_obs <- num$subject_summary$n_obs
  expect_identical(as.numeric(num$sdtab$ID),
                   rep(as.numeric(num$subject_ids), n_obs))
  for (key in c("letters", "padded")) {
    twin <- id_fit(key)
    expect_identical(num$sdtab[names(num$sdtab) != "ID"],
                     twin$sdtab[names(twin$sdtab) != "ID"], label = key)
    expect_identical(match(twin$sdtab$ID, twin$subject_ids),
                     match(num$sdtab$ID, num$subject_ids), label = key)
  }
})

test_that("every ID column survives ferx_save_fit() / ferx_load_fit() as text (#786)", {
  skip_on_cran()
  for (key in id_keys) {
    fit <- id_fit(key)
    path <- tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    back <- ferx_load_fit(path)
    before <- id_column_types(unclass(fit), "fit")
    after <- id_column_types(unclass(back), "fit")
    # A loaded fit carries a subset of the frames (no individual_estimates
    # without the data); each one it has keeps its text IDs.
    expect_true(all(c("fit$sdtab", "fit$ebe_etas") %in% names(after)), label = key)
    expect_identical(after, before[names(after)], label = key)
    for (p in names(after)) {
      get_col <- function(f) eval(parse(text = sub("^fit", "f", p)))$ID
      expect_identical(get_col(unclass(back)), get_col(unclass(fit)), label = paste(key, p))
    }
  }
})

# -- the sweep: every returned frame, one fixture per family -------------------
# The guard (tools/check-glue-id-columns.sh) holds every `("ID", ..)` column;
# the sweep holds the two tables built from parallel name/value vectors, which
# the guard cannot see: `sim_parameters_to_df` (`ferx_sim_parameters()`) and
# `individual_estimates`.

test_that("every ID column is character: fit, simulate, predict, npde (#786)", {
  skip_on_cran()
  fit <- id_fit("letters")
  expect_all_ids_character(unclass(fit), "fit")
  sim <- ferx_simulate(fit$model_path, id_data(fit), n_sim = 1L, seed = 1L, fit = fit)
  expect_all_ids_character(list(sim = sim), "simulate")
  expect_true("simulate$sim@parameters" %in% names(id_column_types(list(sim = sim), "simulate")))
  expect_identical(typeof(ferx_sim_parameters(sim)$ID), "character")
  expect_true("fit$individual_estimates" %in% names(id_column_types(unclass(fit), "fit")))
  pred <- ferx_predict(fit$model_path, id_data(fit), fit = fit)
  expect_all_ids_character(list(pred = pred), "predict")
  post <- suppressWarnings(ferx_calc_npde(fit, nsim = 50L, seed = 1L))
  expect_all_ids_character(unclass(post), "npde")
  expect_true(all(is.finite(post$sdtab$NPDE)))
})

test_that("every ID column is character: covariate model (#786)", {
  skip_on_cran()
  ex <- ferx_example("two_cpt_oral_cov")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE,
                                   settings = list(maxiter = 2L)))
  expect_all_ids_character(unclass(fit), "fit")
})

test_that("every ID column is character: nested levels (#786)", {
  skip_on_cran()
  ex <- ferx_example("nested_levels")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE,
                                   settings = list(maxiter = 0L)))
  expect_all_ids_character(unclass(fit), "fit")
})

test_that("every ID column is character: adaptive dosing (#786)", {
  skip_on_cran()
  ex <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_all_ids_character(unclass(res), "adaptive")
})

test_that("every ID column is character: survival fit and its TTE NPDE (#786)", {
  skip_on_cran()
  ex <- ferx_example("pktte_joint")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE))
  post <- suppressWarnings(ferx_calc_npde(fit, nsim = 50L, seed = 5L))
  expect_all_ids_character(unclass(post), "fit")
  tte <- ferx:::.fitrx_sdtab_blocks(post)$tte
  expect_gt(length(tte), 0L)
  expect_true(all(is.finite(post$sdtab$NPDE[tte])))
})

test_that("every ID column is character: SAEM conditional distributions (#786)", {
  skip_on_cran()
  fit <- warfarin_saem_conddist_fit()
  skip_if(is.null(fit$cond_dist), "no cond_dist result in this fit")
  expect_all_ids_character(unclass(fit), "fit")
  expect_identical(typeof(fit$cond_dist$data$ID), "character")
})

test_that("ferx_calc_npde(data = <frame>) names the ID pair that does not match (#786 review, 1)", {
  skip_on_cran()
  fit <- id_fit("padded")
  # Plain read.csv() reads "001" as the number 1, which is not the fit's text.
  err <- tryCatch(suppressWarnings(ferx_calc_npde(fit, data = utils::read.csv(fit$data_path),
                                                  nsim = 20L, seed = 1L)),
                  error = conditionMessage)
  expect_match(err, 'observation row 1 is subject "001" in fit$sdtab but "1"', fixed = TRUE)
  expect_match(err, 'colClasses = c(ID = "character")', fixed = TRUE)
  as_text <- utils::read.csv(fit$data_path, colClasses = c(ID = "character"))
  post <- suppressWarnings(ferx_calc_npde(fit, data = as_text, nsim = 20L, seed = 1L))
  expect_true(all(is.finite(post$sdtab$NPDE)))
})

test_that("ferx_calc_npde refuses NPDE rows whose text IDs do not line up (#786)", {
  attach_npde <- getFromNamespace(".ferx_attach_npde", "ferx")
  attach_tte <- getFromNamespace(".ferx_attach_tte_npde", "ferx")
  sdtab <- data.frame(ID = c("S1", "S2"), TIME = c(1, 1))
  swapped <- data.frame(ID = c("S2", "S1"), TIME = c(1, 1), NPDE = c(0.1, 0.2))
  # As numbers both IDs read NA, and NA "equals" NA: compared as text they differ.
  expect_error(attach_npde(sdtab, swapped),
               'observation row 1 is subject "S1" in fit$sdtab but "S2"', fixed = TRUE)
  expect_error(attach_tte(sdtab, swapped, 1:2),
               'time-to-event row 1 is subject "S1" in fit$sdtab but "S2"', fixed = TRUE)
  expect_identical(attach_npde(sdtab, swapped[2:1, ])$NPDE, c(0.2, 0.1))
  expect_identical(attach_tte(sdtab, swapped[2:1, ], 1:2)$NPDE, c(0.2, 0.1))
})
