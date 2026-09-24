test_that("ferx_inits_from_nca returns ferx_inits with expected shape", {
  ex <- ferx_example("warfarin")
  inits <- ferx_inits_from_nca(ex$model, ex$data)

  expect_s3_class(inits, "ferx_inits")
  expect_named(inits, c("theta", "omega", "method", "warnings"))
  expect_true(is.numeric(inits$theta))
  expect_false(is.null(names(inits$theta)))
  expect_true(length(inits$theta) > 0L)
  expect_true(is.matrix(inits$omega))
  expect_identical(inits$method, "nca_sweep")
})
test_that("ferx_inits_from_nca respects explicit method", {
  ex <- ferx_example("warfarin")

  inits_nca <- ferx_inits_from_nca(ex$model, ex$data, method = "nca")
  expect_identical(inits_nca$method, "nca")

  inits_ebe <- ferx_inits_from_nca(ex$model, ex$data, method = "nca_ebe")
  expect_true(inits_ebe$method %in% c("nca_ebe", "nca_sweep")) # may fall back for ODE
})
test_that("ferx_inits_from_nca rejects invalid method", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_inits_from_nca(ex$model, ex$data, method = "bogus"),
    "'arg'"
  )
})
test_that("print.ferx_inits emits theta and method", {
  ex <- ferx_example("warfarin")
  inits <- ferx_inits_from_nca(ex$model, ex$data, method = "nca")
  out <- capture.output(print(inits))
  expect_true(any(grepl("ferx_inits", out)))
  expect_true(any(grepl("Theta:", out)))
  expect_true(any(grepl("method = nca", out)))
})

# -- #391: the same reader as every other entry point --------------------------

# The bundled warfarin data with TIME renamed to TAFD, and the warfarin model
# with a `[data]` block mapping it back. `d` defaults to the unmodified rows.
column_mapped_warfarin <- function(d = example_rows(), env = parent.frame()) {
  names(d)[names(d) == "TIME"] <- "TAFD"
  # normalizePath(): on macOS tempdir() can end in "T//", and the model parser
  # reads `//` in a `[data]` path as the start of a comment.
  csv <- normalizePath(write_nonmem_csv(d, env))
  model <- model_with_lines(
    c("[data]", paste0("  path = ", csv), "  time = TAFD"), env
  )
  list(model = model, data = csv)
}

test_that("ferx_inits_from_nca honours the [data] column map (#391)", {
  ex <- ferx_example("warfarin")
  mapped <- column_mapped_warfarin()

  # The other entry points read this model; so must this one.
  expect_no_error(ferx_predict(mapped$model, mapped$data))

  for (method in c("nca", "nca_sweep")) {
    expect_identical(
      ferx_inits_from_nca(mapped$model, mapped$data, method = method),
      ferx_inits_from_nca(ex$model, ex$data, method = method),
      info = method
    )
  }
  # The data path from the [data] block, too.
  expect_identical(
    ferx_inits_from_nca(mapped$model, method = "nca"),
    ferx_inits_from_nca(ex$model, ex$data, method = "nca")
  )
})

test_that("a column-mapped read is not labelled with an unrelated validation code (#391)", {
  # An infusion into CMT = 0 on the column-mapped model. The old reader failed
  # first, on the mapped column, and the one validation finding was attached to
  # that unrelated read error. Now the read succeeds and the call stops at the
  # finding the code names.
  d <- example_rows()
  d$RATE[d$EVID == 1] <- 10
  d$CMT[d$EVID == 1]  <- 0
  mapped <- column_mapped_warfarin(d)

  probe <- engine_error_probe(ferx_inits_from_nca(mapped$model, mapped$data))
  expect_coded_refusal(probe, "infusion into compartment 0",
                       "E_DOSE_CMT_NOT_INFUSABLE")
  expect_false(grepl("Missing TIME column", conditionMessage(probe$cond),
                     fixed = TRUE))
})

test_that("ferx_inits_from_nca applies the model's [data_selection] (#391)", {
  sel_model <- model_with_lines(
    c("[data_selection]", "  ignore = TIME > 24"), environment()
  )
  d <- example_rows()
  kept <- d[d$TIME <= 24, ]
  expect_lt(nrow(kept), nrow(d))
  pre_filtered <- write_nonmem_csv(kept)
  ex <- ferx_example("warfarin")

  from_selection <- ferx_inits_from_nca(sel_model, ex$data, method = "nca")
  expect_identical(
    from_selection,
    ferx_inits_from_nca(ex$model, pre_filtered, method = "nca")
  )
  expect_false(identical(
    from_selection$theta,
    ferx_inits_from_nca(ex$model, ex$data, method = "nca")$theta
  ))
})
