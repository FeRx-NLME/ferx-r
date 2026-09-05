# MFL parsing, `@`-symbol resolution against a base model, and the coverage
# table. The grammar is the engine's; what is asserted here is that R reports
# the expansion the engine produced, and that a gap is a row rather than an
# aborted run.

test_that("a space parses without a model, symbols left as written", {
  sp <- ferx_search_space("COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])")

  expect_s3_class(sp, "ferx_search_space")
  expect_s3_class(sp, "data.frame")
  expect_equal(names(sp), c("feature", "keyword", "optional"))
  expect_equal(nrow(sp), 1L)
  expect_equal(sp$keyword, "COVARIATE")
  expect_true(sp$optional)
  expect_match(sp$feature, "@IIV")
  expect_false(attr(sp, "resolved"))
  expect_equal(nrow(attr(sp, "covariate_effects")), 0L)
})

test_that("several statements come back as several rows", {
  sp <- ferx_search_space("ABSORPTION([INST, FO]); ELIMINATION(FO)")
  expect_equal(nrow(sp), 2L)
  expect_equal(sp$keyword, c("ABSORPTION", "ELIMINATION"))
  expect_false(any(sp$optional))
})

test_that("resolving against a model expands the symbols to ground terms", {
  ex <- ferx_example("two_cpt_oral_cov")
  sp <- ferx_search_space("COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])",
                          model = ex$model, data = ex$data)

  expect_true(attr(sp, "resolved"))
  # The unresolved space is kept alongside the expansion.
  expect_match(attr(sp, "parsed")$feature, "@IIV")
  expect_false(any(grepl("@", sp$feature)))

  effects <- attr(sp, "covariate_effects")
  # 5 parameters with an eta (CL, V1, Q, V2, KA) x 2 declared continuous
  # covariates (WT, CRCL) x 2 effect forms.
  expect_equal(nrow(effects), 20L)
  expect_setequal(unique(effects$parameter), c("CL", "V1", "Q", "V2", "KA"))
  expect_setequal(unique(effects$covariate), c("WT", "CRCL"))
  expect_setequal(unique(effects$effect), c("pow", "lin"))
  expect_true(all(effects$optional))
  expect_true(all(effects$op == "*"))
})

test_that("a configuration resolves against its own base model", {
  ex  <- ferx_example("two_cpt_oral_cov")
  cfg <- ferx_search_config(ex$search)

  from_cfg <- ferx_search_space(cfg)
  from_str <- ferx_search_space(cfg$mfl, model = cfg$base, data = cfg$data)

  # The two entry forms describe the same space.
  expect_equal(as.data.frame(from_cfg), as.data.frame(from_str))
  expect_equal(attr(from_cfg, "covariate_effects"),
               attr(from_str, "covariate_effects"))
})

test_that("a symbol resolving to nothing drops its feature with a note", {
  ex <- ferx_example("two_cpt_oral_cov")
  # The dataset declares no categorical covariates, so this space is empty on
  # this model - Pharmpy's behaviour, and a note rather than an error.
  sp <- ferx_search_space("COVARIATE?(@IIV, @CATEGORICAL, cat)",
                          model = ex$model, data = ex$data)
  expect_equal(nrow(attr(sp, "covariate_effects")), 0L)
  expect_true(length(attr(sp, "notes")) > 0L)
  expect_match(paste(attr(sp, "notes"), collapse = " "), "dropped")
})

test_that("a name that is not a parameter is an error naming it", {
  ex <- ferx_example("two_cpt_oral_cov")
  expect_error(
    ferx_search_space("COVARIATE?(NOPE, WT, pow)",
                      model = ex$model, data = ex$data),
    "NOPE"
  )
})

test_that("caller errors are caught before the engine is called", {
  expect_error(ferx_search_space("COVARIATE?(CL, WT, pow)",
                                 model = tempfile(fileext = ".ferx")),
               "Model file not found")
  expect_error(ferx_search_space("COVARIATE?(CL, WT, pow)",
                                 data = "nowhere.csv"),
               "needs a 'model'")
  expect_error(ferx_search_space(42), "single MFL string")
})

test_that("an unparseable space is an error", {
  expect_error(ferx_search_space("COVARIATE?(CL, WT"))
})

test_that("the space prints its source and its features", {
  sp <- ferx_search_space("ABSORPTION([INST, FO])")
  expect_output(print(sp), "ferx search space")
  expect_output(print(sp), "ABSORPTION")
})

# -- coverage ---------------------------------------------------------------

test_that("a covered space is all TRUE with no reasons", {
  cov <- ferx_search_coverage("ABSORPTION([INST, FO]); ELIMINATION(FO)")
  expect_equal(names(cov), c("feature", "covered", "reason"))
  expect_true(all(cov$covered))
  expect_true(all(is.na(cov$reason)))
})

test_that("an unsupported feature is a row, not an aborted run", {
  cov <- ferx_search_coverage("ELIMINATION([FO, MM])")
  expect_true(any(!cov$covered))
  gap <- cov[!cov$covered, ]
  expect_match(gap$feature, "MM")
  expect_true(all(nzchar(gap$reason)))
})

test_that("coverage accepts a space or a configuration", {
  sp <- ferx_search_space("ELIMINATION([FO, MM])")
  expect_equal(ferx_search_coverage(sp),
               ferx_search_coverage("ELIMINATION([FO, MM])"))

  ex  <- ferx_example("two_cpt_oral_cov")
  cfg <- ferx_search_config(ex$search)
  expect_true(all(ferx_search_coverage(cfg)$covered))

  expect_error(ferx_search_coverage(42), "single MFL string")
})
