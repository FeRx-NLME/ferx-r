# ferx_mbma_data(): build and validate a summary-level dataset (#370, #298).
#
# Each check in the helper catches a data-preparation mistake that otherwise
# fits without complaint (#298's table), so each one has a test here that
# fails when the check is removed. The engine tests at the end are the oracle
# for the level index: the counted form indexed by LEVEL_IDX must fit exactly
# like the data-driven form, and the helper's level table must name the same
# levels, in the same order, as the fit's theta_levels.

# Two studies x two arms x two weeks. Study "B" is listed first and sorts
# second, and "placebo" sorts after "drug", so the output order is checked.
mb_arms <- function() {
  data.frame(
    trial = rep(c("B", "A"), each = 4),
    treatment = rep(rep(c("placebo", "drug"), each = 2), 2),
    week = rep(c(0, 12), 4),
    change = c(0, -1.1, 0, -2.4, 0, -0.8, 0, -2.0),
    se = c(0.20, 0.21, 0.19, 0.22, 0.18, 0.23, 0.20, 0.24),
    sd = c(2.1, 2.3, 2.0, 2.4, 1.9, 2.2, 2.1, 2.5),
    n = c(120, 120, 118, 118, 95, 95, 97, 97),
    region = rep(c("US", "EU"), each = 4),
    stringsAsFactors = FALSE
  )
}

mb <- function(data = mb_arms(), ...) {
  args <- list(
    data = data, study = "trial", arm = "treatment", time = "week",
    mean = "change", n = "n", se = "se"
  )
  extra <- list(...)
  if (!is.null(extra$sd)) args$se <- NULL
  # A NULL in `...` drops the argument, so `se = NULL` passes neither.
  do.call(ferx_mbma_data, utils::modifyList(args, extra))
}

# --- the frame ----------------------------------------------------------------

test_that("the frame is sorted by study, arm and time with ID = STUDY", {
  d <- mb()
  expect_s3_class(d, "ferx_mbma_data")
  expect_identical(
    names(d), c("ID", "STUDY", "OCC", "TIME", "DV", "MDV", "NARM", "SE")
  )
  # A -> 1, B -> 2; within each, drug -> OCC 1, placebo -> OCC 2.
  expect_identical(d$STUDY, rep(1:2, each = 4))
  expect_identical(d$ID, d$STUDY)
  expect_identical(d$OCC, rep(rep(1:2, each = 2), 2))
  expect_identical(d$TIME, rep(c(0, 12), 4))
  expect_identical(d$DV, c(0, -2.0, 0, -0.8, 0, -2.4, 0, -1.1))
  expect_identical(d$SE, c(0.20, 0.24, 0.18, 0.23, 0.19, 0.22, 0.20, 0.21))
  expect_identical(d$NARM, c(97, 97, 95, 95, 118, 118, 120, 120))
  expect_identical(d$MDV, rep(0L, 8))
  expect_identical(attr(d, "codes"), list(study = c(A = 1L, B = 2L)))
  expect_identical(attr(d, "arms"), data.frame(
    STUDY = c(1L, 1L, 2L, 2L), OCC = c(1L, 2L, 1L, 2L),
    study = c("A", "A", "B", "B"), arm = c("drug", "placebo", "drug", "placebo")
  ))
  expect_null(attr(d, "levels"))
  expect_identical(attr(d, "scale"), "continuous")
})

test_that("a numeric study keeps its own values and sorts numerically", {
  x <- mb_arms()
  x$trial <- rep(c(10, 2), each = 4)
  d <- mb(x)
  expect_identical(d$STUDY, rep(c(2, 10), each = 4))
  expect_null(attr(d, "codes")$study)
})

test_that("covariates are carried through; character coded, logical 0/1", {
  x <- mb_arms()
  x$flare <- rep(c(TRUE, FALSE), each = 4)
  x$dose <- rep(c(0, 50), 4)
  d <- mb(x, covariates = c("region", "flare", "dose"))
  expect_identical(d$region, rep(c(1L, 2L), each = 4))
  expect_identical(d$flare, rep(c(0L, 1L), each = 4))
  expect_identical(d$dose, c(0, 50, 0, 50, 0, 50, 0, 50))
  expect_identical(attr(d, "codes")$region, c(EU = 1L, US = 2L))
  expect_null(attr(d, "codes")$flare)
})

test_that("a covariate or index named like an output column is refused", {
  x <- mb_arms()
  x$SE <- 1
  expect_error(
    mb(x, covariates = "SE"),
    "Covariate column(s) `SE` would overwrite an output column; rename them first.",
    fixed = TRUE
  )
  x$TIME <- 1L
  expect_error(
    mb(x, index = "TIME"),
    "`index` column `TIME` would overwrite an output column; rename it first.",
    fixed = TRUE
  )
})

test_that("the arguments themselves are checked", {
  expect_error(mb(list()), "`data` must be a data.frame.", fixed = TRUE)
  expect_error(mb(study = c("a", "b")), "`study` must be a single column name.", fixed = TRUE)
  expect_error(mb(covariates = 1), "`covariates` must be a character vector", fixed = TRUE)
  expect_error(
    mb(n = "N", covariates = "WT"),
    "`data` has no column(s) `N`, `WT`.",
    fixed = TRUE
  )
  x <- mb_arms()
  x$week <- as.character(x$week)
  expect_error(mb(x), "Column `week` (`time`) must be numeric.", fixed = TRUE)
  x <- mb_arms()
  x$treatment[3] <- NA
  expect_error(
    mb(x), "`study`, `arm` and `time` must not be missing; row(s) 3.",
    fixed = TRUE
  )
})

test_that("rows with a missing mean are dropped, with a message", {
  x <- mb_arms()
  x$change[2] <- NA
  x$se[2] <- NA
  expect_message(d <- mb(x), "Dropped 1 row(s) with a missing `mean`.", fixed = TRUE)
  expect_identical(nrow(d), 7L)
})

test_that("a negative mean is accepted (change from baseline)", {
  d <- expect_no_warning(mb())
  expect_true(any(d$DV < 0))
})

# --- row 1: se ---------------------------------------------------------------

test_that("a negative, zero or missing se is refused, naming the cell", {
  for (bad in c(-0.1, 0, NA)) {
    x <- mb_arms()
    x$se[6] <- bad
    expect_error(
      mb(x),
      "`se` must be positive wherever `mean` is present: (study A, arm placebo, time 12).",
      fixed = TRUE
    )
  }
})

test_that("a non-positive n is refused, naming the cell", {
  x <- mb_arms()
  x$n[1] <- 0
  expect_error(
    mb(x),
    "`n` must be positive and present on every row: (study B, arm placebo, time 0).",
    fixed = TRUE
  )
})

# --- row 4: se or sd ----------------------------------------------------------

test_that("exactly one of se and sd", {
  expect_error(mb(sd = "sd", se = "se"), "Give exactly one of `se` and `sd`.", fixed = TRUE)
  expect_error(mb(se = NULL), "Give exactly one of `se` and `sd`.", fixed = TRUE)
})

test_that("sd is converted to se = sd / sqrt(n), and says so", {
  expect_message(
    d <- mb(sd = "sd"),
    "Converted `sd` to a standard error as sd / sqrt(n).",
    fixed = TRUE
  )
  x <- mb_arms()
  ord <- order(x$trial, x$treatment, x$week)
  expect_equal(d$SE, (x$sd / sqrt(x$n))[ord], tolerance = 1e-15)
  expect_no_message(mb())
})

test_that("a negative sd is refused under its own name", {
  x <- mb_arms()
  x$sd[1] <- -2
  expect_error(
    mb(x, sd = "sd"),
    "`sd` must be positive wherever `mean` is present: (study B, arm placebo, time 0).",
    fixed = TRUE
  )
})

# --- rows 2 and 3: scale -----------------------------------------------------

test_that("a proportion outside [0, 1] is refused; the percent hint is gated", {
  x <- mb_arms()
  x$change <- c(35, 30, 36, 28, 40, 33, 41, 30)
  e <- expect_error(mb(x, scale = "proportion"))
  expect_match(
    conditionMessage(e),
    "`mean` is outside [0, 1] for scale = \"proportion\": (study B, arm placebo, time 0),",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(e),
    "The values look like percentages; use scale = \"percent\".",
    fixed = TRUE
  )
  # Above 100 it cannot be a percentage: no hint.
  x$change[1] <- 135
  e <- expect_error(mb(x, scale = "proportion"))
  expect_no_match(conditionMessage(e), "percent", fixed = TRUE)
  # In range, accepted.
  x$change <- x$change / 200
  expect_no_error(mb(x, scale = "proportion"))
})

test_that("a percentage outside [0, 100] is refused, and in range is converted", {
  x <- mb_arms()
  x$change <- c(35, 30, 36, 28, 40, 33, 41, 130)
  expect_error(
    mb(x, scale = "percent"),
    "`mean` is outside [0, 100] for scale = \"percent\": (study A, arm drug, time 12).",
    fixed = TRUE
  )
  x$change[8] <- 30
  expect_message(
    d <- mb(x, scale = "percent"),
    "Converted `mean` and `se` from percent to proportions (divided by 100).",
    fixed = TRUE
  )
  expect_equal(d$DV, c(0.41, 0.30, 0.40, 0.33, 0.36, 0.28, 0.35, 0.30), tolerance = 1e-15)
  expect_equal(d$SE, c(0.20, 0.24, 0.18, 0.23, 0.19, 0.22, 0.20, 0.21) / 100, tolerance = 1e-15)
  expect_identical(attr(d, "scale"), "percent")
})

# --- rows 5, 6, 7: warnings ----------------------------------------------------

test_that("n changing within an arm warns, naming the arm, and NARM keeps each row", {
  x <- mb_arms()
  x$n[2] <- 110
  expect_warning(
    d <- mb(x),
    "`n` changes between timepoints within study B arm placebo; each row keeps its own n in NARM.",
    fixed = TRUE
  )
  expect_identical(d$NARM[7:8], c(120, 110))
})

test_that("a single-arm study warns, naming the study", {
  x <- mb_arms()[-(3:4), ]
  expect_warning(
    mb(x),
    "Study(ies) B contribute a single arm, so no comparison within them informs a treatment effect.",
    fixed = TRUE
  )
})

test_that("an arm without a time 0 record warns without calling it an error", {
  x <- mb_arms()[-5, ]
  w <- expect_warning(mb(x))
  expect_identical(
    conditionMessage(w),
    paste(
      "No time = 0 record for study A arm placebo.",
      "Expected for change-from-baseline data; otherwise check the time axis."
    )
  )
})

# --- row 8: duplicates ---------------------------------------------------------

test_that("a repeated (study, arm, time) is refused, naming it", {
  x <- mb_arms()
  x <- rbind(x, x[2, ])
  expect_error(
    mb(x),
    "Each (study, arm, time) must appear once; repeated: (study B, arm placebo, time 12).",
    fixed = TRUE
  )
})

# --- row 9: a user index ------------------------------------------------------

test_that("a user index is checked and carried through", {
  x <- mb_arms()
  x$IDX <- rep(c(1, 2), 4)
  d <- mb(x, index = "IDX")
  expect_identical(d$IDX, rep(1:2, 4))
  expect_false("LEVEL_IDX" %in% names(d))

  for (bad in list(c(NA, 2), c(1.5, 2), c(0, 2))) {
    x$IDX <- rep(bad, 4)
    expect_error(
      mb(x, index = "IDX"),
      sprintf(
        "`index` column `IDX` must hold whole numbers from 1; offending value(s) %s at (study B, arm placebo, time 0),",
        bad[1]
      ),
      fixed = TRUE
    )
  }
  x$IDX <- rep(c(1, 4), 4)
  expect_error(
    mb(x, index = "IDX"),
    "`index` column `IDX` skips level(s) 2, 3 of 1..4; a counted block's levels must all appear.",
    fixed = TRUE
  )
  x$IDX <- as.character(rep(1:2, 4))
  expect_error(mb(x, index = "IDX"), "`index` column `IDX` must be numeric.", fixed = TRUE)
  expect_error(
    mb(x, index = "IDX", levels = "STUDY"),
    "Give `levels` or `index`, not both.",
    fixed = TRUE
  )
})

# --- levels -------------------------------------------------------------------

test_that("levels numbers the observed combinations in ascending order", {
  x <- mb_arms()
  x$trial <- rep(c(10, 2), each = 4)
  x$week <- rep(c(0, 0.5), 4)
  d <- mb(x, levels = c("STUDY", "TIME"))
  # STUDY 2 sorts before 10 numerically; both arms of a study share a level.
  expect_identical(d$LEVEL_IDX, c(1L, 2L, 1L, 2L, 3L, 4L, 3L, 4L))
  expect_identical(attr(d, "levels"), data.frame(
    index = 1:4,
    label = c("STUDY=2,TIME=0", "STUDY=2,TIME=0.5", "STUDY=10,TIME=0", "STUDY=10,TIME=0.5"),
    STUDY = c(2, 2, 10, 10),
    TIME = c(0, 0.5, 0, 0.5)
  ))
})

test_that("levels are sorted, not numbered in order of appearance", {
  # Drug (OCC 1) is seen at weeks 0 and 12; placebo (OCC 2) adds week 4,
  # which first appears after week 12 in the sorted frame.
  x <- mb_arms()[mb_arms()$trial == "A", ]
  x <- rbind(x, transform(x[1, ], week = 4, change = -0.5))
  d <- mb(x, levels = c("STUDY", "TIME"))
  expect_identical(attr(d, "levels")$label, c("STUDY=1,TIME=0", "STUDY=1,TIME=4", "STUDY=1,TIME=12"))
  expect_identical(d$TIME, c(0, 12, 0, 4, 12))
  expect_identical(d$LEVEL_IDX, c(1L, 3L, 1L, 2L, 3L))
})

test_that("levels must name output columns", {
  expect_error(
    mb(levels = c("study", "TIME")),
    "`levels` must name distinct output columns among `STUDY`, `TIME`.",
    fixed = TRUE
  )
  expect_error(
    mb(levels = c("TIME", "TIME")),
    "`levels` must name distinct output columns",
    fixed = TRUE
  )
  d <- mb(levels = c("region", "TIME"), covariates = "region")
  expect_identical(attr(d, "levels")$label[1], "region=1,TIME=0")
})

test_that("a level value is labelled in its shortest form", {
  lv <- ferx:::.mbma_level_value
  expect_identical(lv(7), "7")
  expect_identical(lv(-3), "-3")
  expect_identical(lv(0.5), "0.5")
  expect_identical(lv(2.25), "2.25")
  expect_identical(lv(0.1), "0.1")
})

# --- print -------------------------------------------------------------------

test_that("print reports studies, arms, timepoints, implied SD and levels", {
  x <- mb_arms()[-8, ]
  d <- suppressWarnings(mb(x, levels = c("STUDY", "TIME")))
  out <- capture.output(print(d, n = 2))
  expect_identical(out[1:5], c(
    "<ferx_mbma_data>  7 records, 2 studies, 4 arms",
    "  arms per study:      2",
    "  timepoints per arm:  1-2",
    "  implied SD (SE * sqrt(NARM)): 1.75 to 2.39",
    "  LEVEL_IDX: 4 levels of (STUDY, TIME); declare theta NAME[4]"
  ))
  expect_identical(out[length(out)], "... 5 more records")
  expect_false(any(grepl("LEVEL_IDX:", capture.output(print(mb())), fixed = TRUE)))
  # Dropping the column keeps the attribute; the line goes with the column.
  d$LEVEL_IDX <- NULL
  expect_false(any(grepl("LEVEL_IDX:", capture.output(print(d)), fixed = TRUE)))
})

# --- the engine reads it ------------------------------------------------------

mb_model <- function(theta, effect) {
  path <- tempfile(fileext = ".ferx")
  writeLines(paste0(
    "[parameters]\n", theta, "\n",
    "  sigma ADD_ERR ~ 1.0 (variance) FIX\n\n",
    "[covariates]\n",
    "  SE continuous\n\n",
    "[individual_parameters]\n",
    "  E = ", effect, "\n\n",
    "[structural_model]\n",
    "  y = E\n\n",
    "[error_model]\n",
    "  DV ~ additive(ADD_ERR) weight = SE\n\n",
    "[fit_options]\n",
    "  method = focei\n",
    "  maxiter = 20\n",
    "  covariance = false\n"
  ), path)
  path
}

test_that("the counted form on LEVEL_IDX fits like the data-driven form", {
  x <- mb_arms()
  x$trial <- rep(c(10, 2), each = 4)
  x$week <- rep(c(0, 0.5), 4)
  d <- mb(x, levels = c("STUDY", "TIME"))
  data <- tempfile(fileext = ".csv")
  utils::write.csv(d, data, row.names = FALSE)

  col <- mb_model("  theta PLACEBO[STUDY, TIME, contrast = none](0.0, -10.0, 10.0)", "PLACEBO")
  counted <- mb_model("  theta PLACEBO[4](0.0, -10.0, 10.0)", "PLACEBO[LEVEL_IDX]")
  # gradient = "fd": on the default gradient a gathered theta never leaves
  # its initial value (FeRx-NLME/ferx-core#1628), and two fits stuck at
  # their inits agree whatever the index says.
  fit <- function(model) {
    suppressWarnings(ferx_fit(model, data, verbose = FALSE, gradient = "fd"))
  }
  f_col <- fit(col)
  f_cnt <- fit(counted)

  # The helper's table names the engine's levels, in the engine's order.
  expect_identical(attr(d, "levels")$label, f_col$theta_levels$label)
  expect_identical(
    names(f_col$theta), paste0("PLACEBO[", attr(d, "levels")$label, "]")
  )
  # Closed form, outside the engine: with sigma fixed to 1 and weight = SE,
  # each level is its cell's inverse-variance weighted mean (0 at TIME 0,
  # -1.3745 and -1.7198 at TIME 0.5). Measured worst error 1.3e-15.
  w <- 1 / d$SE^2
  cell_mean <- unname(tapply(d$DV * w, d$LEVEL_IDX, sum) / tapply(w, d$LEVEL_IDX, sum))
  expect_lt(max(abs(unname(f_cnt$theta) - cell_mean)), 1e-8)
  expect_equal(f_col$ofv, f_cnt$ofv, tolerance = 1e-10)
  expect_equal(unname(f_col$theta), unname(f_cnt$theta), tolerance = 1e-10)
})
