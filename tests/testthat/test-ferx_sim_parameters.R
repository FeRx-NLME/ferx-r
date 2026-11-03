# Per-subject drawn parameters of a simulation (ferx-r #698, ferx-core #2505):
# `attr(sim, "parameters")`, read through `ferx_sim_parameters()`.

# One-compartment oral, log-normal eta on CL, V and KA. `cl` is CL's expression.
.sp_write_model <- function(cl = "TVCL * exp(ETA_CL)", extra = character()) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(2.0, 0.01, 50.0)",
    "  theta TVV(20.0, 0.5, 200.0)",
    "  theta TVKA(1.2, 0.005, 20.0)",
    "  omega ETA_CL ~ 0.09",
    "  omega ETA_V  ~ 0.04",
    "  omega ETA_KA ~ 0.16",
    "  sigma PROP ~ 0.1 (sd)",
    "[individual_parameters]",
    paste0("  CL = ", cl),
    "  V  = TVV * exp(ETA_V)",
    "  KA = TVKA * exp(ETA_KA)",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ proportional(PROP)",
    extra
  ), path)
  path
}

# `n` subjects, 100 mg at t = 0, observations at `times`; `wt(id, t)` is WT.
.sp_write_data <- function(n = 3L, times = c(0.5, 1, 2, 4, 8, 12),
                           wt = function(id, t) 70) {
  rows <- list()
  for (id in seq_len(n)) {
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, DV = 0, AMT = 100, EVID = 1, CMT = 1, MDV = 1,
      WT = wt(id, 0)
    )
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = times, DV = 1, AMT = 0, EVID = 0, CMT = 1, MDV = 0,
      WT = vapply(times, function(t) wt(id, t), numeric(1))
    )
  }
  path <- tempfile(fileext = ".csv")
  utils::write.csv(do.call(rbind, rows), path, row.names = FALSE, quote = FALSE)
  path
}

test_that("ferx_simulate() carries one parameter row per (DRAW, SIM, ID) that joins back", {
  sim <- ferx_simulate(.sp_write_model(), .sp_write_data(n = 4L),
                       n_sim = 3L, seed = 2505L)
  params <- ferx_sim_parameters(sim)

  expect_s3_class(params, "ferx_sim_parameters")
  expect_s3_class(params, "data.frame")
  expect_identical(names(params),
                   c("DRAW", "SIM", "ID", "ETA_CL", "ETA_V", "ETA_KA",
                     "CL", "V", "KA"))
  expect_type(params$DRAW, "integer")
  expect_type(params$SIM, "integer")
  expect_type(params$ID, "character")
  expect_identical(attr(params, "varying"), character(0))

  # One row per key, the same keys the observation rows carry, in their order.
  keys <- unique(sim[c("DRAW", "SIM", "ID")])
  expect_equal(nrow(params), 3L * 4L)
  expect_equal(params[c("DRAW", "SIM", "ID")], keys, ignore_attr = TRUE)

  # The join leaves the observation rows as they were.
  joined <- merge(sim, params, by = c("DRAW", "SIM", "ID"))
  expect_equal(nrow(joined), nrow(sim))

  # The individual parameters are the drawn eta folded into theta.
  expect_equal(params$CL, 2.0 * exp(params$ETA_CL))
  expect_equal(params$V, 20.0 * exp(params$ETA_V))
  expect_equal(params$KA, 1.2 * exp(params$ETA_KA))
  expect_true(all(params$ETA_CL != 0))
})

test_that("F * dose / CL from the table is the steady-state AUC over tau (ferx-core #2522)", {
  # 100 mg every 12 h at steady state, IPRED on a 0.01 h grid over one interval:
  # for a linear one-compartment oral model AUC_tau is exactly F * dose / CL
  # (F = 1). The trapezoid of the simulated IPRED must reproduce it from the
  # returned CL, as in ferx-core's own check.
  tau <- 12
  dose <- 100
  grid <- seq(0, tau, by = 0.01)
  rows <- list()
  for (id in 1:5) {
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, DV = 0, AMT = dose, EVID = 1, CMT = 1, MDV = 1,
      SS = 1, II = tau
    )
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = grid, DV = 1, AMT = 0, EVID = 0, CMT = 1, MDV = 0,
      SS = 0, II = 0
    )
  }
  data <- tempfile(fileext = ".csv")
  utils::write.csv(do.call(rbind, rows), data, row.names = FALSE, quote = FALSE)

  sim <- ferx_simulate(.sp_write_model(), data, n_sim = 2L, seed = 7L)
  params <- ferx_sim_parameters(sim)
  expect_equal(nrow(params), 2L * 5L)

  trapezoid <- function(t, y) sum(diff(t) * (head(y, -1L) + tail(y, -1L)) / 2)
  rel <- vapply(seq_len(nrow(params)), function(i) {
    p <- params[i, ]
    rows <- sim[sim$DRAW == p$DRAW & sim$SIM == p$SIM & sim$ID == p$ID, ]
    rows <- rows[order(rows$TIME), ]
    auc <- trapezoid(rows$TIME, rows$IPRED)
    want <- dose / p$CL
    abs(auc - want) / want
  }, numeric(1))
  # ferx-core measured 1.99e-6 (the trapezoid's O(dt^2) error) and allows 1e-5;
  # another subject's CL, or the typical one, is tens of percent off.
  expect_lt(max(rel), 1e-5)
})

test_that("[derived] and [output] columns come through by name", {
  ex <- ferx_example("warfarin_derived")
  sim <- suppressWarnings(ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 3L))
  params <- ferx_sim_parameters(sim)

  for (col in c("KE", "T_HALF", "CMAX", "AUC_0_72")) {
    expect_true(col %in% names(params), info = col)
    expect_type(params[[col]], "double")
  }
  expect_true(all(c("ETA_CL", "ETA_V", "ETA_KA", "CL", "V", "KA") %in% names(params)))
  expect_equal(params$KE, params$CL / params$V)
  expect_equal(params$T_HALF, 0.6931472 / params$KE)

  # CMAX is the simulated subject's own IPRED maximum, not the typical one.
  key <- paste(sim$DRAW, sim$SIM, sim$ID)
  cmax <- tapply(sim$IPRED, key, max)
  expect_equal(as.vector(cmax[paste(params$DRAW, params$SIM, params$ID)]),
               params$CMAX)

  # Per-row derived quantities differ within a subject and are flagged.
  expect_true(all(c("DAY", "TAU_TIME") %in% attr(params, "varying")))
  expect_false("KE" %in% attr(params, "varying"))
})

test_that("a time-varying covariate lists the affected parameter in `varying`", {
  model <- .sp_write_model("TVCL * (WT / 70) * exp(ETA_CL)")
  # WT 70 at the dose and first observation, 90 from t = 4 on.
  tv <- ferx_simulate(
    model,
    .sp_write_data(n = 2L, wt = function(id, t) if (t < 4) 70 else 90),
    n_sim = 1L, seed = 3L
  )
  params <- ferx_sim_parameters(tv)
  expect_identical(attr(params, "varying"), "CL")
  # The reported CL is the first observation row's, at WT = 70.
  expect_equal(params$CL, 2.0 * exp(params$ETA_CL))
  # Printing says which values are first-row snapshots.
  expect_output(print(params), "Varying within subject.*CL")

  flat <- ferx_simulate(model, .sp_write_data(n = 2L), n_sim = 1L, seed = 3L)
  flat_params <- ferx_sim_parameters(flat)
  expect_identical(attr(flat_params, "varying"), character(0))
  expect_false(any(grepl("Varying within subject",
                         capture.output(print(flat_params)))))
})

test_that("a [mixture] model's table carries MIXNUM, matching the observation rows", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL1(1.0, 0.01, 100.0)",
    "  theta TVCL2(10.0, 0.01, 100.0)",
    "  theta TVV(10.0, 0.1, 1000.0)",
    "  theta MIXL(-0.8472979, -10.0, 10.0)",
    "  omega ETA_CL ~ 0.01",
    "  sigma EPS ~ 0.01",
    "[mixture]",
    "  nsub = 2",
    "  logit(1) = MIXL",
    "[individual_parameters]",
    "  CL = if (MIXNUM == 1) TVCL1 * exp(ETA_CL) else TVCL2 * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  pk one_cpt_iv(cl=CL, v=V)",
    "[error_model]",
    "  DV ~ proportional(EPS)"
  ), path)
  sim <- muffle_mixl_lint(ferx_simulate(path,
                                        test_path("fixtures", "mixture_iv.csv"),
                                        n_sim = 5L, seed = 7L))
  params <- ferx_sim_parameters(sim)

  expect_identical(names(params)[1:4], c("DRAW", "SIM", "ID", "MIXNUM"))
  expect_type(params$MIXNUM, "integer")
  expect_true(all(c(1L, 2L) %in% params$MIXNUM))
  per_run <- unique(sim[c("DRAW", "SIM", "ID", "MIXNUM")])
  expect_equal(params[c("DRAW", "SIM", "ID", "MIXNUM")], per_run,
               ignore_attr = TRUE)
  # CL is the drawn class's.
  scale <- ifelse(params$MIXNUM == 1L, 1.0, 10.0)
  expect_equal(params$CL, scale * exp(params$ETA_CL))

  # A model without [mixture] has no MIXNUM column.
  plain <- ferx_simulate(.sp_write_model(), .sp_write_data(), n_sim = 1L)
  expect_false("MIXNUM" %in% names(ferx_sim_parameters(plain)))
})

test_that("an IOV model's table carries one kappa column per occasion", {
  ex <- ferx_example("warfarin_iov")
  sim <- suppressWarnings(ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 4L))
  params <- ferx_sim_parameters(sim)

  kappa_cols <- grep("^KAPPA_CL_[0-9]+$", names(params), value = TRUE)
  expect_gte(length(kappa_cols), 2L)
  expect_identical(kappa_cols, paste0("KAPPA_CL_", seq_along(kappa_cols)))
  expect_type(params$KAPPA_CL_1, "double")
  expect_true(all(is.finite(params$KAPPA_CL_1)))
  # The kappa is drawn: not all zero, and not shared between occasions.
  expect_true(any(params$KAPPA_CL_1 != 0))
  expect_false(isTRUE(all.equal(params$KAPPA_CL_1, params$KAPPA_CL_2)))
  # CL carries the first observation row's occasion kappa, and IOV makes CL
  # differ within a subject, so it is flagged.
  expect_true("CL" %in% attr(params, "varying"))

  # A model without kappa has no kappa columns.
  plain <- ferx_simulate(.sp_write_model(), .sp_write_data(), n_sim = 1L)
  expect_false(any(grepl("^KAPPA", names(ferx_sim_parameters(plain)))))
})

test_that("ferx_simulate_with_uncertainty() carries the table keyed by DRAW", {
  fit <- suppressWarnings(warfarin_fit_cov())
  ex <- ferx_example("warfarin")
  sim <- suppressWarnings(ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit, n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    seed = 9L
  ))
  params <- ferx_sim_parameters(sim)
  n_id <- length(unique(sim$ID))

  expect_equal(nrow(params), 3L * 2L * n_id)
  expect_setequal(unique(params$DRAW), 1:3)
  keys <- unique(sim[c("DRAW", "SIM", "ID")])
  expect_equal(params[c("DRAW", "SIM", "ID")], keys, ignore_attr = TRUE)
  joined <- merge(sim, params, by = c("DRAW", "SIM", "ID"))
  expect_equal(nrow(joined), nrow(sim))
  # Each DRAW has its own theta: a subject's CL / exp(ETA_CL) (its TVCL) is
  # shared within a draw and differs between draws.
  tvcl <- tapply(params$CL / exp(params$ETA_CL), params$DRAW,
                 function(x) diff(range(x)) / mean(x))
  expect_true(all(tvcl < 1e-10))
  per_draw <- tapply(params$CL / exp(params$ETA_CL), params$DRAW, mean)
  expect_gt(diff(range(per_draw)), 0)
})

test_that("ferx_sim_parameters() refuses an object without the table", {
  expect_error(ferx_sim_parameters(data.frame(ID = 1)),
               "carries no drawn-parameter table")
})
