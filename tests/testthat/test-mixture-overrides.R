# A [mixture] fit's per-class overrides (`omega(k)` / `sigma(k)`) with their
# estimate, standard error and name (ferx-core #2833): the fit fields, the rows
# of `fit$estimates` and the section print() adds. The fit is
# fixtures/mixture_override.ferx, whose class 2 overrides ETA_CL.

quiet <- function(x) suppressWarnings(suppressMessages(x))

label <- "ETA_CL (class 2)"

# The print() rows of the override section: from its header to the next blank
# line.
override_section <- function(fit) {
  out <- utils::capture.output(quiet(print(fit)))
  at <- grep("MIXTURE  (per-class overrides)", out, fixed = TRUE)
  if (length(at) == 0L) return(character())
  expect_length(at, 1L)
  rest <- out[seq(at + 2L, length(out))]
  rest[seq_len(which(c(rest, "") == "")[1L] - 1L)]
}

test_that("an override is named, with its SE, on the fit (ferx-core #2833)", {
  skip_on_cran()
  fit <- mixture_override_fit()
  expect_false(is.null(fit$cov_matrix))
  expect_identical(
    fit$mixture_override_params,
    data.frame(parameter = label, kind = "omega", index = 1L, class = 2L,
               fixed = FALSE, stringsAsFactors = FALSE)
  )
  expect_named(fit$mixture_overrides, label)
  expect_named(fit$se_mixture_overrides, label)
  # The override is its own coordinate of the covariance step, under the same
  # label; its SE is not its base's.
  expect_true(label %in% rownames(fit$cov_matrix))
  se <- fit$se_mixture_overrides[[label]]
  expect_true(is.finite(se) && se > 0)
  # The coordinate packs the log of the class SD, so by the delta method the
  # variance's SE is 2 * variance * the packed SE.
  expect_equal(se, 2 * fit$mixture_overrides[[label]] * sqrt(fit$cov_matrix[label, label]),
               tolerance = 1e-8)
  base_se <- fit$estimates$se[fit$estimates$param == "ETA_CL"]
  expect_false(isTRUE(all.equal(se, base_se)))
  expect_identical(summary(fit)$se_mixture_overrides, fit$se_mixture_overrides)
  expect_identical(summary(fit)$mixture_overrides, fit$mixture_overrides)
})

test_that("fit$estimates has one row per override, last (ferx-core #2833)", {
  skip_on_cran()
  fit <- mixture_override_fit()
  est <- fit$estimates
  row <- est[nrow(est), ]
  expect_identical(rownames(row), label)
  expect_identical(row$param, label)
  expect_identical(row$class, 2L)
  expect_identical(row$transform, "variance")
  expect_false(row$fixed)
  expect_identical(row$estimate, fit$mixture_overrides[[label]])
  expect_identical(row$se, fit$se_mixture_overrides[[label]])
  expect_equal(row$rse_pct, 100 * row$se / row$estimate)
  # On its base ETA's scale, and no other row has a class.
  expect_identical(row$scale, est$scale[est$param == "ETA_CL"])
  expect_true(all(is.na(est$class[-nrow(est)])))
  expect_identical(sum(est$param == "ETA_CL"), 1L)
  # Addressable by name, like any other row.
  expect_identical(unname(ferx_coef(fit, label)), row$estimate)
  expect_identical(unname(ferx_se(fit, label)), row$se)
})

test_that("print() lists the override under its own heading (ferx-core #2833)", {
  skip_on_cran()
  fit <- mixture_override_fit()
  rows <- override_section(fit)
  expect_length(rows, 1L)
  expect_true(startsWith(rows, paste0("  ", label, " ")))
  expect_match(rows, sprintf("= %.6f", fit$mixture_overrides[[label]]), fixed = TRUE)
  expect_match(rows, sprintf("SE = %.6f", fit$se_mixture_overrides[[label]]), fixed = TRUE)
  # The engine's OMEGA row: a log-normal ETA's class variance reads as a CV%.
  expect_match(rows, "CV% = ", fixed = TRUE)
})

test_that("the override rows survive save / load (ferx-core #2833)", {
  skip_on_cran()
  fit <- mixture_override_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$mixture_override_params, fit$mixture_override_params)
  expect_identical(back$mixture_overrides, fit$mixture_overrides)
  expect_identical(back$se_mixture_overrides, fit$se_mixture_overrides)
  expect_identical(back$estimates, fit$estimates)
  expect_identical(override_section(back), override_section(fit))
})

test_that("a FIX override has no SE: `[FIX]` and `---` (ferx-core #2833)", {
  skip_on_cran()
  fit <- quiet(fitrx_edit_fit(mixture_override_fit(), function(w) {
    w$mixture_override_params[[1L]]$fixed <- TRUE
    w
  }))
  expect_true(fit$mixture_override_params$fixed)
  row <- fit$estimates[label, ]
  expect_true(row$fixed)
  expect_true(is.na(row$se))
  expect_true(is.na(row$rse_pct))
  expect_identical(row$estimate, fit$mixture_overrides[[label]])
  rows <- override_section(fit)
  expect_true(startsWith(rows, paste0("  ", label, " [FIX]")))
  expect_match(rows, "SE = ---", fixed = TRUE)
})

test_that("a held base does not hold its override (ferx-core #2833)", {
  skip_on_cran()
  # ETA_CL itself listed as held, with no class: the override names the same
  # kind and index, and only its class tells the two apart.
  full <- mixture_override_fit()
  fit <- quiet(fitrx_edit_fit(full, function(w) {
    w$held_at_bound <- list(list(
      name = "ETA_CL", kind = "Omega", index = 0L,
      packed_index = match("ETA_CL", rownames(full$cov_matrix)) - 1L,
      reason = "collapsed_at_guard"
    ))
    w
  }))
  expect_identical(fit$held_at_bound$parameter, "ETA_CL")
  expect_true(is.na(fit$held_at_bound$class))
  expect_true(is.na(fit$estimates["ETA_CL", "se"]))
  expect_identical(fit$estimates[label, "se"], full$se_mixture_overrides[[label]])
  expect_match(override_section(fit),
               sprintf("SE = %.6f", full$se_mixture_overrides[[label]]), fixed = TRUE)
})

test_that("a sigma override reads like its base sigma (ferx-core #2833)", {
  skip_on_cran()
  # The fixture's one override, relabelled on the wire as a class-2 override
  # of EPS: what a `sigma(2) EPS ~ ...` line records.
  fit <- quiet(fitrx_edit_fit(mixture_override_fit(), function(w) {
    w$mixture_override_params[[1L]]$kind <- "Sigma"
    w$mixture_override_params[[1L]]$name <- "EPS_MIX2"
    w
  }))
  expect_identical(fit$mixture_override_params$parameter, "EPS (class 2)")
  expect_identical(fit$mixture_override_params$kind, "sigma")
  row <- fit$estimates["EPS (class 2)", ]
  expect_identical(row$transform, fit$estimates["EPS", "transform"])
  expect_true(is.na(row$scale))
  expect_identical(row$class, 2L)
  rows <- override_section(fit)
  expect_length(rows, 1L)
  s <- fit$mixture_overrides[["EPS (class 2)"]]
  expect_identical(
    rows,
    sprintf("  %-20s %-14s = %.6f  (var = %.6f, CV%% = %.1f)  SE = %.6f",
            "EPS (class 2)", "[proportional]", s, s * s, s * 100,
            fit$se_mixture_overrides[["EPS (class 2)"]])
  )
})

# An additive ETA_CL, declared second, with a class-2 override, and a class-2
# override of the second sigma: the override rows must read their base's
# scale, linked theta and sigma type by the base's index. Both overrides sit
# on index 2, so the first row of each family is the wrong one to read.
additive_override_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      mod <- tempfile(fileext = ".ferx")
      writeLines(c(
        "[parameters]",
        "  theta TVCL(0.13, 0.001, 10)",
        "  theta TVV(8, 0.1, 100)",
        "  theta TVKA(1, 0.01, 50)",
        "  theta P1(0.5, 0.01, 0.99)",
        "  omega ETA_V ~ 0.02",
        "  omega ETA_CL ~ 0.0004",
        "  sigma PROP_ERR ~ 0.1 (sd)",
        "  sigma ADD_ERR ~ 0.2 (sd)",
        "[mixture]",
        "  nsub = 2",
        "  p(1) = P1",
        "  omega(2) ETA_CL ~ 0.002",
        "  sigma(2) ADD_ERR ~ 0.5 (sd)",
        "[individual_parameters]",
        "  CL = TVCL + ETA_CL",
        "  V  = TVV * exp(ETA_V)",
        "  KA = TVKA",
        "[structural_model]",
        "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
        "[error_model]",
        "  DV ~ combined(PROP_ERR, ADD_ERR)"
      ), mod)
      fit <<- quiet(ferx_fit(mod, ferx_example("warfarin")$data, method = "foce",
                             verbose = FALSE, covariance = TRUE,
                             settings = list(maxiter = 20L)))
    }
    fit
  }
})

test_that("an override row reads its base's scale, linked theta and sigma type", {
  skip_on_cran()
  fit <- additive_override_fit()
  expect_false(is.null(fit$cov_matrix))
  om <- "ETA_CL (class 2)"
  sg <- "ADD_ERR (class 2)"
  expect_identical(fit$mixture_override_params[, c("parameter", "kind", "index")],
                   data.frame(parameter = c(om, sg), kind = c("omega", "sigma"),
                              index = c(2L, 2L), stringsAsFactors = FALSE))
  # The premise: the two ETAs and the two sigmas differ, so reading the wrong
  # base shows.
  expect_identical(unname(fit$eta_param_types), c("log_normal", "additive"))
  expect_identical(unname(fit$eta_linked_theta), c("TVV", "TVCL"))
  expect_identical(unname(fit$sigma_types), c("proportional", "additive"))

  est <- fit$estimates
  expect_identical(est[om, "scale"], "additive")
  expect_identical(est[sg, "transform"], "additive")
  # A sigma override packs the log of the class sigma: SE = sigma * packed SE.
  expect_equal(fit$se_mixture_overrides[[sg]],
               fit$mixture_overrides[[sg]] * sqrt(fit$cov_matrix[sg, sg]),
               tolerance = 1e-8)
  expect_identical(est[sg, "se"], fit$se_mixture_overrides[[sg]])

  rows <- override_section(fit)
  expect_length(rows, 2L)
  # The omega override, on the additive scale of ETA_CL, against ETA_CL's theta.
  expect_true(startsWith(rows[1L], paste0("  ", om, " ")))
  expect_match(rows[1L], "[additive]", fixed = TRUE)
  expect_match(rows[1L], sprintf("at TVCL = %.4f", fit$theta[["TVCL"]]), fixed = TRUE)
  # The sigma override, as the additive ADD_ERR is printed.
  s <- fit$mixture_overrides[[sg]]
  expect_identical(
    rows[2L],
    sprintf("  %-20s %-14s = %.6f  (SD = %.4f, var = %.6f)  SE = %.6f",
            sg, "[additive]", s, s, s * s, fit$se_mixture_overrides[[sg]])
  )
})

test_that("a bundle from before the descriptors names no override row", {
  skip_on_cran()
  full <- mixture_override_fit()
  old <- quiet(fitrx_edit_fit(full, function(w) {
    w$mixture_override_params <- NULL
    w$se_mixture_overrides <- NULL
    w
  }))
  # The values stay, for the from-fit paths; nothing names them.
  expect_identical(unname(old$mixture_overrides), unname(full$mixture_overrides))
  expect_null(names(old$mixture_overrides))
  expect_null(old$se_mixture_overrides)
  expect_identical(nrow(old$mixture_override_params), 0L)
  expect_named(old$mixture_override_params,
               c("parameter", "kind", "index", "class", "fixed"))
  expect_false(label %in% old$estimates$param)
  expect_identical(old$estimates$param, full$estimates$param[-nrow(full$estimates)])
  expect_length(override_section(old), 0L)
})

test_that("a model without overrides has no override surface", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  expect_null(fit$mixture_overrides)
  expect_null(fit$se_mixture_overrides)
  expect_identical(nrow(fit$mixture_override_params), 0L)
  expect_true(all(is.na(fit$estimates$class)))
  expect_length(override_section(fit), 0L)
})
