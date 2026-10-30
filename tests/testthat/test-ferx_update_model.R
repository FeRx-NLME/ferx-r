# ferx_update_model() / ferx_save_fit(write_model = TRUE) (#678): the fitted
# model written back with its final estimates as inits.

# Warfarin with a block omega, a FIX theta, an `(sd)` omega and an `(sd)`
# sigma: every declaration form the round trip has to keep.
update_model_fixture <- function() {
  write_test_model(list(
    parameters = c(
      "  theta TVCL(0.134, 0.001, 10.0)   # L/h",
      "  # volume",
      "  theta TVV(8.1, 0.1, 500.0)",
      "  theta TVKA(1.0, FIX)",
      "",
      "  block_omega (ETA_CL, ETA_V) = [0.07, 0.02, 0.02] # CL-V block",
      "  omega ETA_KA ~ 0.6 (sd)  // absorption",
      "",
      "  sigma PROP_ERR ~ 0.01 (sd)  # proportional"
    ),
    individual_parameters = c(
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV  * exp(ETA_V)",
      "  KA = TVKA * exp(ETA_KA)"
    ),
    structural_model = "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    error_model = "  DV ~ proportional(PROP_ERR)"
  ))
}

update_model_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      model <- update_model_fixture()
      fit <<- ferx_fit(model, ferx_example("warfarin")$data, method = "focei",
                       verbose = FALSE, covariance = FALSE,
                       settings = list(maxiter = 30L))
    }
    fit
  }
})

test_that("a refit from the written model starts at the fit's OFV (#678)", {
  fit <- update_model_fit()
  out <- tempfile(fileext = ".ferx")
  on.exit(unlink(out), add = TRUE)
  expect_identical(ferx_update_model(fit, out), out)

  refit <- ferx_fit(out, ferx_example("warfarin")$data, method = "focei",
                    verbose = FALSE, covariance = FALSE,
                    settings = list(maxiter = 0L))
  expect_equal(refit$ofv, fit$ofv, tolerance = 1e-6)
  expect_equal(refit$theta, fit$theta, tolerance = 1e-10)
  expect_equal(refit$omega, fit$omega, tolerance = 1e-10)
  expect_equal(refit$sigma, fit$sigma, tolerance = 1e-10)
})

test_that("the written model keeps bounds, FIX, (sd), blocks and comments (#678)", {
  fit <- update_model_fit()
  text <- ferx_update_model(fit)
  expect_type(text, "character")
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  par <- trimws(lines)

  # Every other block is untouched.
  src <- readLines(fit$model_path)
  expect_identical(lines[seq(which(lines == "[individual_parameters]"), length(lines))],
                   src[seq(which(src == "[individual_parameters]"), length(src))])

  # The FIX theta is left as declared; the estimated ones keep their bounds,
  # and a whole-line comment stays where it was.
  expect_true("theta TVKA(1.0, FIX)" %in% par)
  expect_true(any(grepl("^theta TVCL\\([^,]+, 0\\.001, 10\\.0\\)   # L/h$", par)))
  expect_true(any(grepl("^theta TVV\\([^,]+, 0\\.1, 500\\.0\\)$", par)))
  expect_identical(par[which(par == "# volume") + 1L],
                   grep("^theta TVV", par, value = TRUE))
  tvcl <- as.numeric(sub("^theta TVCL\\(([^,]+),.*$", "\\1",
                         grep("^theta TVCL", par, value = TRUE)))
  expect_equal(tvcl, fit$theta[["TVCL"]], tolerance = 1e-12)

  # The block keeps its structure and takes the fitted lower triangle.
  blk <- grep("^block_omega \\(ETA_CL, ETA_V\\) = \\[", par, value = TRUE)
  expect_length(blk, 1L)
  expect_match(blk, "\\] # CL-V block$")
  tri <- as.numeric(strsplit(sub("^.*\\[(.*)\\].*$", "\\1", blk), ",")[[1]])
  om <- fit$omega
  expect_equal(tri, c(om[1, 1], om[2, 1], om[2, 2]), tolerance = 1e-12)

  # `(sd)` declarations take the SD.
  ka <- grep("^omega ETA_KA ~ ", par, value = TRUE)
  expect_match(ka, "\\(sd\\)  // absorption$")
  expect_equal(as.numeric(sub("^omega ETA_KA ~ (\\S+) \\(sd\\).*$", "\\1", ka)),
               sqrt(om[3, 3]), tolerance = 1e-12)
  sg <- grep("^sigma PROP_ERR ~ ", par, value = TRUE)
  expect_match(sg, "\\(sd\\)  # proportional$")
  expect_equal(as.numeric(sub("^sigma PROP_ERR ~ (\\S+) \\(sd\\).*$", "\\1", sg)),
               unname(fit$sigma[[1]]), tolerance = 1e-12)
})

test_that("ferx_update_model() takes the fitted text, not the edited file (#678)", {
  fit <- update_model_fit()
  before <- ferx_update_model(fit)
  copy <- fit
  copy$model_path <- tempfile(fileext = ".ferx")
  writeLines("not a model", copy$model_path)
  on.exit(unlink(copy$model_path), add = TRUE)
  expect_identical(ferx_update_model(copy), before)
})

test_that("ferx_update_model() refuses to overwrite unless asked (#678)", {
  fit <- update_model_fit()
  out <- tempfile(fileext = ".ferx")
  on.exit(unlink(out), add = TRUE)
  writeLines("keep me", out)
  expect_error(ferx_update_model(fit, out), "already exists")
  expect_identical(readLines(out), "keep me")
  ferx_update_model(fit, out, overwrite = TRUE)
  expect_identical(paste(readLines(out), collapse = "\n"),
                   sub("\n$", "", ferx_update_model(fit)))
  expect_error(ferx_update_model(list()), "ferx_fit object")
  expect_error(ferx_update_model(fit, NA_character_), "`path`")
})

test_that("ferx_save_fit(write_model = TRUE) writes <name>_final.ferx (#678)", {
  fit <- update_model_fit()
  dir <- tempfile("fitrx_dir_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  ferx_save_fit(fit, file.path(dir, "run1"), write_model = TRUE)
  expect_true(file.exists(file.path(dir, "run1.fitrx")))
  final <- file.path(dir, "run1_final.ferx")
  expect_true(file.exists(final))
  expect_identical(rawToChar(readBin(final, "raw", file.size(final))),
                   ferx_update_model(fit))

  # The default writes no model.
  ferx_save_fit(fit, file.path(dir, "run2.fitrx"))
  expect_false(file.exists(file.path(dir, "run2_final.ferx")))
  expect_error(ferx_save_fit(fit, file.path(dir, "run3"), write_model = NA),
               "write_model")
})

test_that("ferx_save_fit(write_model = TRUE) never overwrites the fit's model (#678)", {
  fit <- update_model_fit()
  dir <- tempfile("fitrx_dir_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  copy <- fit
  copy$model_path <- file.path(dir, "run1_final.ferx")
  file.copy(fit$model_path, copy$model_path)
  expect_error(ferx_save_fit(copy, file.path(dir, "run1.fitrx"), write_model = TRUE),
               "overwrite the fit's own model file")
  expect_false(file.exists(file.path(dir, "run1.fitrx")))
})

test_that("a reloaded fit writes the same model (#678)", {
  fit <- update_model_fit()
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)
  ferx_save_fit(fit, tmp)
  expect_identical(ferx_update_model(ferx_load_fit(tmp)), ferx_update_model(fit))
})

# A free `block_sigma` (ferx-core #847) estimates both SDs and their
# correlation; the engine's edit leaves it alone, so the glue writes it.
block_sigma_case <- function(sigma_line) {
  dir <- tempfile("ferx-update-rho-")
  dir.create(dir)
  model <- file.path(dir, "rho.ferx")
  data  <- file.path(dir, "rho.csv")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(1.0, 0.01, 10.0)",
    "  theta TVV(10.0, 0.1, 100.0)",
    "  omega ETA_CL ~ 0.04",
    sigma_line,
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  pk one_cpt_iv(cl=CL, v=V)",
    "[error_model]",
    "  DV ~ combined(PROP_ERR, ADD_ERR)",
    "[fit_options]",
    "  method = focei"
  ), model)
  write.csv(data.frame(
    ID   = rep(1:3, each = 4),
    TIME = rep(c(0, 1, 3, 8), 3),
    DV   = c(0, 9.1, 6.8, 3.0, 0, 8.4, 6.2, 2.6, 0, 9.7, 7.3, 3.4),
    EVID = rep(c(1, 0, 0, 0), 3),
    AMT  = rep(c(100, 0, 0, 0), 3),
    CMT  = 1,
    MDV  = rep(c(1, 0, 0, 0), 3)
  ), data, row.names = FALSE, quote = FALSE)
  list(model = model, data = data,
       fit = ferx_fit(model, data, verbose = FALSE, covariance = FALSE))
}

test_that("a free block_sigma takes the fitted residual covariance (#678)", {
  case <- block_sigma_case(
    "  block_sigma (PROP_ERR, ADD_ERR) = [0.04, 0.10, 0.30]  # residual block"
  )
  fit <- case$fit
  out <- tempfile(fileext = ".ferx")
  on.exit(unlink(c(out, dirname(case$model)), recursive = TRUE), add = TRUE)
  ferx_update_model(fit, out)

  line <- grep("^ *block_sigma", readLines(out), value = TRUE)
  expect_match(line, "^  block_sigma \\(PROP_ERR, ADD_ERR\\) = \\[.*\\]  # residual block$")
  tri <- as.numeric(strsplit(sub("^.*\\[(.*)\\].*$", "\\1", line), ",")[[1]])
  sd <- unname(fit$sigma)
  rho <- fit$residual_correlations$rho
  expect_equal(tri, c(sd[1]^2, rho * sd[1] * sd[2], sd[2]^2), tolerance = 1e-12)

  refit <- ferx_fit(out, case$data, verbose = FALSE, covariance = FALSE,
                    settings = list(maxiter = 0L))
  expect_equal(refit$ofv, fit$ofv, tolerance = 1e-6)
  expect_equal(refit$sigma, fit$sigma, tolerance = 1e-10)
  expect_equal(refit$residual_correlations$rho, rho, tolerance = 1e-10)
})

test_that("a FIX block_sigma is left as declared (#678)", {
  decl <- "  block_sigma (PROP_ERR, ADD_ERR) = [0.04, 0.10, 0.30] FIX"
  case <- block_sigma_case(decl)
  on.exit(unlink(dirname(case$model), recursive = TRUE), add = TRUE)
  lines <- strsplit(ferx_update_model(case$fit), "\n", fixed = TRUE)[[1]]
  expect_true(decl %in% lines)
})

test_that("comments on a declaration spread over lines are kept (#678)", {
  fit <- update_model_fit()
  copy <- fit
  copy$model_text <- NULL
  copy$model_path <- write_test_model(list(
    parameters = c(
      "  theta TVCL(0.134, 0.001, 10.0)",
      "  theta TVV(8.1, 0.1, 500.0)",
      "  theta TVKA(1.0, FIX)",
      "  block_omega (ETA_CL, ETA_V) = [0.07,   # var CL",
      "                                 0.02, 0.02]  # cov, var V",
      "  omega ETA_KA ~ 0.6 (sd)",
      "  sigma PROP_ERR ~ 0.01 (sd)"
    ),
    individual_parameters = c(
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV  * exp(ETA_V)",
      "  KA = TVKA * exp(ETA_KA)"
    ),
    structural_model = "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    error_model = "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(copy$model_path), add = TRUE)
  copy$model_hash <- NULL
  lines <- strsplit(ferx_update_model(copy), "\n", fixed = TRUE)[[1]]
  blk <- grep("block_omega", lines, value = TRUE)
  expect_length(blk, 1L)
  expect_match(blk, "^  block_omega \\(ETA_CL, ETA_V\\) = \\[[^]]+\\]   # var CL  # cov, var V$")
  expect_false(any(grepl("^ +0\\.02, 0\\.02\\]", lines)))
  expect_identical(lines[which(grepl("block_omega", lines)) + 1L],
                   grep("^  omega ETA_KA", lines, value = TRUE))
})

test_that("a block_sigma written over lines is left as written, with a warning (#678)", {
  case <- block_sigma_case("  block_sigma (PROP_ERR, ADD_ERR) = [0.04, 0.10, 0.30]")
  on.exit(unlink(dirname(case$model), recursive = TRUE), add = TRUE)
  copy <- case$fit
  copy$model_text <- NULL
  copy$model_hash <- NULL
  src <- readLines(case$model)
  at <- grep("block_sigma", src)
  src <- append(src[-at], c("  block_sigma (PROP_ERR, ADD_ERR) = [0.04,",
                            "                                    0.10, 0.30]"),
                after = at - 1L)
  writeLines(src, case$model)
  expect_warning(text <- ferx_update_model(copy), "spans more than one line")
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  expect_true(all(src[at + 0:1] %in% lines))
})
