# ferx_update_model() / ferx_save_fit(write_model = TRUE) (#678): the fitted
# model written back with its final estimates as inits.

# Warfarin with a block omega, a FIX theta, an `(sd)` omega and an `(sd)`
# sigma: every declaration form the round trip has to keep.
update_model_fixture <- function() {
  write_test_model(list(
    parameters = c(
      "  theta TVCL(0.134, 0.001, 10.0)",
      "  # volume",
      "  theta TVV(8.1, 0.1, 500.0)",
      "  theta TVKA(1.0, FIX)",
      "",
      "  block_omega (ETA_CL, ETA_V) = [0.07, 0.02, 0.02]",
      "  omega ETA_KA ~ 0.6 (sd)",
      "",
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
  expect_true(any(grepl("^theta TVCL\\([^,]+, 0\\.001, 10\\.0\\)$", par)))
  expect_true(any(grepl("^theta TVV\\([^,]+, 0\\.1, 500\\.0\\)$", par)))
  expect_identical(par[which(par == "# volume") + 1L],
                   grep("^theta TVV", par, value = TRUE))
  tvcl <- as.numeric(sub("^theta TVCL\\(([^,]+),.*$", "\\1",
                         grep("^theta TVCL", par, value = TRUE)))
  expect_equal(tvcl, fit$theta[["TVCL"]], tolerance = 1e-12)

  # The block keeps its structure and takes the fitted lower triangle.
  blk <- grep("^block_omega \\(ETA_CL, ETA_V\\) = \\[", par, value = TRUE)
  expect_length(blk, 1L)
  tri <- as.numeric(strsplit(sub("^.*\\[(.*)\\].*$", "\\1", blk), ",")[[1]])
  om <- fit$omega
  expect_equal(tri, c(om[1, 1], om[2, 1], om[2, 2]), tolerance = 1e-12)

  # `(sd)` declarations take the SD.
  ka <- grep("^omega ETA_KA ~ ", par, value = TRUE)
  expect_match(ka, "\\(sd\\)$")
  expect_equal(as.numeric(sub("^omega ETA_KA ~ (\\S+) \\(sd\\)$", "\\1", ka)),
               sqrt(om[3, 3]), tolerance = 1e-12)
  sg <- grep("^sigma PROP_ERR ~ ", par, value = TRUE)
  expect_equal(as.numeric(sub("^sigma PROP_ERR ~ (\\S+) \\(sd\\)$", "\\1", sg)),
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
