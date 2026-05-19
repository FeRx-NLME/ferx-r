WARFARIN_MOD <- "$PROB warfarin one-cpt oral
$PARAM TVCL = 0.2, TVV = 10.0, TVKA = 1.5
$CMT GUT CENT
$OMEGA @labels ETA_CL ETA_V ETA_KA
0.09 0.04 0.30
$SIGMA @labels PROP_ERR
0.02
$MAIN
double CL = TVCL * exp(ETA_CL);
double V  = TVV  * exp(ETA_V);
double KA = TVKA * exp(ETA_KA);
$PKMODEL ncmt=1, depot=TRUE
$TABLE
capture IPRED = CENT/V;
capture Y     = IPRED * (1 + EPS(1));
"

write_mod_to_temp <- function(text, ext = ".mod") {
  path <- tempfile(fileext = ext)
  writeLines(text, path)
  path
}

test_that("ferx_translate_mrgsolve() returns translated .ferx text", {
  path <- write_mod_to_temp(WARFARIN_MOD)
  on.exit(unlink(path), add = TRUE)

  out <- ferx_translate_mrgsolve(path)
  expect_type(out, "character")
  expect_length(out, 1L)
  expect_match(out, "[parameters]", fixed = TRUE)
  expect_match(out, "theta TVCL(0.2)", fixed = TRUE)
  expect_match(out, "pk one_cpt_oral(cl=CL, v=V, ka=KA)", fixed = TRUE)
  expect_match(out, "DV ~ proportional(PROP_ERR)", fixed = TRUE)
})

test_that("ferx_translate_mrgsolve() writes to disk when output= is given", {
  in_path  <- write_mod_to_temp(WARFARIN_MOD)
  out_path <- tempfile(fileext = ".ferx")
  on.exit(unlink(c(in_path, out_path)), add = TRUE)

  res <- ferx_translate_mrgsolve(in_path, output = out_path)
  expect_true(file.exists(out_path))
  written <- paste(readLines(out_path), collapse = "\n")
  expect_match(written, "pk one_cpt_oral", fixed = TRUE)
  # The return value is the translated text, returned invisibly.
  expect_type(res, "character")
})

test_that("ferx_translate_mrgsolve() refuses to overwrite without overwrite=TRUE", {
  in_path  <- write_mod_to_temp(WARFARIN_MOD)
  out_path <- tempfile(fileext = ".ferx")
  writeLines("# placeholder", out_path)
  on.exit(unlink(c(in_path, out_path)), add = TRUE)

  expect_error(
    ferx_translate_mrgsolve(in_path, output = out_path),
    "already exists"
  )
  expect_silent(ferx_translate_mrgsolve(in_path, output = out_path, overwrite = TRUE))
})

test_that("ferx_translate_mrgsolve() errors for non-existent file", {
  expect_error(
    ferx_translate_mrgsolve("/no/such/file.mod"),
    "File not found"
  )
})

test_that("ferx_translate_mrgsolve() surfaces translator errors", {
  bad_mod <- "$PARAM TVCL = 1\n$MAIN\nfor (int i = 0; i < 3; i++) {}\n"
  path <- write_mod_to_temp(bad_mod)
  on.exit(unlink(path), add = TRUE)

  expect_error(
    ferx_translate_mrgsolve(path),
    "mrgsolve translation failed"
  )
})

test_that("ferx_translate_mrgsolve() warns on unexpected extension", {
  path <- write_mod_to_temp(WARFARIN_MOD, ext = ".txt")
  on.exit(unlink(path), add = TRUE)

  expect_warning(
    ferx_translate_mrgsolve(path),
    "does not end in .cpp or .mod"
  )
})
