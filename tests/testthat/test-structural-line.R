# The `Structural:` line of print.ferx_fit() and the `Structure:` line of
# print(summary(fit)) name each theta block once with its free-coefficient
# count, whatever its size (#538). The line lists the model's declarations;
# the THETA table, not this line, is where ferx-core's
# THETA_BLOCK_COMPACT_MIN decides between per-level rows and a summary.

structural_labels <- getFromNamespace("ferx_rust_structural_theta_labels", "ferx")

# The theta part of each line: the text inside `compartment-free  (...)`.
structural_from_print <- function(fit) {
  line <- grep("^  Structural:", utils::capture.output(print(fit)), value = TRUE)
  expect_length(line, 1L)
  sub("^  Structural:  ", "", line)
}
structural_from_summary <- function(fit) {
  line <- grep("^Structure:", utils::capture.output(print(summary(fit))), value = TRUE)
  expect_length(line, 1L)
  sub("  \\|  IIV: .*$", "", sub("^Structure: ", "", line))
}

level_names <- function(block, n) sprintf("%s[STUDY=%d,TIME=0]", block, seq_len(n))

fake_structured_fit <- function(theta_names) {
  make_fake_fit(
    theta = stats::setNames(seq_along(theta_names) / 10, theta_names),
    omega = matrix(0.1, 1, 1),
    model_structure = list(model_type = "compartment-free", theta_names = theta_names,
                           iiv = "ETA_E0", iov = character(), residual = "additive")
  )
}

test_that("mbma_placebo names PLACEBO once, with its free-coefficient count", {
  skip_on_cran()
  ex <- ferx_example("mbma_placebo")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE))
  n <- sum(grepl("^PLACEBO\\[", names(fit$theta)))
  expect_identical(n, 18L)
  s <- structural_from_print(fit)
  expect_identical(s, "compartment-free  (TVE0, EMAX, ED50, ET50, PLACEBO[STUDY, TIME] (18 free coefficients))")
  expect_false(grepl("STUDY=", s, fixed = TRUE))
  expect_identical(structural_from_summary(fit), s)
  # Still a reload away: the bundle carries the names as a list.
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  loaded <- suppressWarnings(ferx_load_fit(f))
  expect_identical(structural_from_print(loaded), s)
})

test_that("a block is named once whether it has 5 or 50 coefficients", {
  for (n in c(5L, 50L)) {
    fit <- fake_structured_fit(c("TVE0", level_names("PLACEBO", n), "EMAX"))
    s <- structural_from_print(fit)
    expect_identical(s, sprintf(
      "compartment-free  (TVE0, PLACEBO[STUDY, TIME] (%d free coefficients), EMAX)", n))
    expect_identical(structural_from_summary(fit), s, info = n)
  }
})

test_that("scalar thetas are unchanged; keyless and adjacent blocks are named apart", {
  expect_identical(structural_labels(c("TVCL", "TVV", "TVKA")), c("TVCL", "TVV", "TVKA"))
  expect_identical(structural_labels(character()), character())
  expect_identical(
    structural_labels(c("W[1]", "W[2]", "W[3]", "A[STUDY=1]", "B[STUDY=1]", "B[STUDY=2]")),
    c("W (3 free coefficients)", "A[STUDY] (1 free coefficient)",
      "B[STUDY] (2 free coefficients)")
  )
})
