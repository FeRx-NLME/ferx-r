test_that("optional scalar wrappers collapse missing values to NULL", {
  for (f in list(ferx:::.fitrx_opt_num, ferx:::.fitrx_opt_int, ferx:::.fitrx_opt_chr,
                 ferx:::.fitrx_unwrap_opt_num, ferx:::.fitrx_unwrap_opt_int,
                 ferx:::.fitrx_unwrap_opt_chr)) {
    expect_null(f(NULL))
    expect_null(f(character(0)))
  }
  expect_identical(ferx:::.fitrx_opt_num(2.5), 2.5)
  expect_null(ferx:::.fitrx_opt_num(NA))
  expect_identical(ferx:::.fitrx_opt_int(5L), 5L)
  expect_null(ferx:::.fitrx_opt_int(NA))
  expect_identical(ferx:::.fitrx_opt_chr("x"), "x")
  expect_null(ferx:::.fitrx_opt_chr(""))     # empty string is "missing"
  expect_null(ferx:::.fitrx_opt_chr(NA))
})
test_that("optional vector wrappers keep non-empty vectors and drop empties", {
  expect_identical(ferx:::.fitrx_opt_num_vec(c(1, 2, 3)), c(1, 2, 3))
  expect_null(ferx:::.fitrx_opt_num_vec(numeric(0)))
  expect_null(ferx:::.fitrx_opt_num_vec(NULL))
  expect_identical(ferx:::.fitrx_unwrap_opt_num_vec(list(1, 2)), c(1, 2))
  expect_null(ferx:::.fitrx_unwrap_opt_num_vec(NULL))
})

test_that(".fitrx_unwrap_nullable_vec() maps JSON nulls to NA instead of dropping them", {
  # A Rust `Vec<Option<T>>` round-trips as a list with NULL holes. unlist()
  # would drop them and shift every later element onto the wrong kappa - which
  # for a two-kappa model would silently report the second kappa's weight on
  # the first.
  expect_equal(.fitrx_unwrap_opt_chr_vec(list(NULL, "NARM")),
               c(NA_character_, "NARM"))
  expect_equal(.fitrx_unwrap_nullable_num_vec(list(NULL, 400)),
               c(NA_real_, 400))
  expect_null(.fitrx_unwrap_opt_chr_vec(NULL))
  expect_null(.fitrx_unwrap_opt_chr_vec(list()))
  # A scalar (auto_unbox on a length-1 field written by an older writer).
  expect_equal(.fitrx_unwrap_opt_chr_vec("NARM"), "NARM")
})
