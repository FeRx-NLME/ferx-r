# The compiled package statically links ferx-search (LGPL-3.0-or-later) and
# NLopt's Luksan routines (LGPL-2.1-or-later), so the installed package must
# carry their license texts and every linked crate's notices (#773).

notice_file <- function(...) {
  path <- system.file(..., package = "ferx")
  expect_true(nzchar(path) && file.exists(path), label = file.path(...))
  path
}

test_that("the installed package carries the license texts and notices", {
  for (f in c("COPYRIGHTS", "THIRD-PARTY-NOTICES.txt")) {
    expect_gt(file.size(notice_file(f)), 0)
  }
  lgpl3 <- readLines(notice_file("licenses", "LGPL-3.0.txt"))
  expect_true(any(grepl("GNU LESSER GENERAL PUBLIC LICENSE", lgpl3, fixed = TRUE)))
  expect_true(any(grepl("Version 3, 29 June 2007", lgpl3, fixed = TRUE)))
  gpl3 <- readLines(notice_file("licenses", "GPL-3.0.txt"))
  expect_true(any(grepl("GNU GENERAL PUBLIC LICENSE", gpl3, fixed = TRUE)))
  lgpl21 <- readLines(notice_file("licenses", "LGPL-2.1.txt"))
  expect_true(any(grepl("Version 2.1, February 1999", lgpl21, fixed = TRUE)))
})

test_that("COPYRIGHTS names each license and the rebuild recipe", {
  text <- paste(readLines(notice_file("COPYRIGHTS")), collapse = "\n")
  for (needle in c("MIT License", "LGPL-3.0-or-later", "LGPL-2.1-or-later",
                   "Pharmpy", "src/algs/luksan", "crates/ferx-search",
                   "R CMD INSTALL", "THIRD-PARTY-NOTICES.txt")) {
    expect_true(grepl(needle, text, fixed = TRUE), label = needle)
  }
})

test_that("THIRD-PARTY-NOTICES lists ferx-search, NLopt and the pinned commit", {
  text <- readLines(notice_file("THIRD-PARTY-NOTICES.txt"))
  expect_true(any(grepl("^  ferx-search [0-9.]+  LGPL-3\\.0-or-later$", text)))
  expect_true(any(grepl("^  ferx-core [0-9.]+  MIT$", text)))
  expect_true(any(grepl("^  nlopt ", text)))
  expect_true(any(grepl("src/algs/luksan", text, fixed = TRUE)))
  expect_true(any(grepl(
    "^  https://github.com/FeRx-NLME/ferx-core/tree/[0-9a-f]{40}$", text
  )))
})

test_that("DESCRIPTION points at inst/COPYRIGHTS", {
  copyright <- utils::packageDescription("ferx")$Copyright
  expect_false(is.null(copyright))
  expect_match(copyright, "inst/COPYRIGHTS", fixed = TRUE)
  expect_match(copyright, "LGPL-3.0-or-later", fixed = TRUE)
})
