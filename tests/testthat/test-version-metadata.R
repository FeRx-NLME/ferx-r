# `DESCRIPTION`'s `Version` and `NEWS.md`'s top heading are two statements of
# the same fact, and they disagreed for ten commits (#296): DESCRIPTION said
# `0.3.0` while NEWS.md already said `0.3.0.9000`, so `packageVersion("ferx")`
# returned one string for the v0.3.0 tag and for every development commit
# after it - including the ones that moved the pinned ferx-core revision
# across a breaking parser change.
#
# Nothing could see that. The two banner tests in test-ferx_fit.R build their
# expectation from `packageVersion("ferx")` itself, so they pass against
# whatever DESCRIPTION happens to say; they check the *rendering*, not the
# value. This file compares the two files to each other.
#
# Both are read from the installed package, so this holds for whatever
# `R CMD check` actually built, not for whatever is in the working tree.

# The version token from a `# ferx <version>` heading, with any trailing
# parenthetical ("(development version)") removed. NA when the line is not a
# ferx heading at all.
.news_heading_version <- function(line) {
  m <- regmatches(line, regexec("^#\\s+ferx\\s+([^\\s(]+)", line, perl = TRUE))[[1]]
  if (length(m) < 2L) return(NA_character_)
  trimws(m[[2]])
}

test_that(".news_heading_version reads the forms NEWS.md uses", {
  # The normalisation itself, so a failure below is about the metadata rather
  # than about this helper.
  expect_identical(.news_heading_version("# ferx 0.4.0"), "0.4.0")
  expect_identical(
    .news_heading_version("# ferx 0.3.0.9000 (development version)"),
    "0.3.0.9000"
  )
  expect_identical(.news_heading_version("## Bug fixes"), NA_character_)
  expect_identical(.news_heading_version("# something else"), NA_character_)
})

test_that("NEWS.md's top heading is the DESCRIPTION version", {
  news_path <- system.file("NEWS.md", package = "ferx")
  skip_if(!nzchar(news_path), "NEWS.md is not part of the installed package")

  headings <- Filter(
    function(v) !is.na(v),
    vapply(readLines(news_path, warn = FALSE), .news_heading_version,
           character(1L), USE.NAMES = FALSE)
  )
  expect_gt(length(headings), 0L)

  desc_path <- system.file("DESCRIPTION", package = "ferx")
  desc_version <- unname(read.dcf(desc_path, fields = "Version")[1L, 1L])

  expect_identical(
    headings[[1L]], desc_version,
    info = paste0(
      "NEWS.md's top heading says ", headings[[1L]],
      " and DESCRIPTION says ", desc_version,
      ". They name the same release; bump both (#296)."
    )
  )
})
