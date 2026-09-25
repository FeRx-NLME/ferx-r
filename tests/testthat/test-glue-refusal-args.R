# A refused call must not keep its argument vectors alive (#394).
#
# extendr converts every argument SEXP to an `Robj` - protecting it - in a
# frame of its own, outside the glue function. While `entry()` raised with
# `Rf_error()`, that frame was longjmp'd over and the protection never came
# off: every refused call that was handed fresh 200,000-element vectors kept
# them for the rest of the session (~6 MB per call here). `entry()` now unwinds
# into extendr's wrapper instead, which drops the `Robj`s before it raises.
#
# The measurement is gc()'s Vcells, not RSS: R's own count of live vector
# memory, which a `gc(full = TRUE)` makes exact. The first batch is allowed a
# one-off rise (the allocator settling on its high-water mark); what must not
# happen is growth from one batch to the next. The yardstick is what one
# argument vector costs when it *is* kept alive, measured the same way, so the
# bound does not depend on the platform's string sizes: leaking, 20 refusals
# retain 20 or more of them (130 MB, ~87 vectors, on a 0.3.0 build); not
# leaking, a busy session may still drift by a few MB (4 MB, under 3 vectors,
# was seen under the full suite). The bound is 10.

refuse_with_big_args <- function(i) {
  big <- big_arg(i)
  tryCatch(
    ferx:::ferx_rust_fit(
      "no-such-model.ferx", "no-such-data.csv", character(), "", "", "",
      1L, "", "", "", big, big
    ),
    error = function(e) conditionMessage(e)
  )
}

vcells_mb <- function() gc(full = TRUE)["Vcells", 2]

big_arg <- function(i, n = 2e5) as.character(seq_len(n) + i * n)

test_that("a refused call does not retain its argument vectors", {
  skip_on_cran()
  msg <- refuse_with_big_args(0)
  expect_match(msg, "^Error parsing model: ")

  before <- vcells_mb()
  kept <- lapply(1:4, function(i) big_arg(1000 + i))
  one_vector <- (vcells_mb() - before) / 4
  rm(kept)
  expect_gt(one_vector, 1)

  after <- vapply(1:3, function(batch) {
    for (i in seq_len(10)) refuse_with_big_args(batch * 100 + i)
    vcells_mb()
  }, numeric(1))

  expect_lt(
    after[3] - after[1], 10 * one_vector,
    label = sprintf("Vcells growth over 20 refusals (MB; one vector = %.1f MB)", one_vector)
  )
})
