# A refused call must not keep its argument vectors alive (#394).
#
# extendr converts every argument SEXP to an `Robj` - protecting it - in a
# frame of its own, outside the glue function. While `entry()` raised with
# `Rf_error()`, that frame was longjmp'd over and the protection never came
# off: every refused call that was handed fresh 200,000-element vectors kept
# them for the rest of the session (~5 MB per call here). `entry()` now unwinds
# into extendr's wrapper instead, which drops the `Robj`s before it raises.
#
# The measurement is gc()'s Vcells, not RSS: R's own count of live vector
# memory, which a `gc(full = TRUE)` makes exact. The first batch is allowed a
# one-off rise (the allocator settling on its high-water mark); what must not
# happen is growth from one batch to the next.

refuse_with_big_args <- function(i, n = 2e5) {
  big <- as.character(seq_len(n) + i * n)
  tryCatch(
    ferx:::ferx_rust_fit(
      "no-such-model.ferx", "no-such-data.csv", character(), "", "", "",
      1L, "", "", "", big, big
    ),
    error = function(e) conditionMessage(e)
  )
}

vcells_mb <- function() gc(full = TRUE)["Vcells", 2]

test_that("a refused call does not retain its argument vectors", {
  skip_on_cran()
  msg <- refuse_with_big_args(0)
  expect_match(msg, "^Error parsing model: ")

  baseline <- vcells_mb()
  after <- vapply(1:3, function(batch) {
    for (i in seq_len(10)) refuse_with_big_args(batch * 100 + i)
    vcells_mb()
  }, numeric(1))

  # Leaking, each batch of 10 adds ~50 MB; not leaking, batches 2 and 3 add
  # nothing measurable.
  expect_lt(after[3] - after[1], 2, label = "Vcells growth, batch 1 -> 3 (MB)")
  expect_lt(after[1] - baseline, 20, label = "Vcells after first batch (MB)")
})

