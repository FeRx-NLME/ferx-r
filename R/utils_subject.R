# Stable per-subject index for datasets that reuse a subject ID.
#
# ferx-core (like NONMEM, which processes records sequentially) starts a new
# subject whenever the ID differs from the immediately preceding record. A
# dataset can therefore contain the same textual ID in two non-contiguous
# blocks - e.g. a second cohort that reuses IDs 12/13/14 - and those are two
# distinct subjects that happen to share an ID. The per-subject tables
# (`fit$ebe_etas`, `fit$individual_estimates`) then carry one row per subject in
# subject order, while `fit$sdtab` / the raw data carry each subject's records
# contiguously in that same order. Keying joins on the raw ID collapses the
# duplicates; keying on this ordinal index does not.
#
# Returns a 1-based integer, one per input record, incrementing at every change
# of ID relative to the previous record. Empty input returns integer(0).
.ferx_subject_index <- function(ids) {
  ids <- as.character(ids)
  n <- length(ids)
  if (n == 0L) return(integer(0))
  brk <- c(TRUE, ids[-1L] != ids[-n])
  # An NA ID compares NA to its neighbour; treat any NA-involved comparison as a
  # break so a lone NA record is never silently merged into an adjacent block.
  brk[is.na(brk)] <- TRUE
  cumsum(brk)
}

# Number of distinct subjects implied by a per-record ID vector.
.ferx_n_subjects <- function(ids) {
  s <- .ferx_subject_index(ids)
  if (length(s)) max(s) else 0L
}
