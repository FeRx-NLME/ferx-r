# -- progress modes (#677) ---------------------------------------------------
#
# `progress` on the bootstrap and on every search takes TRUE / FALSE or one of
# two named forms:
#
#   "bar"   the interactive display: a progress bar for the bootstrap, the
#           engine's per-step summary on stderr for a search.
#   "lines" plain text on stdout, one record per line, written and flushed as
#           each fit finishes - one line per bootstrap replicate, one line per
#           search candidate - so it streams into a log viewer or a job card
#           that shows a script's stdout.
#
# TRUE picks between them by `interactive()`: a console gets the bar, an
# `Rscript` run gets the lines, since a bar redrawn with carriage returns is
# noise in a log. FALSE (the non-interactive default) prints nothing.
.ferx_progress_mode <- function(progress, what, interactive = base::interactive()) {
  if (is.logical(progress) && length(progress) == 1L && !is.na(progress)) {
    if (!progress) return("off")
    return(if (isTRUE(interactive)) "bar" else "lines")
  }
  if (is.character(progress) && length(progress) == 1L && !is.na(progress) &&
      progress %in% c("bar", "lines")) {
    return(progress)
  }
  stop(sprintf('%s: `progress` must be TRUE, FALSE, "bar" or "lines"', what),
       call. = FALSE)
}

# The integer the search bindings take: 0 off, 1 the per-step summary on
# stderr, 2 per-step and per-candidate lines on stdout.
.ferx_search_progress <- function(progress, what) {
  switch(.ferx_progress_mode(progress, what),
         off = 0L, bar = 1L, lines = 2L)
}

# Write one progress line to stdout and push it out at once. Under `Rscript`
# stdout is a pipe and block-buffered, so without the flush a log viewer would
# see the lines in bursts, long after the fits they report.
.ferx_progress_line <- function(...) {
  cat(..., "\n", sep = "")
  flush(stdout())
  invisible(NULL)
}

# A number for a progress line: three decimals, or "-" when there is none.
.ferx_progress_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("-")
  formatC(x, format = "f", digits = 3)
}
