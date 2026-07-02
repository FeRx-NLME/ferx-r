# trailing "\r" stripped to tolerate Windows line endings.
# NOTE: "/.+\.csv$" assumes the path is the last token on the line. If
# ferx-core ever appends text after the path (e.g. "(appended)"), this
# match will fail silently and fall through to the fallback.
.ferx_find_trace_from_lines <- function(lines) {
  if (length(lines) == 0L) return(NULL)
  parse_one <- function(line) {
    line <- sub("\r$", "", line)
    m <- regexpr("/.+\\.csv$", line, perl = TRUE)
    if (m == -1L) return(NA_character_)
    regmatches(line, m)
  }
  # Primary: lines that look like the trace announcement.
  cand <- grep("optimizer trace", lines, value = TRUE, fixed = TRUE)
  if (length(cand) > 0L) {
    parsed <- vapply(cand, parse_one, character(1L), USE.NAMES = FALSE)
    parsed <- parsed[!is.na(parsed)]
    if (length(parsed) > 0L) return(parsed[[length(parsed)]])
  }
  # Fallback: any line that mentions "trace" and ends in ".csv".
  cand2 <- grep("trace", lines, value = TRUE, ignore.case = TRUE)
  if (length(cand2) == 0L) return(NULL)
  parsed <- vapply(cand2, parse_one, character(1L), USE.NAMES = FALSE)
  parsed <- parsed[!is.na(parsed)]
  if (length(parsed) == 0L) NULL else parsed[[length(parsed)]]
}
