# A fit's population fingerprint under another ferx version's encoding scheme
# (#559). `PopulationFingerprint::verify` cannot compare one, so it checks
# nothing and returns a single note saying so (ferx-core
# STALE_FINGERPRINT_WARNING); this is the only way a live fit produces that
# note. `scheme` is a plain field of the fingerprint's JSON, so an edit of it
# stays a valid fingerprint. Scheme 1 is long retired, never the current one.
foreign_scheme_fit <- function(fit) {
  fitrx_edit_fit(fit, function(w) {
    stopifnot(!identical(w$population_fingerprint$scheme, 1L))
    w$population_fingerprint$scheme <- 1L
    w
  })
}

# The note verify returns for it, word for word.
foreign_scheme_note <- paste0(
  "this fit's population fingerprint is from another ferx version's scheme, ",
  "so the population was not verified against the one the fit was given"
)
