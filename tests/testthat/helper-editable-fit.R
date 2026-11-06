# Fits are read-only (#713), and a test that exercises a damaged, old or
# otherwise impossible fit has to make one. These are the two ways to do that.

# `fit` without its seal: every assignment works, as on a hand-built list, and
# `ferx_save_fit()` refuses it.
editable_fit <- function(fit) {
  attr(fit, "ferx_seal") <- NULL
  fit
}

# An edited fit made sealed again, for a test of what the writer or a post-hoc
# step does with a fit no user can now produce: an older one, say.
resealed_fit <- function(fit) ferx:::.ferx_fit_seal(editable_fit(fit))
