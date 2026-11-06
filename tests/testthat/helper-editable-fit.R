# Fits are read-only (#713), and a test that exercises a damaged, old or
# otherwise impossible fit has to make one. These are the two ways to do that.

# `fit` without its seal, under a subclass whose assignment methods exist only
# here: every assignment works, and `ferx_save_fit()` refuses it. Anything but a
# fit passes through untouched.
editable_fit <- function(fit) {
  if (!inherits(fit, "ferx_fit")) return(fit)
  attr(fit, "ferx_seal") <- NULL
  class(fit) <- c("ferx_editable_fit", "ferx_fit")
  fit
}

local({
  unguarded <- function(x, assign) {
    cls <- oldClass(x)
    x <- assign(unclass(x))
    class(x) <- cls
    x
  }
  ns <- asNamespace("ferx")
  registerS3method("$<-", "ferx_editable_fit", function(x, name, value) {
    unguarded(x, function(l) { l[[name]] <- value; l })
  }, envir = ns)
  registerS3method("[[<-", "ferx_editable_fit", function(x, i, value) {
    unguarded(x, function(l) { l[[i]] <- value; l })
  }, envir = ns)
  registerS3method("[<-", "ferx_editable_fit", function(x, i, value) {
    unguarded(x, function(l) { l[i] <- value; l })
  }, envir = ns)
})

# An edited fit made sealed again, for a test of what the writer or a post-hoc
# step does with a fit no user can now produce: an older one, say.
resealed_fit <- function(fit) ferx:::.ferx_fit_seal(editable_fit(fit))
