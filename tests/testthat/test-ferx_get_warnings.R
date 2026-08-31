test_that("ferx_get_warnings(as_df = TRUE) returns the underlying data frame", {
  fit <- warfarin_fit_cov()
  df <- ferx_get_warnings(fit, as_df = TRUE)
  expect_identical(df, fit$warnings_structured)
})
test_that("ferx_get_warnings() prints a grouped summary", {
  fit <- warfarin_fit_cov()
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("ferx fit warnings", out)))
  # Footer tallies are always present
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})
test_that("ferx_get_warnings() rejects non-fit input", {
  expect_error(ferx_get_warnings(list()), "ferx_fit")
})
test_that("ferx_get_warnings() falls back to flat warnings when structured is absent", {
  fake <- structure(
    list(
      model_name = "legacy",
      warnings = c("Outer optimization did not converge", "something else"),
      warnings_structured = NULL
    ),
    class = "ferx_fit"
  )
  df <- ferx_get_warnings(fake, as_df = TRUE)
  expect_equal(nrow(df), 2L)
  expect_true(all(df$severity == "warning"))
  expect_true(all(df$category == "general"))
})
test_that(".ferx_warning_guidance dispatches the covariance family by message content", {
  # Each message is paired with the category ferx-core's `classify_warning`
  # actually assigns it, not with a hand-picked one. Every message below carries
  # "covariance step failed" or "covariance step" + "not positive definite", so
  # core codes it `covariance_failed`; the regularisation messages carry
  # "covariance step regularized", so core codes those `covariance_regularized`.
  # Neither reaches the code `covariance_step`, which core reserves for its
  # Info-level cost note -- so pinning these to "covariance_step" (as this test
  # did until the routing fix) asserted a pairing production never emits, and
  # every branch under test was dead in the real call path.
  g <- function(msg, category = "covariance_failed") {
    ferx:::.ferx_warning_guidance(category, message = msg)
  }

  # NonPdHessian path: eigenvalue list in message.
  msg_npd <- paste0(
    "Covariance step: Hessian is not positive definite. ",
    "Eigenvalues: [8.4000, 2.1000, -0.0100]. SE estimates not available."
  )
  expect_match(g(msg_npd), "eigenvalue", ignore.case = TRUE)
  expect_match(g(msg_npd), "near-zero|negative", ignore.case = TRUE, perl = TRUE)

  # Ill-conditioned Hessian entries: names a parameter.
  # Every parameter in this message carries one of ferx-core's two cause
  # labels; "non-finite diagonal" was never one of them.
  msg_ic <- paste0(
    "Covariance step failed: Hessian has ill-conditioned entries for the ",
    "following parameter(s) \u2014 theta[CL] (FD stencil non-finite; model may ",
    "overflow at perturbation \u2014 try tuning fd_hessian_step). ",
    "SE estimates not available."
  )
  expect_match(g(msg_ic), "fd_hessian_step", ignore.case = TRUE)

  # Omega non-PD.
  # A tiny POSITIVE minimum eigenvalue: ferx-core writes "near-singular" for
  # this case, not "not positive definite" (covariance.rs picks the descriptor
  # from the sign), so the fixture must too.
  msg_omega <- paste0(
    "Covariance step failed: Omega matrix is not positive definite at ",
    "convergence (min eigenvalue = 1.2e-10; eigenvalues: [0.5000, 1.2e-10]). ",
    "SE estimates not available."
  )
  expect_match(g(msg_omega), "near-singular", ignore.case = TRUE)
  expect_false(grepl("eigenvalue list", g(msg_omega), ignore.case = TRUE))

  # Non-finite OFV.
  msg_ofv <- paste0(
    "Covariance step failed: base OFV is non-finite at convergence ",
    "(likely numerical overflow or underflow in model evaluation). ",
    "SE estimates not available."
  )
  expect_match(g(msg_ofv), "overflow|underflow", ignore.case = TRUE, perl = TRUE)

  # Regularisation: minor, moderate, severe.
  base_reg <- function(sev) paste0(
    "Covariance step regularized: eigenvalue floor applied to FD Hessian ",
    "(1 of 3 free-block eigenvalues clipped; min eig = 1.2e-6, floor = 8.4e-14; ",
    "severity: ", sev, "). Standard errors are likely reliable."
  )
  reg <- function(sev) g(base_reg(sev), category = "covariance_regularized")
  # Anchored on the tier word itself. The previous alternations did not
  # discriminate: "unreliable" contains "reliable", so the severe text matched
  # the minor pattern, and the moderate text matched the severe pattern --
  # swapping the moderate and severe return values kept the suite green.
  expect_match(reg("minor"),    "Minor Hessian regularisation",    fixed = TRUE)
  expect_match(reg("moderate"), "Moderate Hessian regularisation", fixed = TRUE)
  expect_match(reg("severe"),   "Severe Hessian regularisation",   fixed = TRUE)
  # Each tier's own word, and nobody else's.
  expect_false(grepl("Moderate|Severe", reg("minor")))
  expect_false(grepl("Minor|Severe",    reg("moderate")))
  expect_false(grepl("Minor|Moderate",  reg("severe")))
  # Only minor is benign enough to omit the SIR cross-check.
  expect_false(grepl("ferx_sir", reg("minor"), ignore.case = TRUE))
  expect_match(reg("moderate"), "ferx_sir", fixed = TRUE)
  expect_match(reg("severe"),   "ferx_sir", fixed = TRUE)
  # The tier is chosen by what FRACTION of the free-block eigenvalues had to be
  # clipped (covariance.rs), not by the size of the floor, and core's own
  # minor-tier interpretation is just "standard errors are likely reliable" --
  # so the guidance must not volunteer how common or how benign it is.

  # Generic fallback for unrecognised message.
  expect_match(g("Covariance step failed"), "identifiability", ignore.case = TRUE)

})

test_that("the informational covariance_step note does not read as a failure", {
  # Core codes only its Info-level cost note `covariance_step`. It matches none
  # of the targeted branches, so before this arm existed it fell through to the
  # "standard errors unavailable" fallback and told the user the step had failed
  # when it had not yet run.
  # ferx-core's primary cost note, verbatim (api/fit.rs).
  msg <- paste0("Covariance step: 35 parameters \u2192 1225 OFV evaluations ",
                "(finite-difference Hessian). This may take several minutes on ",
                "complex models.")
  g <- ferx:::.ferx_warning_guidance("covariance_step", message = msg)
  expect_match(g, "Informational", ignore.case = TRUE)
  expect_false(grepl("unavailable|identifiability", g, ignore.case = TRUE))
})

test_that("an unrecognised covariance_step message is not called benign", {
  # The informational arm is matched on the message, not on the code alone.
  # Core's `CovarianceStep` classification ends in a catch-all -- any otherwise
  # unmatched message carrying "covariance step:" and "parameters" lands there
  # -- so the code is not proof the message is the cost note. Answering an
  # unknown one "no action needed" fails in the dangerous direction, so it must
  # fall through to the failure fallback instead.
  g <- ferx:::.ferx_warning_guidance(
    "covariance_step",
    message = "Covariance step: something new about 12 parameters."
  )
  # Not the benign cost note (the message is not the cost note) and not the
  # failure fallback either: `covariance_step` is Info by construction, so
  # "standard errors unavailable" would contradict the row it prints under.
  # Not the benign cost note ...
  # Not the benign cost note: that arm matches on the message, never on the
  # code alone, because `covariance_step`'s classification ends in a catch-all.
  # An unrecognised message under it falls through to the generic fallback.
  expect_false(grepl("No action needed", g, fixed = TRUE))
  expect_match(g, "identifiability", ignore.case = TRUE)
})

test_that("a covariance step that succeeded is not reported as unavailable", {
  # ferx-core emits this on the CovarianceStepResult::Success path: the matrix
  # was produced and the SEs exist, they are merely missing their cross-partial
  # terms. `classify_warning` codes it `covariance_regularized` (the
  # "off-diagonal fd stencil" arm), so the routing fix admits it to the
  # covariance block for the first time -- where, before this arm, it matched
  # nothing and inherited "Standard errors unavailable", contradicting the
  # message it prints under.
  msg <- paste0(
    "Covariance step: off-diagonal FD stencil(s) non-finite for theta[CL], ",
    "sigma[1]. Cross-partial correlation set to 0; SE for these parameter(s) ",
    "may be over-optimistic. Try tuning fd_hessian_step."
  )
  g <- ferx:::.ferx_warning_guidance("covariance_regularized", message = msg)
  expect_match(g, "over-optimistic", ignore.case = TRUE)
  expect_match(g, "fd_hessian_step", fixed = TRUE)
  expect_false(grepl("unavailable", g, ignore.case = TRUE))
})

test_that("ferx_covariance()'s `covariance` category reaches the guidance", {
  # `ferx_covariance()` folds the engine's flat covariance warnings into the
  # structured table under `category = "covariance"` (R/ferx_covariance.R) --
  # a fourth token, not one of core's three. Omitting it left the entire
  # post-hoc covariance surface without guidance even after the routing fix.
  msg <- paste0(
    "Covariance step failed: Hessian has ill-conditioned entries for the ",
    "following parameter(s) \u2014 theta[CL] (zero diagonal \u2014 flat objective)."
  )
  g <- ferx:::.ferx_warning_guidance("covariance", message = msg)
  expect_match(g, "Hessian diagonal", fixed = TRUE)
})

test_that("a legacy fit's `general` covariance message still gets guidance", {
  # `ferx_get_warnings()` stamps `general` on every message of a fit saved
  # before structured warnings existed. Every branch is keyed on the message,
  # so such a message can be answered exactly as well as one arriving under its
  # own code -- routing it around the block re-created, for those fits, the
  # silence this whole change removes.
  msg <- paste0(
    "Covariance step failed: Omega matrix is not positive definite at ",
    "convergence."
  )
  g <- ferx:::.ferx_warning_guidance("general", message = msg)
  expect_match(g, "Omega is near-singular", fixed = TRUE)
  # A `general` message that is not about the covariance step stays unhandled.
  expect_null(ferx:::.ferx_warning_guidance("general", message = "Something else"))
})
test_that("ferx_get_warnings() shows guidance for an unused declaration", {
  # End to end through the printer, with the category ferx-core actually
  # assigns: its unused-declaration messages classify to `general`, not to
  # `unused_parameter`. Pinning this to the hand-made pairing (as it did until
  # the routing fix) asserted a row production never produces.
  msg <- paste0(
    "theta 'TVCL' is declared in [parameters] but not referenced in any model ",
    "expression -- it will not affect predictions or be meaningfully estimated"
  )
  fake <- structure(
    list(
      model_name = "m",
      warnings_structured = data.frame(
        severity      = "warning",
        category      = "general",
        message       = msg,
        source_method = "",
        stringsAsFactors = FALSE
      ),
      warnings = msg,
      condition_number = NULL,
      eta_normality = NULL,
      uses_sde = FALSE
    ),
    class = "ferx_fit"
  )
  # Collapsed and whitespace-normalised before matching: the printer wraps at
  # 70 columns, so asserting on a phrase in a single output line passes or fails
  # on where the wrap happens to fall rather than on the behaviour under test.
  out <- capture.output(ferx_get_warnings(fake))
  flat <- gsub("[[:space:]]+", " ", paste(out, collapse = " "))
  expect_true(
    grepl("Remove it from", flat, fixed = TRUE),
    info = paste("Expected guidance text not found in output:\n",
                 paste(out, collapse = "\n"))
  )
})

# ---- header from test-diagnostics-helpers.R ----
# Unit tests for the pure-R diagnostic/formatting helpers in diagnostics.R:
# the warning-guidance lookup and the cli styling shims. No model fit required.

# ---------------------------------------------------------------------------
# .ferx_warning_guidance
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# .ferx_use_cli / .ferx_style
# ---------------------------------------------------------------------------




test_that(".ferx_warning_guidance toggles the SDE hint for positive autocorrelation", {
  pos_ode <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive", uses_sde = FALSE)
  expect_match(pos_ode, "Positive IWRES autocorrelation")
  expect_match(pos_ode, "SDE process noise")            # hint shown when not already using SDE
  pos_sde <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive", uses_sde = TRUE)
  expect_match(pos_sde, "Positive IWRES autocorrelation")
  expect_false(grepl("SDE process noise", pos_sde))     # hint suppressed
})
test_that(".ferx_warning_guidance returns negative-autocorrelation guidance", {
  neg <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "Negative autocorrelation")
  expect_match(neg, "Negative IWRES autocorrelation")
})
# Every token ferx-core's `WarningCode::as_str()` can emit, minus `general`.
# Hand-transcribed, and deliberately so: `.core_warning_cats()` below
# checks it against the real thing whenever a sibling ferx-core checkout is
# present, which is the only mechanism here that can notice core growing a code.
.core_warning_cats <- function() {
  c(
    "absorption_twin_declined", "bloq_method", "boundary_estimate", "cancelled",
    "condition_number", "convergence", "covariance_failed",
    "covariance_regularized", "covariance_step", "data_quality",
    "dw_autocorrelation", "eps_shrinkage", "eta_normality", "eta_shrinkage",
    "experimental", "flat_parameter", "flip_flop", "gradient_fallback",
    "high_correlation", "importance_sampling", "inflated_rse", "mu_referencing",
    "multi_start", "ode_solver", "omega_structure", "optimizer_config",
    "optimizer_health", "parameter_at_runaway_guard", "simulation", "sir",
    "threads", "vi_bad_basin"
  )
}

# Categories that reach `.ferx_warning_guidance()` without being ferx-core
# `WarningCode` tokens.
#
#   covariance        - assigned by ferx_covariance() (R/ferx_covariance.R) to
#                       the engine's flat covariance warnings.
#   unused_parameter  - no core `WarningCode`; both of core's unused-declaration
#                       messages classify to `general` and are rerouted here by
#                       message, so the arm is reachable (see the test below
#                       that drives it with core's real message texts).
#
# `ebe_convergence` used to sit here too and has been removed: nothing emits it.
# Core has no such `WarningCode`, and ferx-r reports EBE convergence as the
# integer field `ebe_convergence_warnings`, never as a warning row. Keeping a
# guidance arm for it in a change whose premise is that unreachable branches are
# a defect would have been the same mistake in miniature.
.extra_warning_cats <- function() {
  c("covariance", "unused_parameter")
}

# ferx-core tokens this package deliberately answers with nothing. Guidance text
# for them is a separate change; what matters here is that the set is explicit,
# so a code ferx-core adds later lands in neither list and fails the drift test
# rather than silently printing nothing.
.unanswered_warning_cats <- function() {
  c("absorption_twin_declined", "boundary_estimate", "eps_shrinkage",
    "eta_shrinkage", "experimental", "flat_parameter", "flip_flop",
    "high_correlation", "inflated_rse", "simulation")
}

test_that("every category in the guidance table returns guidance", {
  # A completeness walk, and nothing more. Two things it deliberately does NOT
  # establish, both covered elsewhere:
  #   - it is not a drift guard. The loop visits only the tokens listed above,
  #     so a code core added and nobody transcribed is invisible to it; that is
  #     what `matches ferx-core's WarningCode vocabulary` is for.
  #   - it says nothing about the covariance family. Called with an empty
  #     message, all four covariance tokens land on the generic fallback, so
  #     this would still pass if `.ferx_covariance_guidance()` were reduced to
  #     its last line. The message inventory below is what pins those branches.
  answered <- setdiff(c(.core_warning_cats(), .extra_warning_cats()),
                      .unanswered_warning_cats())
  for (cat in answered) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
  # ... and the deliberately unanswered ones really are unanswered, so the list
  # cannot rot into a place where handled codes hide.
  for (cat in .unanswered_warning_cats()) {
    expect_null(ferx:::.ferx_warning_guidance(cat), info = cat)
  }
})

# Locate a sibling ferx-core checkout. Walks upward rather than using a fixed
# relative depth: testthat's wd is tests/testthat, so `../../../ferx-core` only
# resolves for a plain checkout. ferx-r's CLAUDE.md mandates working in a git
# worktree under <repo>/.claude/worktrees/<name>, where that path points at
# <repo>/.claude/worktrees/ferx-core and never exists -- so the guards that
# depend on it silently skipped in exactly the workflow they were written for.
.find_core_src <- function() {
  d <- normalizePath(".", mustWork = FALSE)
  for (i in seq_len(10)) {
    cand <- file.path(d, "ferx-core", "src")
    if (dir.exists(cand)) return(cand)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NULL
}

# Materialise ferx-core's `src/` at the revision this package is PINNED to, and
# return that directory (or NULL).
#
# The sibling checkout's working tree is NOT the right source to check against:
# `Cargo.toml` tracks `branch = "main"`, but CI builds the revision recorded in
# `src/rust/Cargo.lock`, and a local checkout can sit either side of it. Reading
# the pin means these guards assert against the vocabulary and the message texts
# that actually ship. Falls back to NULL (skip) rather than silently reading a
# different revision.
.core_src_at_pin <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(if (identical(cached, NA_character_)) NULL else cached)
    fail <- function() { cached <<- NA_character_; NULL }
    core_src <- .find_core_src(); if (is.null(core_src)) return(fail())
    repo <- dirname(core_src)
    lock <- file.path("..", "..", "src", "rust", "Cargo.lock")
    if (!file.exists(lock)) return(fail())
    ln <- readLines(lock, warn = FALSE)
    i <- grep('^name = "ferx-core"$', ln)
    if (!length(i)) return(fail())
    src_line <- grep("ferx-core.*#", ln[i[1]:min(length(ln), i[1] + 5L)], value = TRUE)
    if (!length(src_line)) return(fail())
    rev <- sub('".*$', "", sub("^.*#", "", src_line[1]))
    if (!grepl("^[0-9a-f]{7,40}$", rev)) return(fail())
    have <- suppressWarnings(system2("git", c("-C", repo, "cat-file", "-e",
                                              paste0(rev, "^{commit}")),
                                     stdout = FALSE, stderr = FALSE))
    if (!identical(have, 0L)) return(fail())
    dest <- file.path(tempdir(), paste0("ferx-core-", substr(rev, 1, 12)))
    if (!dir.exists(dest)) {
      dir.create(dest, recursive = TRUE)
      tar <- file.path(tempdir(), paste0(substr(rev, 1, 12), ".tar"))
      ok <- suppressWarnings(system2("git", c("-C", repo, "archive", "--format=tar",
                                              "-o", tar, rev, "src"),
                                     stdout = FALSE, stderr = FALSE))
      if (!identical(ok, 0L)) return(fail())
      untar(tar, exdir = dest)
    }
    out <- file.path(dest, "src")
    if (!dir.exists(out)) return(fail())
    cached <<- out
    out
  }
})

test_that(".ferx_warning_guidance matches ferx-core's WarningCode vocabulary", {
  # The real drift guard. Reads `WarningCode::as_str()` out of a sibling
  # ferx-core checkout instead of trusting a second hand-maintained copy, so a
  # code added in core fails here even though nobody edited this file.
  #
  # Skipped when the sibling is absent: ferx-r CI builds the pinned crate, not a
  # checkout. That is the right trade -- drift is introduced on a developer
  # machine, which per CLAUDE.md always has ../ferx-core, and this fails there.
  core_src <- .core_src_at_pin()
  skip_if(is.null(core_src), "cannot materialise ferx-core at the pinned revision")
  types_rs <- file.path(core_src, "types.rs")
  skip_if(!file.exists(types_rs), "sibling ferx-core has no src/types.rs")

  src <- readLines(types_rs, warn = FALSE)
  hits <- regmatches(src, regexpr('WarningCode::[A-Za-z]+ => "[a-z_]+"', src))
  tokens <- gsub('.*=> "|"$', "", hits)
  expect_gt(length(tokens), 20L)          # the file was found and understood

  # `general` is core's unrecognised-message bucket and is deliberately absent
  # from the table -- the message text is the only guidance there is.
  expect_setequal(setdiff(tokens, "general"), .core_warning_cats())
  for (cat in setdiff(tokens, c("general", .unanswered_warning_cats()))) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
})

test_that(".ferx_warning_guidance gives `general` no category-level guidance", {
  # `general` is core's bucket for a message its classifier did not recognise,
  # so there is no category-level remediation to give -- the message text is the
  # guidance. Asserted so it reads as a decision rather than an omission. The
  # one exception is a covariance-step message arriving under `general`, which
  # the message-keyed covariance family still answers (tested above).
  expect_null(ferx:::.ferx_warning_guidance("general"))
  expect_null(ferx:::.ferx_warning_guidance("general", message = "Anything at all."))
})
test_that("core's unused-declaration messages reach the unused_parameter guidance", {
  # ferx-core has no `unused_parameter` WarningCode -- both messages fall
  # through `classify_warning` to `general`, so the arm was unreachable in
  # production for exactly the reason the covariance block was. The pre-existing
  # test passed only because it fed the category in by hand.
  theta_msg <- paste0(
    "theta 'TVCL' is declared in [parameters] but not referenced in any model ",
    "expression -- it will not affect predictions or be meaningfully estimated"
  )
  ip_msg <- paste0(
    "[individual_parameters] `KA` is computed but never used \u2014 not mapped ",
    "into the `pk(...)` model and not referenced in any other block, so it has ",
    "no effect. Map `KA` in [structural_model] (e.g. `f=F`) or remove `KA`."
  )
  for (msg in c(theta_msg, ip_msg)) {
    g <- ferx:::.ferx_warning_guidance("general", message = msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = msg)
    expect_match(g, "never referenced", fixed = TRUE, info = msg)
  }
  # An unrelated `general` message is still left to its own text.
  expect_null(ferx:::.ferx_warning_guidance("general", message = "Something else."))
})

test_that("every covariance message ferx-core emits gets non-contradictory guidance", {
  # The inventory test. Every message below is copied verbatim from a ferx-core
  # emit site and paired with the category `classify_warning` assigns it -- not
  # with a category chosen by hand. The previous tests enumerated only messages
  # the block was already known to handle, which is why a whole dead code path
  # had green tests, and why the routing fix could newly admit messages the
  # block answers wrongly without anything failing.
  #
  # `ok` is what the guidance must say; `never` is what it must not. The
  # standing invariant is the `never` column on the two success-path messages:
  # guidance must not tell a user their standard errors are unavailable when the
  # engine produced them.
  # Each `msg` is ferx-core's own text with its `{}` placeholders filled in --
  # em dashes, arrows and the superscript are the real characters, not ASCII
  # stand-ins. `ok` is what the guidance must say; `never` is what it must not.
  # The standing invariant is the `never` column on the success-path messages:
  # guidance must not report standard errors unavailable when the engine
  # produced them.
  # Each `msg` is ferx-core's own text with its `{}` placeholders filled in.
  # `ok` is what the guidance must say; `never` what it must not. The standing
  # invariant is the `never` column on the success-path messages: guidance must
  # not report standard errors unavailable when the engine produced them.
  cases <- list(
    list(cat = "covariance_failed",
         msg = paste0("Covariance step: Hessian is not positive definite. ",
                      "Eigenvalues: [-0.0012, 0.4500]. SE estimates not available."),
         ok = "eigenvalue list", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: Hessian has ill-conditioned entries ",
                      "for the following parameter(s) \u2014 theta[CL] (FD stencil ",
                      "non-finite; model may overflow at perturbation \u2014 try ",
                      "tuning fd_hessian_step). SE estimates not available."),
         ok = "Hessian diagonal", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: Omega matrix is not positive definite ",
                      "at convergence (min eigenvalue = -1.000e-9; eigenvalues: ",
                      "[-1.000e-9, 0.0900]). SE estimates not available."),
         ok = "Omega is near-singular", never = NULL),
    list(cat = "covariance_failed",
         # Same emit site, the other descriptor. ferx-core picks it from the
         # sign of the smallest eigenvalue, so both must reach the omega arm
         # rather than one of them falling to the generic fallback.
         msg = paste0("Covariance step failed: Omega matrix is near-singular at ",
                      "convergence (min eigenvalue = 1.000e-12; eigenvalues: ",
                      "[1.000e-12, 0.0900]). SE estimates not available."),
         ok = "Omega is near-singular", never = "identifiability"),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: base OFV is non-finite at convergence ",
                      "(likely numerical overflow or underflow in model evaluation). ",
                      "SE estimates not available."),
         ok = "overflow", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: could not compute eigenvalues of the ",
                      "FD Hessian (Hessian may contain NaN or Inf). SE estimates ",
                      "not available."),
         ok = "identifiability", never = NULL),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (1 of 3 free-block eigenvalues clipped; min eig = ",
                      "1.200e-6, floor = 8.400e-14; severity: minor). Standard ",
                      "errors are likely reliable."),
         ok = "Minor Hessian regularisation", never = "ferx_sir"),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (2 of 3 free-block eigenvalues clipped; min eig = ",
                      "1.200e-8, floor = 8.400e-14; severity: moderate). Standard ",
                      "errors should be interpreted with caution."),
         ok = "Moderate Hessian regularisation", never = "Severe"),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (3 of 3 free-block eigenvalues clipped; min eig = ",
                      "-1.200e-8, floor = 8.400e-14; severity: severe). SIR-based ",
                      "confidence intervals are recommended."),
         ok = "Severe Hessian regularisation", never = "Moderate"),
    list(cat = "covariance_regularized",
         # SUCCESS path -- the covariance matrix exists.
         msg = paste0("Covariance step: off-diagonal FD stencil(s) non-finite for ",
                      "theta[CL], sigma[1]. Cross-partial correlation set to 0; SE ",
                      "for these parameter(s) may be over-optimistic. Try tuning ",
                      "fd_hessian_step."),
         ok = "over-optimistic", never = "unavailable"),
    list(cat = "covariance_step",
         # Info-level cost note, primary form -- the step has not even run.
         msg = paste0("Covariance step: 35 parameters \u2192 1225 OFV evaluations ",
                      "(finite-difference Hessian). This may take several minutes ",
                      "on complex models."),
         ok = "Informational", never = "unavailable"),
    list(cat = "covariance_step",
         # Second cost-note form: the evaluation count overflows usize.
         msg = paste0("Covariance step: 4294967296 parameters \u2192 n\u00b2 OFV ",
                      "evaluations (finite-difference Hessian). Estimate exceeds ",
                      "usize range; expect this to be very slow."),
         ok = "Informational", never = "unavailable"),
    list(cat = "general",
         # COV_CANCELLED_MSG, which classifies to `general` and so reaches the
         # block only because `general` is admitted on the message.
         msg = paste0("Covariance step cancelled before completion; standard ",
                      "errors not available."),
         ok = "nothing was diagnosed", never = "identifiability")
  )
  for (case in cases) {
    g <- ferx:::.ferx_warning_guidance(case$cat, message = case$msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = case$msg)
    expect_match(g, case$ok, ignore.case = TRUE, info = case$msg)
    if (!is.null(case$never)) {
      expect_false(grepl(case$never, g, ignore.case = TRUE), info = case$msg)
    }
    # Same message under ferx_covariance()'s token must give the same answer --
    # true for every case here because all of them hit a targeted branch. It is
    # deliberately NOT true in the fallback, where the code carries information
    # the message does not; see "the fallback does not report failure under a
    # code that means success".
    expect_identical(ferx:::.ferx_warning_guidance("covariance", message = case$msg),
                     ferx:::.ferx_warning_guidance(case$cat, message = case$msg),
                     info = case$msg)
  }
})


test_that("a flat THETA is not answered as an unused parameter, by either door", {
  # ferx-core's flat-theta message contains the literal phrase "computed but
  # never used", so the unused-parameter reroute captured it. Gating the reroute
  # on `general` did NOT fix that: `ferx_load_fit()` does not restore
  # `warnings_structured`, so on a loaded fit every row arrives under `general`
  # and the gate becomes the admission criterion rather than an exclusion. The
  # message must be answered the same way through both doors.
  msg <- paste0("[parameters] `TVCL` has no effect on the objective (gradient ",
                "\u2248 0 at the initial estimate) \u2014 it is likely computed but ",
                "never used (unmapped, or dropped from the structural / scaling ",
                "model). Freezing it at its initial value (1.5) so the remaining ",
                "parameters can be estimated; map or remove `TVCL` to silence this.")
  # This package has no `flat_parameter` guidance (that text is a separate
  # change), so the requirement here is only that the message is NOT captured by
  # the unused-parameter reroute and answered with its text.
  loaded <- ferx:::.ferx_warning_guidance("general", message = msg)
  expect_null(loaded)
  expect_false(identical(loaded, ferx:::.ferx_warning_guidance("unused_parameter")))

  # The genuine unused-declaration messages still resolve, and a live fit's own
  # typed code still wins over anything the text might suggest.
  unused <- paste0("theta 'TVX' is declared in [parameters] but not referenced in ",
                   "any model expression \u2014 it will not affect predictions")
  expect_match(ferx:::.ferx_warning_guidance("general", message = unused),
               "never referenced", fixed = TRUE)
  expect_null(ferx:::.ferx_warning_guidance("flat_parameter", message = unused))
})

test_that("`general` does not claim messages that only mention the covariance step", {
  # `ferx_load_fit()` does not restore `warnings_structured`, so EVERY row of a
  # fit read back from disk arrives under `general`. Admitting any message that
  # merely contains "covariance step" therefore handed SIR's own diagnostics the
  # covariance fallback -- which tells the user to pass covariance = FALSE, the
  # one setting that removes the matrix SIR needs.
  # Verbatim from sir.rs and postfit.rs -- both clauses these previously dropped
  # ("so draws mostly stay inside the parameter bounds", "see the covariance
  # warning above for the cause") are the parts that mention the covariance
  # step, so trimming them weakened the very thing under test.
  shrunk <- paste0("SIR: proposal was shrunk in 1 direction(s) so draws mostly stay ",
                   "inside the parameter bounds [KA +1.00]. Those directions come ",
                   "from eigenvalue-floored (non-identified) curvature in the ",
                   "covariance step; the SIR CIs along them understate the true ",
                   "uncertainty.")
  requested <- paste0("SIR requested but the covariance step did not succeed and no ",
                      "usable SIR proposal could be built from it, so SIR could not ",
                      "run \u2014 see the covariance warning above for the cause.")
  for (msg in c(shrunk, requested)) {
    expect_null(ferx:::.ferx_warning_guidance("general", message = msg), info = msg)
  }
  # ... while a real covariance message under `general` is still answered, with
  # or without a [METHOD] chain prefix.
  cancelled <- "Covariance step cancelled before completion; standard errors not available."
  expect_match(ferx:::.ferx_warning_guidance("general", message = cancelled),
               "cancelled", ignore.case = TRUE)
  expect_match(ferx:::.ferx_warning_guidance("general", message = paste0("[FOCEI] ", cancelled)),
               "cancelled", ignore.case = TRUE)
})

test_that("no covariance message ferx-core emits is missing from the inventory", {
  # The completeness guard for the inventory test above.
  #
  # That test is only as good as the list of messages someone hand-collected for
  # it, and a hand-collected list cannot notice a message added later. This one
  # derives the list from ferx-core's source instead: every string literal that
  # *begins* with "Covariance step" -- which is the shape of every covariance
  # warning the engine emits -- must be accounted for, either by a case in the
  # inventory or by an explicit exemption naming why it is not a warning.
  #
  # Skipped without a sibling ferx-core checkout, for the same reason as the
  # WarningCode vocabulary test: ferx-r CI builds the pinned crate rather than a
  # checkout, and drift is introduced on a developer machine.
  core_src <- .core_src_at_pin()
  skip_if(is.null(core_src), "cannot materialise ferx-core at the pinned revision")

  files <- list.files(core_src, pattern = "[.]rs$",
                      recursive = TRUE, full.names = TRUE)
  # Test-only sources: sibling `*_tests.rs` files and the `src/api/tests/` dir.
  files <- files[!grepl("_tests[.]rs$|/tests/", files)]
  lits <- character(0)
  for (f in files) {
    ln <- readLines(f, warn = FALSE)
    # Drop whole-line `//` comments before joining: ferx-core quotes its own
    # message text inside comments, and those quotes are not emit sites.
    ln <- ln[!grepl("^[[:space:]]*//", ln)]
    # Rust continues a string literal with a trailing backslash; rejoin so a
    # wrapped message is recovered as one literal.
    joined <- gsub("\\\\\n[[:space:]]*", "", paste(ln, collapse = "\n"))
    m <- gregexpr('"(?:[^"\\\\]|\\\\.)*"', joined, perl = TRUE)
    for (s in regmatches(joined, m)[[1]]) {
      body <- substr(s, 2, nchar(s) - 1)
      if (grepl("^Covariance step", body)) lits <- c(lits, body)
    }
  }
  lits <- sort(unique(lits))
  expect_gt(length(lits), 8L)          # the source was found and understood

  # Not warnings, so no guidance is owed. Keyed on a distinctive fragment.
  exempt <- c(
    # Inline `#[cfg(test)]` fixture in run_sir.rs (the file itself is not a
    # test-only source, so the path filter above cannot drop it).
    "matrix was not positive definite",
    # Answered by the generic fallback rather than a targeted arm: correct as
    # far as it goes (the standard errors really are unavailable), just not
    # specific. Targeted text for these is a separate change.
    "fd_hessian_step must be",
    "score cross-product matrix"
  )
  # Fragments the inventory test covers. Not one per case: the two Omega
  # descriptors share "Omega matrix is" and the two cost-note forms share
  # "OFV evaluations", so 11 fragments cover 15 cases.
  covered <- c(
    "cancelled before completion",
    "ill-conditioned entries",
    "Omega matrix is",
    "base OFV is non-finite",
    "could not compute eigenvalues",
    "regularized: eigenvalue floor",
    "Hessian is not positive definite",
    "off-diagonal FD stencil",
    "OFV evaluations"
  )
  accounted <- function(lit) any(vapply(c(exempt, covered),
                                        function(k) grepl(k, lit, fixed = TRUE),
                                        logical(1)))
  unaccounted <- lits[!vapply(lits, accounted, logical(1))]
  expect_identical(
    unaccounted, character(0),
    info = paste0(
      "ferx-core emits covariance message(s) with no case in the inventory test ",
      "and no exemption. Add a case to \"every covariance message ferx-core ",
      "emits gets non-contradictory guidance\", or exempt it here with the ",
      "reason it is not a warning:\n  ",
      paste(unaccounted, collapse = "\n  ")
    )
  )
})

# ---------------------------------------------------------------------------
# End-to-end anchor: the engine's own warnings, through the real classifier
# ---------------------------------------------------------------------------
# Every other test here feeds `.ferx_warning_guidance()` a message typed by
# hand. That exercises the dispatch but not the seam that actually matters --
# ferx-core emitting a message, `classify_warning` assigning it a category, and
# this package answering the pair. Nothing pinned that seam, and it is exactly
# where the routing defect this file exists to fix lived: the categories were
# hand-paired in the tests, so a pairing production never emits looked correct.
#
# The model is deliberately pathological: one theta declared and never
# referenced makes ferx-core emit both a real unused-declaration warning and a
# real covariance-step failure, on a fit that settles in a couple of seconds.
warning_anchor_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      mf <- tempfile(fileext = ".ferx")
      writeLines(c(
        "[parameters]",
        "  theta TVCL(0.134, 0.001, 10.0)",
        "  theta TVV(8.1, 0.1, 500.0)",
        "  theta TVKA(1.0, 0.01, 50.0)",
        "  theta TVUNUSED(1.0, 0.1, 10.0)",
        "  omega ETA_CL ~ 0.07",
        "  omega ETA_V  ~ 0.02",
        "  omega ETA_KA ~ 0.40",
        "  sigma PROP_ERR ~ 0.01 (sd)",
        "",
        "[individual_parameters]",
        "  CL = TVCL * exp(ETA_CL)",
        "  V  = TVV  * exp(ETA_V)",
        "  KA = TVKA * exp(ETA_KA)",
        "",
        "[structural_model]",
        "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
        "",
        "[error_model]",
        "  DV ~ proportional(PROP_ERR)"
      ), mf)
      ex <- ferx_example("warfarin")
      fit <<- ferx_fit(mf, ex$data, method = "gn", verbose = FALSE,
                       covariance = TRUE, settings = list(maxiter = 30L))
    }
    fit
  }
})

test_that("the engine's own warnings reach the guidance they were routed for", {
  ws <- warning_anchor_fit()$warnings_structured
  expect_true(is.data.frame(ws) && nrow(ws) > 0L)
  guide <- function(i) ferx:::.ferx_warning_guidance(ws$category[i], message = ws$message[i])

  # 1. The unused theta. ferx-core has no `unused_parameter` code, so this
  #    arrives under `general` and is recovered from the message. If core gains
  #    a code for it, or rewords it, this fails instead of silently printing
  #    nothing.
  i_unused <- which(grepl("TVUNUSED", ws$message, fixed = TRUE) &
                    !grepl("^Covariance step", ws$message))
  expect_length(i_unused, 1L)
  expect_identical(ws$category[i_unused], "general")
  expect_match(guide(i_unused), "never referenced", fixed = TRUE)

  # 2. The covariance failure: a CRITICAL row that printed no guidance at all
  #    before the routing fix, because the block was gated on `covariance_step`
  #    while `classify_warning` codes this `covariance_failed`.
  i_cov <- which(ws$category == "covariance_failed")
  expect_length(i_cov, 1L)
  expect_match(ws$message[i_cov], "ill-conditioned entries", fixed = TRUE)
  # The engine labels TVUNUSED "(zero diagonal -- flat objective)", so the
  # guidance must answer THAT cause -- not the finite-difference one, which
  # would be wrong here twice over: the label says the stencil succeeded, and
  # this fit ran on the analytic R-matrix path where nothing is differenced.
  g_cov <- guide(i_cov)
  expect_match(ws$message[i_cov], "zero diagonal", fixed = TRUE)
  expect_match(g_cov, "Hessian diagonal", fixed = TRUE)
  # ... and specifically not the generic fallback, which is where a message the
  # targeted branches failed to recognise would land.
  expect_false(grepl("Check identifiability", g_cov, fixed = TRUE))

  # 3. No row in a real table may produce an empty guidance block.
  for (i in seq_len(nrow(ws))) {
    g <- guide(i)
    expect_true(is.null(g) || (is.character(g) && nzchar(g)), info = ws$message[i])
  }
})

test_that(".ferx_warning_guidance returns NULL for an unknown category", {
  expect_null(ferx:::.ferx_warning_guidance("not_a_real_category"))
})

# ---- header from test-diagnostics-more.R ----
# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   ferx_estimates(), .ferx_est_row(), ferx_eta_cov(), ferx_cor_matrix(),
#   ferx_get_warnings(), and .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row       <- getFromNamespace(".ferx_est_row",            "ferx")
.compute_norm  <- getFromNamespace(".ferx_compute_eta_normality", "ferx")

# ---------------------------------------------------------------------------
# ferx_estimates — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_eta_cov — message branches plus the correlation path
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# ferx_cor_matrix — non-positive diagonal warning
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_get_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------




test_that("ferx_get_warnings prints 'No warnings.' for an empty structured table", {
  empty <- data.frame(severity = character(0), category = character(0),
                      message = character(0), source_method = character(0),
                      stringsAsFactors = FALSE)
  fit <- make_fake_fit(warnings_structured = empty)
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("No warnings", out)))
})
test_that("ferx_get_warnings labels each severity level", {
  df <- data.frame(
    severity      = c("critical", "warning", "info"),
    category      = c("convergence", "condition_number", "mu_referencing"),
    message       = c("did not converge", "ill-conditioned", "mu detected"),
    source_method = c("focei", "focei", "focei"),
    stringsAsFactors = FALSE
  )
  fit <- make_fake_fit(warnings_structured = df)
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})

# ferx-core#1021: the `sir` category now dispatches on message content. The
# proposal-conditioning diagnostics are a statement about the model (parameters
# the data do not identify), not about SIR tuning, so they must not be answered
# with "raise sir_samples".
test_that(".ferx_warning_guidance dispatches sir by message content", {
  g <- function(msg) ferx:::.ferx_warning_guidance("sir", message = msg)

  rank_def <- g(paste0(
    "SIR: proposal covariance is rank-deficient beyond the FIX-ed parameters: ",
    "1 direction(s) carry no uncertainty [TVCL +0.71, TVV -0.70]."
  ))
  expect_match(rank_def, "not identified", fixed = TRUE)
  expect_false(grepl("sir_samples", rank_def, fixed = TRUE))

  shrunk <- g(paste0(
    "SIR: proposal was shrunk in 3 direction(s) so draws stay inside the ",
    "parameter bounds [TVKA -0.74 (sd 6.15e3 -> 2.89e0)]."
  ))
  expect_match(shrunk, "understate", fixed = TRUE)
  expect_false(grepl("sir_samples", shrunk, fixed = TRUE))

  # Availability / tuning messages keep the original advice.
  tuning <- g("SIR requested but covariance matrix is not available")
  expect_match(tuning, "sir_samples", fixed = TRUE)
})
