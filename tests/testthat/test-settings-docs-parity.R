# Every `settings` key the engine accepts must be documented in ?ferx_fit
# (#509). The engine's list comes from `ferx_rust_fit_option_keys()`
# (ferx-core's `framework_keys()`, every method's `method_specific_keys()` and
# the data-selection keys), as `ferx_rust_known_blocks()` does for blocks, so
# a key the engine gains without R docs fails here instead of going quietly
# undocumented.

# Keys the engine accepts that ?ferx_fit deliberately has no `settings` entry
# for, each with the reason. Every other key is either an entry of the
# `settings` argument or a `ferx_fit()` argument of its own.
sd_exempt <- c(
  frem_predictions = "structural map written by ferx_model_to_frem()",
  frem_sigma = "structural map written by ferx_model_to_frem()",
  bloq = "alias of the bloq_method argument",
  gradient_method = "alias of the gradient argument",
  ignore_subjects = "what the ignore_ids argument sends"
)

sd_rd <- function() {
  src <- testthat::test_path("..", "..", "man", "ferx_fit.Rd")
  if (file.exists(src)) tools::parse_Rd(src) else tools::Rd_db("ferx")[["ferx_fit.Rd"]]
}

sd_tag <- function(x) attr(x, "Rd_tag") %||% ""

sd_text <- function(x) paste(unlist(x), collapse = "")

# The `\item{name}{...}` entries of the Rd's \arguments section, named.
sd_arguments <- function(rd) {
  args <- Filter(function(x) identical(sd_tag(x), "\\arguments"), rd)[[1L]]
  items <- Filter(function(x) identical(sd_tag(x), "\\item"), args)
  stats::setNames(lapply(items, `[[`, 2L), vapply(items, function(x) sd_text(x[[1L]]), ""))
}

# Every `\code{key}` in the label of an `\item` nested anywhere in `x`: the
# keys a `settings` entry documents, aliases included.
sd_item_keys <- function(x) {
  if (!is.list(x)) return(character(0))
  here <- character(0)
  if (identical(sd_tag(x), "\\item") && length(x) == 2L) {
    codes <- Filter(function(y) identical(sd_tag(y), "\\code"), x[[1L]])
    here <- vapply(codes, sd_text, "")
  }
  c(here, unlist(lapply(x, sd_item_keys), use.names = FALSE))
}

# `x` with every nested `\item` whose label names `key` removed.
sd_drop_item <- function(x, key) {
  if (!is.list(x)) return(x)
  keep <- vapply(x, function(y) {
    !(identical(sd_tag(y), "\\item") && length(y) == 2L && key %in% sd_item_keys(list(y)))
  }, logical(1))
  out <- lapply(x[keep], sd_drop_item, key = key)
  attributes(out) <- attributes(x)[setdiff(names(attributes(x)), "names")]
  names(out) <- names(x)[keep]
  out
}

# The engine keys `settings_rd` (the `settings` argument's Rd) leaves
# undocumented, apart from `ferx_fit()` arguments and the exemptions.
sd_undocumented <- function(settings_rd, arguments) {
  keys <- ferx:::ferx_rust_fit_option_keys()
  setdiff(keys, c(sd_item_keys(settings_rd), arguments, names(sd_exempt)))
}

test_that("the engine reports its settings keys", {
  keys <- ferx:::ferx_rust_fit_option_keys()
  expect_type(keys, "character")
  expect_identical(keys, sort(unique(keys)))
  # One key from each source the glue unions: framework, a method arm, and the
  # data-selection keys.
  expect_true(all(c("covariance_method", "vi_lr", "imp_samples", "ignore_subjects") %in% keys))
})

test_that("?ferx_fit documents every settings key the engine accepts (#509)", {
  args <- sd_arguments(sd_rd())
  expect_identical(sd_undocumented(args[["settings"]], names(args)), character(0))
})

test_that("a deleted settings entry fails the parity check (#509)", {
  args <- sd_arguments(sd_rd())
  expect_identical(
    sd_undocumented(sd_drop_item(args[["settings"]], "vi_lr"), names(args)),
    "vi_lr"
  )
  # An entry naming two keys documents both, so dropping it loses both.
  expect_identical(
    sd_undocumented(sd_drop_item(args[["settings"]], "iscale_max"), names(args)),
    c("iscale_max", "iscale_min")
  )
})

test_that("?ferx_fit's settings entries name only keys the engine accepts (#509)", {
  args <- sd_arguments(sd_rd())
  documented <- sd_item_keys(args[["settings"]])
  expect_identical(setdiff(documented, ferx:::ferx_rust_fit_option_keys()), character(0))
})

test_that("the exemptions are engine keys with no settings entry (#509)", {
  args <- sd_arguments(sd_rd())
  expect_true(all(names(sd_exempt) %in% ferx:::ferx_rust_fit_option_keys()))
  expect_identical(intersect(names(sd_exempt), sd_item_keys(args[["settings"]])), character(0))
})

test_that("a settings key with a dedicated argument is refused by name (#509)", {
  ex <- ferx_example("warfarin")
  reserved <- c("method", "covariance", "verbose", "bloq_method", "bloq",
                "threads", "mu_referencing", "sir", "gradient", "gradient_method")
  for (key in reserved) {
    # A key the model file also sets warns about the override before the
    # engine refuses it.
    expect_error(
      suppressWarnings(
        ferx_fit(ex$model, ex$data, settings = stats::setNames(list(TRUE), key))
      ),
      paste0("setting `", key, "` conflicts with a dedicated ferx_fit() argument"),
      fixed = TRUE,
      info = key
    )
  }
})

test_that("optimizer_trace / scale_params: an explicit argument wins over settings (#509)", {
  ex <- ferx_example("warfarin")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data,
                  optimizer_trace = FALSE, scale_params = FALSE,
                  settings = list(optimizer_trace = TRUE, scale_params = TRUE),
                  covariance = FALSE, verbose = FALSE))
  expect_false(fit$call_settings$optimizer_trace)
  expect_false(fit$call_settings$scale_params)
  expect_null(fit$trace_path)

  # Left at its default, the argument does not shadow a settings value.
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data,
                  settings = list(scale_params = TRUE),
                  covariance = FALSE, verbose = FALSE))
  expect_true(fit$call_settings$scale_params)
})

test_that("an explicit FALSE overrides the model file's [fit_options] (#509)", {
  # optimizer_trace, scale_params and inits_from_nca, which ferx_fit() merges
  # into the settings it sends.
  ex <- ferx_example("warfarin")
  lines <- readLines(ex$model)
  at <- grep("^\\[fit_options\\]", lines)
  expect_length(at, 1L)
  model <- tempfile(fileext = ".ferx")
  on.exit(unlink(model), add = TRUE)
  # `parameter_scaling = none`, since the default `auto` supersedes
  # scale_params and would hide whether it reached the engine.
  scaled <- c("  parameter_scaling = none", "  optimizer_trace = true", "  scale_params = true",
              "  inits_from_nca = nca_sweep")
  writeLines(append(lines, scaled, after = at), model)
  plain <- tempfile(fileext = ".ferx")
  on.exit(unlink(plain), add = TRUE)
  writeLines(append(lines, scaled[1L], after = at), plain)
  run <- function(...) {
    suppressWarnings(ferx_fit(model, ex$data, covariance = FALSE, verbose = FALSE, ...))
  }

  # Left at their defaults, the arguments leave the model file's values on.
  on_fit <- run()
  expect_false(is.null(on_fit$trace_path))

  # Explicit FALSE turns both off, also when settings repeats the file's value:
  # no trace, and the fit of the model without either key.
  baseline <- suppressWarnings(ferx_fit(plain, ex$data, covariance = FALSE, verbose = FALSE))
  for (settings in list(NULL, list(optimizer_trace = TRUE, scale_params = TRUE,
                                   inits_from_nca = TRUE))) {
    off_fit <- run(optimizer_trace = FALSE, scale_params = FALSE, inits_from_nca = FALSE,
                   settings = settings)
    expect_false(off_fit$call_settings$inits_from_nca)
    expect_null(off_fit$trace_path)
    expect_identical(off_fit$theta, baseline$theta)
    expect_identical(off_fit$ofv, baseline$ofv)
  }
  # The file's keys do change this fit, so the identity above shows the
  # arguments reached the engine.
  expect_false(identical(on_fit$theta, baseline$theta))
})
