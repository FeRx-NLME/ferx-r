#' Load and validate a .ferxsearch search configuration
#'
#' Reads a \code{.ferxsearch} file with the engine's own loader and returns
#' everything it validated, so a search space is checked before the first
#' candidate is fitted rather than while watching fits fail. The loader is
#' strict: an unknown section, an unparseable \code{[space] mfl}, an empty
#' space, a feature the engine cannot express (a coverage gap), an unrecognised
#' \code{[rank] type} or an invalid \code{[rank.penalties]} charge are all
#' errors here, each naming the offender.
#'
#' The file format (ferx-core's search tooling) is TOML:
#'
#' \preformatted{
#' base = "../models/two_cpt_oral_cov.ferx"
#' data = "../data/two_cpt_oral_cov.csv"
#'
#' [space]
#' mfl = "COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])"
#'
#' [rank]
#' type = "bic"
#'
#' [strictness]
#' require_converged = true
#'
#' [run]
#' threads = 4
#' retries = 2
#' }
#'
#' \code{base}, \code{data} and \code{[run] cache_dir} are resolved relative to
#' the configuration file's own directory, so the file travels with its model.
#' A tool section (\code{[covsearch]}, \code{[modelsearch]},
#' \code{[iivsearch]}, \code{[iovsearch]}, \code{[ruvsearch]},
#' \code{[structsearch]}, \code{[allometry]}) is kept for the tool that owns it
#' and reported in \code{$tools}; any other section name is an error, so a
#' misspelt \code{[strictnes]} cannot silently leave the gate at its defaults.
#'
#' Two of the sections the engine accepts have no R tool to run them:
#' \code{[globalsearch]} and \code{[structsearch]}. A file carrying one loads
#' cleanly and then has that section ignored by whichever R tool it is handed
#' to - a \code{[globalsearch]} space run through \code{\link{ferx_covsearch}}
#' is a stepwise search, not the global one asked for - so loading such a file
#' warns, here and in every tool that takes \code{config =}.
#'
#' \code{[rank.penalties]} (the schedule a \code{[rank] type = "penalized"}
#' charges) comes back in \code{$rank$penalties} as the \emph{effective}
#' schedule: the file's keys overlaid on the engine's pyDarwin defaults, with
#' \code{$rank$penalties_set} naming the charges the file changed. A file that
#' restates a default is not distinguishable from one that omits it, since the
#' two describe the same run.
#'
#' @param path Path to a \code{.ferxsearch} file.
#'
#' @return An S3 object of class \code{ferx_search_config}: a list with
#'   \code{path}, \code{dir}, \code{base}, \code{data} (\code{NULL} when the
#'   file defers to the model's \code{[data]} block), \code{mfl} (the space
#'   source, verbatim), \code{space} (a data frame with one row per feature:
#'   \code{feature}, \code{keyword}, \code{optional}), \code{rank}
#'   (\code{type}, \code{cutoff}, \code{penalties} - a named numeric of the
#'   effective \code{[rank.penalties]} schedule - and \code{penalties_set},
#'   the charges the file changed), \code{strictness} (the \emph{effective}
#'   gate, the file's keys overlaid on the engine's defaults),
#'   \code{strictness_set} (the keys the file stated explicitly), \code{run}
#'   (\code{threads}, \code{retries}, \code{cache_dir}, \code{resume}) and
#'   \code{tools} (the tool sections the file carries).
#'
#' @examples
#' cfg <- ferx_search_config(ferx_example("two_cpt_oral_cov")$search)
#' cfg$space
#' cfg$rank$type
#'
#' @seealso \code{\link{ferx_search_space}} to expand the space against a
#'   model, \code{\link{ferx_search_coverage}} for the coverage table as data,
#'   and \code{\link{ferx_search_results}} to read a run's candidate table.
#' @family search
#' @export
ferx_search_config <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("'path' must be a single file path to a .ferxsearch file")
  }
  if (!file.exists(path)) stop("Search configuration not found: ", path)

  raw <- ferx_rust_search_config_load(normalizePath(path))

  cfg <- structure(
    list(
      path = normalizePath(path),
      dir  = as.character(raw$dir),
      base = as.character(raw$base),
      data = .ferx_chr_or_null(raw$data),
      mfl  = as.character(raw$mfl),
      space = .ferx_space_frame(raw$feature, raw$keyword, raw$optional),
      rank = list(
        type   = .ferx_chr_or_na(raw$rank_type),
        cutoff = .ferx_na_if_nan(raw$rank_cutoff),
        # The effective schedule, not the file's keys: what a `penalized`
        # criterion would charge. `penalties_set` is the overlay, computed in
        # the glue because only the engine holds the defaults.
        penalties     = stats::setNames(as.numeric(raw$penalty_value),
                                        as.character(raw$penalty_name)),
        penalties_set = as.character(raw$penalty_set)
      ),
      strictness = list(
        require_converged    = isTRUE(raw$require_converged),
        require_covariance   = isTRUE(raw$require_covariance),
        max_condition_number = .ferx_na_if_nan(raw$max_condition_number),
        max_correlation      = .ferx_na_if_nan(raw$max_correlation),
        reject_on_boundary   = isTRUE(raw$reject_on_boundary),
        reject_init_stall    = isTRUE(raw$reject_init_stall)
      ),
      strictness_set = as.character(raw$strictness_set),
      run = list(
        # 0 is the engine's "let the runner choose" sentinel for threads.
        threads   = if (as.integer(raw$threads) > 0L) as.integer(raw$threads) else NULL,
        retries   = as.integer(raw$retries),
        cache_dir = .ferx_chr_or_null(raw$cache_dir),
        resume    = isTRUE(raw$resume)
      ),
      tools = as.character(raw$tools)
    ),
    class = "ferx_search_config"
  )

  .ferx_search_warn_unconsumed(cfg$tools, "ferx_search_config")
  cfg
}

#' @param x A \code{ferx_search_config} object.
#' @param ... Ignored.
#' @rdname ferx_search_config
#' @export
print.ferx_search_config <- function(x, ...) {
  cat("<ferx search configuration>\n")
  cat("  file:  ", basename(x$path), "\n", sep = "")
  cat("  base:  ", x$base, "\n", sep = "")
  cat("  data:  ", x$data %||% "(from the model's [data] block)", "\n", sep = "")
  cat("\n")

  cat("Search space (", nrow(x$space), " feature",
      if (nrow(x$space) == 1L) "" else "s", "):\n", sep = "")
  cat("  mfl = ", x$mfl, "\n", sep = "")
  for (i in seq_len(nrow(x$space))) {
    cat(sprintf("  %-46s %s\n", x$space$feature[i],
                if (x$space$optional[i]) "(exploratory)" else "(structural)"))
  }
  cat("\n")

  cutoff <- if (is.na(x$rank$cutoff)) "(tool default)" else format(x$rank$cutoff)
  type   <- if (is.na(x$rank$type)) "(tool default)" else x$rank$type
  cat("Rank: ", type, "   cutoff: ", cutoff, "\n", sep = "")
  # The schedule only steers a run that ranks on it, or a global search's
  # search-level charges - so a file that neither asks for `penalized` nor
  # changes a charge gets no block, and the common case stays quiet.
  if (identical(x$rank$type, "penalized") || length(x$rank$penalties_set)) {
    cat("Penalties (* = set by the file):\n")
    for (nm in names(x$rank$penalties)) {
      mark <- if (nm %in% x$rank$penalties_set) "*" else " "
      # 15 digits, so the crash charge reads 99999999 rather than a rounded
      # 1e+08 - a penalty a user set is a number they should recognise.
      cat(sprintf("  %s %-22s %s\n", mark, nm,
                  format(x$rank$penalties[[nm]], digits = 15, trim = TRUE)))
    }
  }

  cat("Strictness (* = set by the file):\n")
  for (nm in names(x$strictness)) {
    mark <- if (nm %in% x$strictness_set) "*" else " "
    val <- x$strictness[[nm]]
    shown <- if (is.numeric(val) && is.na(val)) "off" else format(val)
    cat(sprintf("  %s %-22s %s\n", mark, nm, shown))
  }

  cat("Run: threads = ", x$run$threads %||% "auto",
      ", retries = ", x$run$retries,
      ", resume = ", x$run$resume,
      if (is.null(x$run$cache_dir)) "" else paste0(", cache_dir = ", x$run$cache_dir),
      "\n", sep = "")
  if (length(x$tools)) {
    cat("Tool sections: ", paste(x$tools, collapse = ", "), "\n", sep = "")
    orphans <- intersect(x$tools, .FERX_SEARCH_SECTIONS_WITHOUT_TOOL)
    if (length(orphans)) {
      cat("  no R tool runs: ", paste(orphans, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}

# -- Internal helpers shared by the search surface --

# NaN is the glue's "no value" sentinel for an optional engine numeric; R's own
# missing marker is NA, and `is.nan(NA_real_)` is FALSE so this is idempotent.
.ferx_na_if_nan <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.nan(x)) NA_real_ else x
}

# Same sentinel, for a field that stays in a fixed-shape list: an unstated
# `[rank] type` is NA rather than a dropped element, so `cfg$rank$type` is
# always there to test.
.ferx_chr_or_na <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L || !nzchar(x)) NA_character_ else x
}

# "" is the glue's sentinel for an absent path / string.
.ferx_chr_or_null <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L || !nzchar(x)) NULL else x
}

.ferx_space_frame <- function(feature, keyword, optional) {
  data.frame(
    feature  = as.character(feature),
    keyword  = as.character(keyword),
    optional = as.logical(optional),
    stringsAsFactors = FALSE
  )
}
