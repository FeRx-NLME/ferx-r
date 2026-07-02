#' Retrieve the result of a background fit
#'
#' Blocks until the background fit started by \code{\link{ferx_fit_async}}
#' completes, optionally printing live optimizer trace progress to the console.
#' Pressing \code{Ctrl+C} cancels the wait and kills the background process;
#' \code{NULL} is returned.
#'
#' @param handle A \code{ferx_job} handle returned by
#'   \code{\link{ferx_fit_async}}.
#' @param verbose Logical. Whether to print optimizer trace progress while
#'   waiting. Default \code{TRUE}.
#'
#' @return The \code{ferx_fit} object produced by the background fit.
#'   Returns \code{NULL} if the user interrupts.
#'
#' @section Errors:
#' If the background fit itself errors (e.g. bad model file, engine failure),
#' \code{ferx_collect()} re-throws the error in the calling session.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#' handle <- ferx_fit_async(ex$model, ex$data, method = "focei")
#' fit <- ferx_collect(handle)
#' summary(fit)
#' }
#'
#' @family fitting
#' @export
ferx_collect <- function(handle, verbose = TRUE) {
  if (!inherits(handle, "ferx_job")) {
    stop("`handle` must be a ferx_job object returned by ferx_fit_async()")
  }
  if (handle$backend == "callr") {
    result <- .ferx_collect_callr(handle, verbose)
  } else {
    result <- .ferx_collect_rstudio(handle, verbose)
  }
  if (is.null(result)) return(invisible(NULL))
  if (!handle$user_wanted_trace && !is.null(result$trace_path) &&
        file.exists(result$trace_path)) {
    unlink(result$trace_path)
    result$trace_path <- NULL
  }
  result
}

.ferx_collect_callr <- function(handle, verbose) {
  bg            <- handle$bg
  sidecar_path  <- handle$sidecar_path
  tail_n        <- handle$tail_n
  poll_interval <- handle$poll_interval
  trace_path    <- NULL
  state         <- list(last_iter = -1L, lines_printed = 0L)

  tryCatch({
    while (bg$is_alive()) {
      # Drain the outer wrapper's streams to prevent pipe buffer blocking.
      # The inner process's stderr was already consumed by the wrapper, so
      # the trace path comes from the sidecar file, not these streams.
      tryCatch(bg$read_output_lines(), error = function(e) NULL)
      tryCatch(bg$read_error_lines(),  error = function(e) NULL)
      if (verbose && is.null(trace_path))
        trace_path <- .ferx_read_sidecar(sidecar_path)
      if (verbose && !is.null(trace_path) && file.exists(trace_path)) {
        state <- .ferx_print_trace_tail(trace_path, tail_n, state)
      }
      Sys.sleep(poll_interval)
    }
    result <- bg$get_result()
    if (!is.null(sidecar_path) && file.exists(sidecar_path))
      unlink(sidecar_path)
    result
  }, interrupt = function(e) { # nocov start
    if (bg$is_alive()) bg$kill()
    if (!is.null(sidecar_path) && file.exists(sidecar_path))
      unlink(sidecar_path)
    message("\nFit interrupted by user.")
    NULL
  }) # nocov end
}

.ferx_collect_rstudio <- function(handle, verbose) {
  job_id        <- handle$job_id
  rds_path      <- handle$rds_path
  sidecar_path  <- handle$sidecar_path
  tail_n        <- handle$tail_n
  poll_interval <- handle$poll_interval
  trace_path    <- NULL
  state         <- list(last_iter = -1L, lines_printed = 0L)

  on_done <- function() {
    if (!is.null(handle$script_path) && file.exists(handle$script_path))
      unlink(handle$script_path)
  }

  # Declare outside the repeat so it's always in scope after the loop,
  # including when the loop body runs zero times (job already done).
  job_state <- "unknown"

  tryCatch({
    repeat {
      job_state <- tryCatch(
        rstudioapi::jobGetState(job_id),
        error = function(e) "unknown"
      )
      if (job_state %in% c("succeeded", "failed")) break

      if (verbose) {
        if (is.null(trace_path))
          trace_path <- .ferx_read_sidecar(sidecar_path)
        if (!is.null(trace_path) && file.exists(trace_path)) {
          state <- .ferx_print_trace_tail(trace_path, tail_n, state)
        }
      }
      Sys.sleep(poll_interval)
    }
    on_done()
    if (job_state == "failed" || !file.exists(rds_path)) {
      stop("Background fit failed. Check the Jobs pane for details.")
    }
    result <- readRDS(rds_path)
    unlink(rds_path)
    if (file.exists(sidecar_path)) unlink(sidecar_path)
    result
  }, interrupt = function(e) {
    tryCatch(rstudioapi::jobRemove(job_id), error = function(e) NULL)
    on_done()
    if (file.exists(sidecar_path)) unlink(sidecar_path)
    message("\nFit interrupted by user.")
    NULL
  })
}

#' @rdname ferx_fit_async
#' @param x A `ferx_job` handle returned by `ferx_fit_async()`.
#' @export
print.ferx_job <- function(x, ...) {
  status <- if (x$backend == "callr") {
    if (tryCatch(x$bg$is_alive(), error = function(e) FALSE)) "running" else "done"
  } else {
    tryCatch(rstudioapi::jobGetState(x$job_id), error = function(e) "unknown")
  }
  backend_note <- if (x$backend == "rstudio") {
    "rstudio (visible in Jobs pane)"
  } else {
    "callr"
  }
  cat(sprintf("<ferx_job> [%s] %s\n", status, x$model_label))
  cat(sprintf("Backend : %s\n", backend_note))
  cat("Call ferx_collect(handle) to retrieve the result.\n")

  # Show current trace progress when available.
  # Both backends write a sidecar file with the trace CSV path once the
  # first iteration is recorded; read it non-destructively here.
  trace_path <- .ferx_read_sidecar(x$sidecar_path)
  if (!is.null(trace_path)) {
    tr <- tryCatch(ferx_trace(trace_path), error = function(e) NULL)
    if (!is.null(tr) && nrow(tr) > 0L) {
      cat("\n")
      tbl <- .runlog_iter_table(tr, truncate = TRUE,
                                trunc_total = 20L,
                                trunc_head  = 5L,
                                trunc_tail  = 5L)
      cat(paste(tbl, collapse = "\n"), "\n")
      cat(sprintf("  %d iteration(s) so far\n", nrow(tr)))
    }
  }

  invisible(x)
}

# Parse a trace CSV path from a character vector of stdout/stderr lines.
# ferx-core emits "[ferx] optimizer trace -> /path/to/file.csv" on stderr.
#
# The path may contain spaces (e.g. when TMPDIR resolves under
# "/Users/Some User/...") so we anchor on ".csv" at end-of-line rather
# than refusing any whitespace inside the match. Greedy capture from the
# first "/" through the trailing ".csv" on the announcement line, with a
# Read the trace-CSV path from a sidecar file.
# Returns the path string when the sidecar exists and contains a valid
# path to an existing file; returns NULL in every other case (missing
# sidecar, empty/blank content, or the referenced trace CSV not yet on
# disk).  Wraps readLines in tryCatch so a race between file.exists() and
# readLines() (sidecar deleted by the background process at collect time)
# cannot propagate as an unhandled error.
.ferx_read_sidecar <- function(sidecar_path) {
  if (is.null(sidecar_path) || !file.exists(sidecar_path)) return(NULL)
  sp <- tryCatch(trimws(readLines(sidecar_path, warn = FALSE)[1L]),
                 error = function(e) "")
  if (nzchar(sp) && file.exists(sp)) sp else NULL
}

# Keep the old name as an internal alias so existing tests that mock the
# bg-object interface still work; .ferx_find_trace_from_bg() now delegates
# to the lines-based helper.
.ferx_find_trace_from_bg <- function(bg) {
  stdout_lines <- tryCatch(bg$read_output_lines(), error = function(e) character(0))
  stderr_lines <- tryCatch(bg$read_error_lines(), error = function(e) character(0))
  .ferx_find_trace_from_lines(c(stdout_lines, stderr_lines))
}

# Tail the trace CSV in place: each poll reads the file, prints any new rows
# beyond `state$last_iter`, and updates the state. Append-only display so it
# works across terminals without relying on ANSI cursor movement.
.ferx_print_trace_tail <- function(path, tail_n, state) {
  tr <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(tr) || nrow(tr) == 0L) return(state)
  last <- tr$iter[nrow(tr)]
  if (is.null(last) || !is.finite(last) || last <= state$last_iter) return(state)

  if (state$last_iter < 0L) {
    cat(sprintf("  %5s  %12s  %12s  %12s\n",
                "iter", "OFV", "grad_norm", "optimizer"))
  }
  new_rows <- tr[tr$iter > state$last_iter, , drop = FALSE]
  if (nrow(new_rows) > tail_n) {
    new_rows <- utils::tail(new_rows, tail_n)
  }
  for (i in seq_len(nrow(new_rows))) {
    row <- new_rows[i, ]
    gn <- row$grad_norm
    gn_str <- if (is.null(gn) || is.na(gn)) "         -  " else sprintf("%12.4g", gn)
    opt_str <- row$optimizer %||% ""
    cat(sprintf("  %5d  %12.4f  %s  %12s\n",
                as.integer(row$iter), row$ofv, gn_str, opt_str))
  }
  state$last_iter <- last
  state$lines_printed <- state$lines_printed + nrow(new_rows)
  state
}
