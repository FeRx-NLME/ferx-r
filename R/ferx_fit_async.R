#' Fit a model asynchronously
#'
#' Submits \code{\link{ferx_fit}} to a background process and returns
#' immediately, leaving the R session free. Call \code{\link{ferx_collect}}
#' on the returned handle to block-wait for the result and see live trace
#' progress.
#'
#' In RStudio the fit is submitted via \code{rstudioapi::jobRunScript()} and
#' appears in the Jobs pane. In Positron, plain R, or any other interactive
#' session it uses \code{callr::r_bg()}. Both backends return the same
#' \code{ferx_job} handle and the same \code{ferx_collect()} call retrieves
#' the result.
#'
#' The function forces \code{optimizer_trace = TRUE} internally so it has a
#' parseable progress channel. When you did not explicitly request the trace,
#' the temp CSV is deleted after \code{ferx_collect()} returns; pass
#' \code{optimizer_trace = TRUE} explicitly to keep it.
#'
#' @param model,data,... Forwarded to \code{\link{ferx_fit}}. See its
#'   documentation for the full argument list.
#' @param tail_n Integer. Number of recent trace rows to keep visible per
#'   poll in \code{\link{ferx_collect}}. Default 6.
#' @param poll_interval Numeric. Seconds between polls in
#'   \code{\link{ferx_collect}}. Default 0.5.
#'
#' @return A \code{ferx_job} handle. Pass it to \code{\link{ferx_collect}}
#'   to retrieve the \code{ferx_fit} result. In non-interactive sessions
#'   (\code{Rscript}, knitr, batch mode) the function falls back to a
#'   synchronous \code{\link{ferx_fit}} call and returns the result directly.
#'
#' @section Non-interactive use:
#' When called in a non-interactive session, \code{ferx_fit_async()} silently
#' falls back to \code{\link{ferx_fit}} since background jobs only add value
#' in an interactive console.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#' handle <- ferx_fit_async(ex$model, ex$data, method = "focei")
#' # R prompt is free immediately
#' print(handle)
#' fit <- ferx_collect(handle)
#' summary(fit)
#' }
#'
#' @family fitting
#' @export
ferx_fit_async <- function(model, data = NULL, ...,
                           tail_n = 6L,
                           poll_interval = 0.5) {
  if (!interactive()) {
    return(ferx_fit(model, data, ...))
  }
  if (!is.numeric(tail_n) || length(tail_n) != 1L || tail_n < 1L) {
    stop("`tail_n` must be a positive integer scalar")
  }
  if (!is.numeric(poll_interval) || length(poll_interval) != 1L ||
        poll_interval <= 0) {
    stop("`poll_interval` must be a positive numeric scalar")
  }
  tail_n <- as.integer(tail_n)
  dots <- list(...)
  user_wanted_trace <- isTRUE(dots$optimizer_trace)
  model_label <- if (is.character(model)) basename(model) else "model"

  use_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    rstudioapi::hasFun("jobRunScript")

  if (use_rstudio) {
    handle <- .ferx_async_rstudio(model, data, dots, model_label,
                                   tail_n, poll_interval, user_wanted_trace)
  } else {
    if (!requireNamespace("callr", quietly = TRUE)) {
      stop(
        "ferx_fit_async() requires the `callr` package. ",
        "Install with: install.packages(\"callr\")"
      )
    }
    handle <- .ferx_async_callr(model, data, dots, model_label,
                                 tail_n, poll_interval, user_wanted_trace)
  }

  cat(sprintf(
    "Fitting %s in background [%s]. Call ferx_collect(handle) when ready.\n",
    model_label, handle$backend
  ))
  invisible(handle)
}

.ferx_async_callr <- function(model, data, dots, model_label,
                               tail_n, poll_interval, user_wanted_trace) {
  sidecar_path <- tempfile("ferx_trace_loc_", fileext = ".txt")

  # Outer process monitors the inner ferx_fit process, drains its stderr to
  # extract the trace path, and writes it to sidecar_path so print.ferx_job
  # can read it non-destructively at any time.
  bg <- callr::r_bg(
    func = function(model, data, dots, sidecar_path) { # nocov start
      .find_trace <- utils::getFromNamespace(".ferx_find_trace_from_lines", "ferx")
      dots$optimizer_trace <- TRUE
      inner <- callr::r_bg(
        func = function(model, data, dots) {
          do.call(ferx::ferx_fit, c(list(model = model, data = data), dots))
        },
        args = list(model = model, data = data, dots = dots),
        package = TRUE, supervise = TRUE
      )
      trace_written <- FALSE
      while (inner$is_alive()) {
        if (!trace_written) {
          # Drain and scan for the trace announcement.
          lines <- c(
            tryCatch(inner$read_output_lines(), error = function(e) character(0)),
            tryCatch(inner$read_error_lines(),  error = function(e) character(0))
          )
          path <- .find_trace(lines)
          if (!is.null(path)) {
            writeLines(path, sidecar_path)
            trace_written <- TRUE
          }
        } else {
          # Trace found; keep draining to prevent the inner process from
          # blocking on a full OS pipe buffer (typically 64 KB on Linux/macOS).
          tryCatch(inner$read_output_lines(), error = function(e) NULL)
          tryCatch(inner$read_error_lines(),  error = function(e) NULL)
        }
        Sys.sleep(0.2)
      }
      # Drain remaining output after the process exits (fast fits may finish
      # before the first poll fires, leaving buffered stderr unread).
      if (!trace_written) {
        lines <- c(
          tryCatch(inner$read_output_lines(), error = function(e) character(0)),
          tryCatch(inner$read_error_lines(),  error = function(e) character(0))
        )
        path <- .find_trace(lines)
        if (!is.null(path)) writeLines(path, sidecar_path)
      }
      # Re-signal inner failures as plain conditions so the caller does not
      # see a double-nested "callr subprocess failed: callr subprocess failed"
      # message.  The original error message is preserved.
      tryCatch(
        inner$get_result(),
        error = function(e) stop(conditionMessage(e), call. = FALSE)
      )
    }, # nocov end
    args = list(model = model, data = data, dots = dots,
                sidecar_path = sidecar_path),
    package = TRUE,
    supervise = TRUE
  )
  structure(
    list(
      backend           = "callr",
      model_label       = model_label,
      bg                = bg,
      sidecar_path      = sidecar_path,
      tail_n            = tail_n,
      poll_interval     = poll_interval,
      user_wanted_trace = user_wanted_trace
    ),
    class = "ferx_job"
  )
}

.ferx_async_rstudio <- function(model, data, dots, model_label,
                                 tail_n, poll_interval, user_wanted_trace) {
  rds_path     <- tempfile("ferx_result_", fileext = ".rds")
  sidecar_path <- tempfile("ferx_trace_loc_", fileext = ".txt")

  script_path <- tempfile("ferx_job_", fileext = ".R")
  args_path   <- tempfile("ferx_args_", fileext = ".rds")
  saveRDS(list(model = model, data = data, dots = dots), args_path)

  # The job script deletes args_path itself immediately after readRDS so the
  # file is gone before callr::r_bg() is called - no race between cleanup and
  # the inner bg process starting up.
  writeLines(
    c(
      "local({",
      sprintf("  saved <- readRDS(%s)", deparse(args_path)),
      sprintf("  unlink(%s)", deparse(args_path)),
      "  model        <- saved$model",
      "  data         <- saved$data",
      "  dots         <- saved$dots",
      sprintf("  rds_path     <- %s", deparse(rds_path)),
      sprintf("  sidecar_path <- %s", deparse(sidecar_path)),
      "  dots$optimizer_trace <- TRUE",
      "  bg <- callr::r_bg(",
      "    func = function(model, data, dots) {",
      "      do.call(ferx::ferx_fit, c(list(model = model, data = data), dots))",
      "    },",
      "    args = list(model = model, data = data, dots = dots),",
      "    package = TRUE, supervise = TRUE",
      "  )",
      "  trace_written <- FALSE",
      "  while (bg$is_alive()) {",
      "    if (!trace_written) {",
      "      lines <- c(",
      "        tryCatch(bg$read_output_lines(), error = function(e) character(0)),",
      "        tryCatch(bg$read_error_lines(),  error = function(e) character(0))",
      "      )",
      "      path <- .find_trace(lines)",
      "      if (!is.null(path)) {",
      "        writeLines(path, sidecar_path)",
      "        trace_written <- TRUE",
      "      }",
      "    }",
      "    Sys.sleep(0.2)",
      "  }",
      "  # Drain remaining output after the process exits (fast fits may finish",
      "  # before the first poll fires, leaving buffered stderr unread).",
      "  if (!trace_written) {",
      "    lines <- c(",
      "      tryCatch(bg$read_output_lines(), error = function(e) character(0)),",
      "      tryCatch(bg$read_error_lines(),  error = function(e) character(0))",
      "    )",
      "    path <- .find_trace(lines)",
      "    if (!is.null(path)) writeLines(path, sidecar_path)",
      "  }",
      "  result <- bg$get_result()",
      "  saveRDS(result, rds_path)",
      "})"
    ),
    script_path
  )

  # If jobRunScript() throws, clean up both temp files before propagating.
  job_id <- tryCatch(
    rstudioapi::jobRunScript(
      path       = script_path,
      name       = model_label,
      workingDir = getwd()
    ),
    error = function(e) {
      unlink(c(script_path, args_path))
      stop(e)
    }
  )

  structure(
    list(
      backend           = "rstudio",
      model_label       = model_label,
      job_id            = job_id,
      rds_path          = rds_path,
      sidecar_path      = sidecar_path,
      script_path       = script_path,
      tail_n            = tail_n,
      poll_interval     = poll_interval,
      user_wanted_trace = user_wanted_trace
    ),
    class = "ferx_job"
  )
}
