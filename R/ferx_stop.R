#' Stop a background fit
#'
#' Terminates a background job started by \code{\link{ferx_fit_async}}
#' without waiting for it to finish. Unlike pressing \code{Ctrl+C} during
#' \code{\link{ferx_collect}} (which cancels an active wait and kills the
#' background process), \code{ferx_stop()} can be called at any time to
#' terminate the job directly.
#'
#' @param handle A \code{ferx_job} handle returned by
#'   \code{\link{ferx_fit_async}}.
#'
#' @return Invisibly, \code{TRUE} if a running job was stopped, or
#'   \code{FALSE} if the job had already finished (nothing to stop).
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#' handle <- ferx_fit_async(ex$model, ex$data, method = "focei")
#' ferx_stop(handle)
#' }
#'
#' @family fitting
#' @export
ferx_stop <- function(handle) {
  if (!inherits(handle, "ferx_job")) {
    stop("`handle` must be a ferx_job object returned by ferx_fit_async()")
  }
  stopped <- if (handle$backend == "callr") {
    .ferx_stop_callr(handle)
  } else {
    .ferx_stop_rstudio(handle)
  }
  if (stopped) {
    message(sprintf("Stopped background fit for %s.", handle$model_label))
  } else {
    message(sprintf("%s already finished; nothing to stop.", handle$model_label))
  }
  invisible(stopped)
}

# The outer callr process supervises the inner ferx_fit process
# (`supervise = TRUE`, see .ferx_async_callr()), so killing the outer
# process is enough for callr's supervisor daemon to also terminate the
# inner one - the same mechanism .ferx_collect_callr() already relies on
# when a Ctrl+C interrupts an active wait.
.ferx_stop_callr <- function(handle) {
  bg <- handle$bg
  alive <- tryCatch(bg$is_alive(), error = function(e) FALSE)
  if (alive) bg$kill()
  if (!is.null(handle$sidecar_path) && file.exists(handle$sidecar_path)) {
    unlink(handle$sidecar_path)
  }
  alive
}

.ferx_stop_rstudio <- function(handle) {
  state <- tryCatch(rstudioapi::jobGetState(handle$job_id),
                     error = function(e) "unknown")
  running <- !identical(state, "succeeded") && !identical(state, "failed")
  if (running) {
    tryCatch(rstudioapi::jobRemove(handle$job_id), error = function(e) NULL)
  }
  if (!is.null(handle$script_path) && file.exists(handle$script_path)) {
    unlink(handle$script_path)
  }
  if (!is.null(handle$sidecar_path) && file.exists(handle$sidecar_path)) {
    unlink(handle$sidecar_path)
  }
  running
}
