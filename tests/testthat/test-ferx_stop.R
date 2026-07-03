# Tests for ferx_stop() and its per-backend helpers.

.stop_callr   <- getFromNamespace(".ferx_stop_callr",   "ferx")
.stop_rstudio <- getFromNamespace(".ferx_stop_rstudio", "ferx")
.async_callr  <- getFromNamespace(".ferx_async_callr",  "ferx")

make_fake_bg <- function(alive = TRUE) {
  state <- new.env()
  state$alive <- alive
  state$killed <- FALSE
  list(
    is_alive = function() state$alive,
    kill = function() {
      state$killed <- TRUE
      state$alive <- FALSE
      invisible(TRUE)
    },
    .state = state
  )
}

make_callr_handle <- function(alive = TRUE, sidecar_path = tempfile(fileext = ".txt")) {
  writeLines("dummy", sidecar_path)
  structure(
    list(backend = "callr", model_label = "warfarin", bg = make_fake_bg(alive),
         sidecar_path = sidecar_path, tail_n = 6L, poll_interval = 0.01,
         user_wanted_trace = FALSE),
    class = "ferx_job"
  )
}

make_rstudio_handle <- function(job_id = "job_1",
                                 sidecar_path = tempfile(fileext = ".txt"),
                                 script_path = tempfile(fileext = ".R")) {
  writeLines("dummy", sidecar_path)
  writeLines("dummy", script_path)
  structure(
    list(backend = "rstudio", model_label = "warfarin", job_id = job_id,
         rds_path = tempfile(fileext = ".rds"), sidecar_path = sidecar_path,
         script_path = script_path, tail_n = 6L, poll_interval = 0.01,
         user_wanted_trace = FALSE),
    class = "ferx_job"
  )
}

test_that("ferx_stop rejects a non-ferx_job handle", {
  expect_error(ferx::ferx_stop(list()), "ferx_job")
})

test_that(".ferx_stop_callr kills a running job and cleans up the sidecar", {
  h <- make_callr_handle(alive = TRUE)
  expect_true(.stop_callr(h))
  expect_true(h$bg$.state$killed)
  expect_false(file.exists(h$sidecar_path))
})

test_that(".ferx_stop_callr is a no-op on an already-finished job", {
  h <- make_callr_handle(alive = FALSE)
  expect_false(.stop_callr(h))
  expect_false(h$bg$.state$killed)
})

test_that(".ferx_stop_rstudio removes a running job and cleans up temp files", {
  skip_if_not_installed("mockery")
  h <- make_rstudio_handle()
  fn <- .stop_rstudio
  mockery::stub(fn, "rstudioapi::jobGetState", function(...) "running")
  removed <- FALSE
  mockery::stub(fn, "rstudioapi::jobRemove", function(...) removed <<- TRUE)

  expect_true(fn(h))
  expect_true(removed)
  expect_false(file.exists(h$sidecar_path))
  expect_false(file.exists(h$script_path))
})

test_that(".ferx_stop_rstudio is a no-op on an already-finished job", {
  skip_if_not_installed("mockery")
  h <- make_rstudio_handle()
  fn <- .stop_rstudio
  mockery::stub(fn, "rstudioapi::jobGetState", function(...) "succeeded")
  removed <- FALSE
  mockery::stub(fn, "rstudioapi::jobRemove", function(...) removed <<- TRUE)

  expect_false(fn(h))
  expect_false(removed)
})

test_that("ferx_stop reports whether it stopped a running or finished job", {
  running_handle  <- make_callr_handle(alive = TRUE)
  finished_handle <- make_callr_handle(alive = FALSE)

  expect_message(ferx::ferx_stop(running_handle), "Stopped background fit")
  expect_message(ferx::ferx_stop(finished_handle), "already finished")
})

test_that("ferx_stop kills a real background fit", {
  skip_if_not_installed("callr")
  ex <- ferx::ferx_example("warfarin")
  h  <- .async_callr(ex$model, ex$data,
                     list(settings = list(maxiter = 500L)),
                     "warfarin", 6L, 0.5, FALSE)
  on.exit({
    if (h$bg$is_alive()) h$bg$kill()
    if (file.exists(h$sidecar_path)) unlink(h$sidecar_path)
  }, add = TRUE)

  expect_true(h$bg$is_alive())
  expect_message(ferx::ferx_stop(h), "Stopped background fit")
  expect_false(h$bg$is_alive())
})
