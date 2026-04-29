# Build a minimal trace CSV mirroring the columns the Rust backend writes,
# so trace tests don't need a real fit.
write_fake_trace <- function(path = tempfile(fileext = ".csv"),
                             n = 5L,
                             method = "focei") {
  tr <- data.frame(
    iter             = seq_len(n),
    method           = method,
    phase            = "",
    ofv              = seq(100, 100 - n + 1),
    wall_ms          = seq(10, by = 10, length.out = n),
    grad_norm        = seq(1.0, by = -0.1, length.out = n),
    step_norm        = seq(0.5, by = -0.05, length.out = n),
    inner_iter_count = NA_integer_,
    optimizer        = "slsqp",
    lm_lambda        = NA_real_,
    ofv_delta        = NA_real_,
    step_accepted    = NA_integer_,
    cond_nll         = NA_real_,
    gamma            = NA_real_,
    mh_accept_rate   = NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(tr, path, row.names = FALSE)
  path
}

# Reset the package-internal last-trace state. Tests touch internals via :::
# because that's the contract the no-argument route depends on.
ferx_state <- function() getFromNamespace(".ferx_state", "ferx")

reset_ferx_state <- function() {
  st <- ferx_state()
  st$last_trace_path  <- NULL
  st$last_trace_time  <- NULL
  st$last_trace_model <- NULL
}
