# Build a minimal ferx_fit S3 object for unit tests that don't need a real fit.
make_fake_fit <- function(...) {
  defaults <- list(
    model_name    = "test_model",
    data_name     = "test_data",
    converged     = TRUE,
    method        = "focei",
    method_chain  = NULL,
    ofv           = -100.0,
    aic           = -90.0,
    bic           = -80.0,
    n_subjects    = 10L,
    n_obs         = 100L,
    n_parameters  = 5L,
    n_iterations  = 50L,
    theta         = c(CL = 1.0, V = 10.0),
    se_theta      = NULL,
    omega         = NULL,
    se_omega      = NULL,
    sigma         = NULL,
    se_sigma      = NULL,
    shrinkage_eta = NULL,
    shrinkage_eps = NULL,
    warnings      = character(0)
  )
  structure(modifyList(defaults, list(...)), class = "ferx_fit")
}

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
