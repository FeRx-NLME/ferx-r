## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ferx)


## ----example------------------------------------------------------------------
ex <- ferx_example("warfarin")
ex$model  # path to warfarin.ferx
ex$data   # path to warfarin.csv


## ----inspect-pre--------------------------------------------------------------
ferx_model_inspect(ex$model)


## ----inspect-programmatic-----------------------------------------------------
s <- ferx_model_inspect(ex$model)
s$theta_names   # population parameter names
s$model_type    # short label for the PK model
s$iiv           # IIV parameter names
s$residual      # error model type


## ----fit, eval = FALSE--------------------------------------------------------
# fit <- ferx_fit(ex$model, ex$data)
# fit


## ----inspect-post, eval = FALSE-----------------------------------------------
# ferx_model_inspect(fit)


## ----ferx-model---------------------------------------------------------------
ex <- ferx_example("warfarin")
m <- ferx_model(ex$data, ex$model)
m   # prints path, data path, and structure summary


## ----minimal-pipe, eval = FALSE-----------------------------------------------
# ferx_model(ex$data, ex$model) |>
#   ferx_fit(method = "focei", covariance = TRUE) |>
#   summary()


## ----get-section--------------------------------------------------------------
ferx_model(ex$data, ex$model) |>
  ferx_get_section("parameters")


## ----get-section-fit, eval = FALSE--------------------------------------------
# # Print [parameters] then fit without interrupting the chain
# fit <- ferx_model(ex$data, ex$model) |>
#   ferx_get_section("parameters") |>
#   ferx_fit(method = "focei")


## ----set-section, eval = FALSE------------------------------------------------
# model_copy <- file.path(tempdir(), "warfarin.ferx")
# file.copy(ex$model, model_copy)
# 
# fit <- ferx_model(ex$data, model_copy) |>
#   ferx_set_section("fit_options", c(
#     "  method     = focei",
#     "  maxiter    = 500",
#     "  covariance = true"
#   )) |>
#   ferx_fit()
# 
# summary(fit)
# ferx_model_inspect(fit)    # no path needed after fitting
# ferx_cor_matrix(fit)       # parameter correlation matrix
# ferx_plot_trace(fit)       # OFV + gradient norm convergence plot


## ----check-init, eval = FALSE-------------------------------------------------
# chk <- ferx_check_init(ex$model, ex$data, method = "focei")
# chk$summary      # ofv_start, ofv_end, ofv_drop — is OFV decreasing?
# ferx_plot_trace(chk$fit)   # visual check: first few iterations


## ----method-chain, eval = FALSE-----------------------------------------------
# fit <- ferx_model(ex$data, ex$model) |>
#   ferx_fit(method = c("saem", "focei"), covariance = TRUE)


## ----sim-pred, eval = FALSE---------------------------------------------------
# # Simulate 100 replicates using individual parameters from the fit
# sim <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42, fit = fit)
# 
# # Population predictions using fitted theta (eta = 0)
# pred <- ferx_predict(ex$model, ex$data, fit = fit)


## ----diagnostics, eval = FALSE------------------------------------------------
# head(fit$sdtab)               # PRED, IPRED, CWRES, IWRES per observation
# head(fit$ebe_etas)            # empirical Bayes ETAs per subject
# head(fit$individual_estimates) # individual PK parameters (CL, V, KA, ...)
# fit$eta_normality             # Shapiro-Wilk normality test per ETA

