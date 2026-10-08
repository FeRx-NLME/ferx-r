# Writes fixtures/theta_levels_pre466.fitrx: a level-block fit saved by a ferx
# that predates #466, so its level layout is under `r_extras$theta_levels` only
# (and, predating #430 and #412, without `value` and without `data_bindings`).
# test-theta-levels.R (T17e) loads it with today's ferx_load_fit().
#
# Committed as written by ferx 0.4.0.9000 at 6c7d02a (FeRx-NLME/ferx-r#427).
# Run it with that build, never with a current one, from tests/testthat:
#   Rscript fixtures/make-theta-levels-pre466.R
# The model and data are test-theta-levels.R's tl_col_model() / tl_data.
library(ferx)

model <- tempfile(fileext = ".ferx")
writeLines("[parameters]
  theta TVCL(2.0, 0.001, 20.0)
  theta PLACEBO[STUDY, TIME](0.0, -5.0, 5.0)
  theta TVV(10.0, 0.1, 500.0)
  omega ETA_V ~ 0.04
  sigma PROP_ERR ~ 0.05

[individual_parameters]
  CL = TVCL + PLACEBO
  V  = TVV
  Z  = TVV * exp(ETA_V)

[structural_model]
  pk one_cpt_iv(cl=CL, v=V)

[derived]
  Z_OUT = Z

[error_model]
  DV ~ proportional(PROP_ERR)

[fit_options]
  maxiter = 2
  inner_maxiter = 3
  covariance = false
", model)

data <- tempfile(fileext = ".csv")
writeLines("ID,TIME,DV,EVID,AMT,CMT,RATE,MDV,STUDY,PLA_IDX
1,0,.,1,100,1,0,1,1,1
1,1,8.1,0,.,1,0,0,1,1
1,4,6.2,0,.,1,0,0,1,2
1,12,3.1,0,.,1,0,0,1,3
2,0,.,1,100,1,0,1,2,4
2,1,7.4,0,.,1,0,0,2,4
2,4,5.5,0,.,1,0,0,2,5
2,12,2.8,0,.,1,0,0,2,6", data)

fit <- ferx_fit(model, data, verbose = FALSE)
ferx_save_fit(fit, "fixtures/theta_levels_pre466.fitrx", include_data = TRUE)
