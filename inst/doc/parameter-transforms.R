## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ferx)


## ----logn-inspect-------------------------------------------------------------
ex <- ferx_example("warfarin")
ferx_model_inspect(ex$model)


## ----logn-fit, eval = FALSE---------------------------------------------------
# fit <- ferx_fit(ex$model, ex$data)
# fit   # omega block shows "[log-normal]  CV% = ..."


## ----logn-estimates, eval = FALSE---------------------------------------------
# ferx_estimates(fit)


## ----logit-inspect------------------------------------------------------------
ex_logit <- ferx_example("warfarin_logit_f")
ferx_model_inspect(ex_logit$model)


## ----logit-fit, eval = FALSE--------------------------------------------------
# fit_logit <- ferx_fit(ex_logit$model, ex_logit$data)
# fit_logit


## ----logit-estimates, eval = FALSE--------------------------------------------
# est <- ferx_estimates(fit_logit)
# est[est$param == "THETA_F", ]
# # estimate         — logit scale
# # estimate_natural — (0, 1) scale (inv_logit back-transform)
# # lower_95_natural / upper_95_natural — asymmetric CI on (0, 1) scale


## ----add-inspect--------------------------------------------------------------
ex_add <- ferx_example("warfarin_additive_eta")
ferx_model_inspect(ex_add$model)


## ----add-fit, eval = FALSE----------------------------------------------------
# fit_add <- ferx_fit(ex_add$model, ex_add$data)
# fit_add


## ----add-estimates, eval = FALSE----------------------------------------------
# ferx_estimates(fit_add)

