## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(ferx)


## ----new-print----------------------------------------------------------------
ferx_model_new(print = TRUE)


## ----new-print-2cpt, eval = FALSE---------------------------------------------
# ferx_model_new(template = "2cpt_iv", print = TRUE)


## ----new-write, eval = FALSE--------------------------------------------------
# ferx_model_new("my_model.ferx", edit = FALSE)


## ----section------------------------------------------------------------------
ex <- ferx_example("warfarin")
ferx_model_section(ex$model, "parameters")


## ----section-capture----------------------------------------------------------
params <- ferx_model_section(ex$model, "parameters")
params


## ----read-modify-write--------------------------------------------------------
# work in a fresh temp file so re-knitting the vignette never reuses stale state
model_path <- tempfile(fileext = ".ferx")
stopifnot(file.copy(ex$model, model_path))

# read the section
lines <- ferx_model_section(model_path, "parameters")

# modify: change the TVCL initial value from 0.134 to 0.5
lines <- sub("TVCL\\([^,]+,", "TVCL(0.5,", lines)

# write back
ferx_model_set_section(model_path, "parameters", lines)

# verify
ferx_model_section(model_path, "parameters")


## ----console-only, eval = FALSE-----------------------------------------------
# model_path <- file.path(tempdir(), "run1.ferx")
# 
# # 1. Write a skeleton
# ferx_model_new(model_path, template = "1cpt_oral", edit = FALSE)
# 
# # 2. Switch method to focei
# ferx_model_set_section(model_path, "fit_options", c(
#   "  method     = focei",
#   "  maxiter    = 300",
#   "  covariance = true"
# ))
# 
# # 3. Fit (uses the bundled warfarin data for illustration)
# ex   <- ferx_example("warfarin")
# fit  <- ferx_fit(model_path, ex$data)
# print(fit)


## ----overwrite-error----------------------------------------------------------
model_path <- tempfile(fileext = ".ferx")
ferx_model_new(model_path, edit = FALSE)            # creates the file

tryCatch(
  ferx_model_new(model_path, edit = FALSE),          # would overwrite — error
  error = function(e) message(conditionMessage(e))
)


## ----overwrite-ok-------------------------------------------------------------
ferx_model_new(model_path, template = "1cpt_iv", overwrite = TRUE, edit = FALSE)
ferx_model_section(model_path, "structural_model")

