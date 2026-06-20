## Initialize
library(ferx)
library(vpc)
library(dplyr)
library(ggplot2)
# load_all()

ex <- ferx_example("warfarin_bloq")

## Fit: bloq_method="m3" enables Beal's M3 likelihood for observations
## flagged in the dataset. CENS=1 means DV is an LLOQ; CENS=-1 means DV
## is a ULOQ. This example uses lower-limit censoring. The model file also
## sets this, so the argument below is optional; we pass it explicitly for
## clarity.
result <- ferx_fit(
  model = ex$model,
  data = ex$data,
  method = "focei",
  bloq_method = "m3"
)

data <- read.csv(file = ex$data)
data |>
  filter(EVID == 0) |>
  mutate(DV = as.numeric(DV)) |>
  ggplot(aes(x = TIME, y = DV, group = ID)) +
    geom_line() +
    geom_point()

# Simulate at the fitted estimates (posterior-predictive / VPC flow).
sim <- ferx_simulate(
  ex$model,
  ex$data,
  n_sim = 200,
  seed = 42,
  fit = result
)
obs <- read.csv(file = ex$data)
vpc(
  obs = obs |>
    mutate(DV = as.numeric(DV)),
  sim = sim,
  sim_cols = list(dv = "DV_SIM")
)

## VPC for lower-limit censoring
vpc_cens(
  obs = obs |>
    mutate(DV = as.numeric(DV)),
  sim = sim,
  sim_cols = list(dv = "DV_SIM"),
  lloq = 2
)
