library(ferx)

model_file <- file.path(
  system.file("examples/transit_2cpt", package = "ferx"),
  "transit_2cpt.ferx"
)
data_file <- file.path(
  system.file("examples/transit_2cpt", package = "ferx"),
  "transit_2cpt.csv"
)

# If running from the source tree directly:
if (!file.exists(model_file)) {
  model_file <- "transit_2cpt.ferx"
  data_file  <- "transit_2cpt.csv"
}

fit <- ferx_fit(model_file, data_file)
print(fit)
