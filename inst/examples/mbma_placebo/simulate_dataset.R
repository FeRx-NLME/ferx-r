# Simulate the mbma_placebo example dataset: a synthetic, arm-level,
# longitudinal model-based meta-analysis (MBMA).
#
# Independent of the ferx engine: the truth is a closed-form equation drawn
# with base R, so the bundled CSV is an external reference that does not move
# when the engine changes. ferx is used only at the end, to shape the summary
# table with ferx_mbma_data().
#
# The data are synthetic. The published naproxen MBMA data that ferx-core
# fits is licensed CC BY-NC, which is why it is not bundled here.
#
# Run from this directory:  Rscript simulate_dataset.R
# Writes ../data/mbma_placebo.csv. To write elsewhere, pass a path:
#   Rscript simulate_dataset.R /tmp/mbma_placebo.csv
# or set `out_csv` before sourcing this file.
#
# Truth (natural scale, a symptom score, lower = better):
#
#   y(s, a, t) = TVE0 + eta_s + P(s, t)
#                - EMAX * DOSE / (ED50 + DOSE) * t / (t + ET50) + kappa_sa
#
#   eta_s    ~ N(0, OMEGA_SD^2)        between-study variability on the baseline
#   P(s, t)  study-specific placebo time course, 0 at t = 0
#   kappa_sa ~ N(0, GAMMA2 / N_sa)     between-arm variability (BTAV)
#   arm mean ~ N(y, SD_PT^2 / N_sa),   reported SD from a chi-square draw
#
# Design: 6 studies, each with a placebo arm and 2-3 active dose arms, observed
# at week 0 and three of weeks 2, 4, 8, 12. That is 24 (study, week) cells: 18
# free placebo levels plus TVE0, EMAX, ED50, ET50, omega and gamma^2 is 24
# estimated parameters.
#
# TVE0 below is the nominal 50. Under contrast = sum_to_zero_within the model's
# TVE0 + ETA_E0 is each study's MEAN placebo level over its own visits, so the
# truth a fit recovers is the realised 50 + mean over studies of
# mean_t P(s, t). This script prints it (about 47.5).

library(ferx)

TVE0     <- 50     # nominal placebo level at week 0
EMAX     <- 12     # maximum drug effect (score points)
ED50     <- 40     # dose at half-maximal effect (mg)
ET50     <- 3      # time to half-maximal drug effect (weeks)
OMEGA_SD <- 3      # between-study SD of the baseline
GAMMA2   <- 100    # between-arm variance, for a 1-patient arm
SD_PT    <- 10     # patient-level SD

N_STUDY    <- 6
TIMES_POOL <- c(2, 4, 8, 12)
DOSE_POOL  <- c(10, 25, 50, 100, 200)

set.seed(370)

rows <- list()
placebo_mean <- numeric(0)
for (s in seq_len(N_STUDY)) {
  weeks <- sort(unique(c(0, sample(TIMES_POOL, 3))))
  doses <- c(0, sort(sample(DOSE_POOL, sample(2:3, 1))))
  eta <- rnorm(1, 0, OMEGA_SD)
  # Placebo: a study-specific improving trend plus visit noise, 0 at week 0.
  slope <- runif(1, 2, 8)
  P <- -slope * weeks / (weeks + 2) + c(0, rnorm(length(weeks) - 1, 0, 1))
  placebo_mean <- c(placebo_mean, mean(P))
  for (a in seq_along(doses)) {
    N <- sample(40:250, 1)
    kappa <- rnorm(1, 0, sqrt(GAMMA2 / N))
    D <- doses[a]
    mu <- TVE0 + eta + P - EMAX * D / (ED50 + D) * weeks / (weeks + ET50) + kappa
    arm_mean <- rnorm(length(weeks), mu, SD_PT / sqrt(N))
    arm_sd <- SD_PT * sqrt(rchisq(length(weeks), N - 1) / (N - 1))
    rows[[length(rows) + 1]] <- data.frame(
      trial = sprintf("T%02d", s), arm = a, DOSE = D, week = weeks,
      mean = round(arm_mean, 3), sd = round(arm_sd, 3), n = N)
  }
}
arms <- do.call(rbind, rows)

cat(sprintf("Realised TVE0 truth (nominal %g + mean placebo deviation): %.4f\n",
            TVE0, TVE0 + mean(placebo_mean)))

d <- ferx_mbma_data(arms, study = "trial", arm = "arm", time = "week",
                    mean = "mean", sd = "sd", n = "n", covariates = "DOSE")

if (!exists("out_csv", inherits = FALSE)) {
  args <- commandArgs(trailingOnly = TRUE)
  out_csv <- if (length(args)) args[1] else file.path("..", "data", "mbma_placebo.csv")
}
utils::write.csv(d, out_csv, row.names = FALSE)
cat("Wrote", nrow(d), "rows to", out_csv, "\n")
