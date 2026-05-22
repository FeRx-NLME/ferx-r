library(deSolve)

set.seed(42)

# True population parameters
TVCL  <- 5.0;  TVV1 <- 50.0; TVQ <- 10.0; TVV2 <- 100.0
TVMTT <- 1.0;  TVLAG <- 0.25
KTR   <- 4.0 / TVMTT

omega_CL <- 0.09; omega_V1 <- 0.09; omega_KA <- 0.09
sigma_prop <- 0.30

DOSE <- 250   # mg, oral, CMT=1
obs_times <- c(0.5, 1, 2, 4, 6, 8, 12, 18, 24)  # hours
n_subj <- 20

wt_vec <- round(runif(n_subj, 50, 100), 1)

# ODE system (amounts for transit, concentration for central/peripheral)
transit_2cpt_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dT1 <- -KTR * T1
    dT2 <-  KTR * T1 - KTR * T2
    dT3 <-  KTR * T2 - KTR * T3
    dC1 <-  KA  * T3 / V1 - (CL / V1 + Q / V1) * C1 + Q / V2 * C2
    dC2 <-  Q   * C1 / V1 - Q / V2 * C2
    list(c(dT1, dT2, dT3, dC1, dC2))
  })
}

records <- list()

for (i in seq_len(n_subj)) {
  wt <- wt_vec[i]

  # Individual parameters with IIV
  CL_i <- TVCL * (wt / 70)^0.75 * exp(rnorm(1, 0, sqrt(omega_CL)))
  V1_i <- TVV1 * (wt / 70)^1.00 * exp(rnorm(1, 0, sqrt(omega_V1)))
  Q_i  <- TVQ
  V2_i <- TVV2 * (wt / 70)^1.00
  KA_i <- KTR  * exp(rnorm(1, 0, sqrt(omega_KA)))

  parms <- c(CL = CL_i, V1 = V1_i, Q = Q_i, V2 = V2_i,
             KTR = KTR, KA = KA_i)

  # Initial state: dose added to transit1 at t = TVLAG
  state0 <- c(T1 = DOSE, T2 = 0, T3 = 0, C1 = 0, C2 = 0)
  t_solve <- c(0, obs_times - TVLAG)   # solve from lag-adjusted time zero
  t_solve <- pmax(t_solve, 0)
  t_solve <- sort(unique(t_solve))

  sol <- ode(y = state0, times = t_solve, func = transit_2cpt_ode,
             parms = parms, method = "lsoda")
  sol_df <- as.data.frame(sol)

  # Interpolate IPRED at observation times (adjusted for lagtime)
  C1_at_obs <- approx(sol_df$time, sol_df$C1,
                      xout = pmax(obs_times - TVLAG, 0),
                      rule = 2)$y

  DV <- C1_at_obs * (1 + rnorm(length(obs_times), 0, sigma_prop))
  DV <- pmax(DV, 0)   # no negative concentrations

  # Dose record (EVID=1)
  dose_row <- data.frame(
    ID   = i, TIME = 0, DV = NA, EVID = 1,
    AMT  = DOSE, CMT = 1, MDV = 1, WT = wt
  )

  # Observation records (EVID=0)
  obs_rows <- data.frame(
    ID   = i, TIME = obs_times,
    DV   = round(DV, 4),
    EVID = 0, AMT = NA, CMT = 4, MDV = 0, WT = wt
  )

  records[[i]] <- rbind(dose_row, obs_rows)
}

dat <- do.call(rbind, records)
dat$DV[dat$EVID == 1] <- "."
dat$AMT[dat$EVID == 0] <- "."

cat("Dataset: n =", n_subj, "subjects,", nrow(dat), "rows\n")
print(head(dat, 12))

write.csv(dat, "transit_2cpt.csv", row.names = FALSE, quote = FALSE, na = ".")
cat("Written: transit_2cpt.csv\n")
