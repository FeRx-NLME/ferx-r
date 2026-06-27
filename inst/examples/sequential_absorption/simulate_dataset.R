# Simulate the sequential (zero-order -> first-order) absorption example dataset.
#
# Independent of the ferx engine: integrates the depot/central ODEs with deSolve,
# so the bundled dataset is an external reference (mirrors transit_2cpt's
# simulate_dataset.R). The zero-order input fills the depot at a constant rate
# DOSE/DUR over (0, DUR], then the depot empties to central by first-order KA.
# The hard cutoff at t = DUR is handled by integrating in two exact segments
# ([0, DUR] with the input on, [DUR, end] with it off) rather than relying on the
# adaptive solver to land on the discontinuity.
#
# Run from this directory:  Rscript simulate_dataset.R
# Writes sequential_absorption.csv (copy to ../data/).

library(deSolve)

set.seed(7)

# True population parameters (match the sequential_absorption.ferx header).
TVCL <- 5.0; TVV <- 50.0; TVKA <- 1.0; TVDUR <- 3.0
omega_CL <- 0.09; omega_V <- 0.09
sigma_prop <- 0.15

DOSE <- 100            # mg, into the depot (CMT 1)
obs_times <- c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
n_subj <- 12

# depot/central ODE; `Rin` is the (constant) zero-order input into the depot.
seq_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    ddepot   <- Rin - KA * depot
    dcentral <- KA * depot - (CL / V) * central
    list(c(ddepot, dcentral))
  })
}

records <- list()

for (i in seq_len(n_subj)) {
  CL_i <- TVCL * exp(rnorm(1, 0, sqrt(omega_CL)))
  V_i  <- TVV  * exp(rnorm(1, 0, sqrt(omega_V)))
  KA_i <- TVKA
  DUR_i <- TVDUR

  base <- c(CL = CL_i, V = V_i, KA = KA_i)
  tmax <- max(obs_times)

  # Segment 1: input ON over [0, DUR].
  g1 <- sort(unique(c(seq(0, DUR_i, length.out = 200), obs_times[obs_times <= DUR_i])))
  s1 <- ode(y = c(depot = 0, central = 0), times = g1, func = seq_ode,
            parms = c(base, Rin = DOSE / DUR_i, DUR = DUR_i), method = "lsoda")
  end1 <- s1[nrow(s1), c("depot", "central")]

  # Segment 2: input OFF over [DUR, tmax], carrying the state.
  g2 <- sort(unique(c(seq(DUR_i, tmax, length.out = 400), obs_times[obs_times >= DUR_i])))
  s2 <- ode(y = c(depot = end1[["depot"]], central = end1[["central"]]),
            times = g2, func = seq_ode,
            parms = c(base, Rin = 0, DUR = DUR_i), method = "lsoda")

  sol <- rbind(as.data.frame(s1), as.data.frame(s2))
  sol <- sol[!duplicated(sol$time), ]

  central_obs <- approx(sol$time, sol$central, xout = obs_times, rule = 2)$y
  conc <- central_obs / V_i
  DV <- pmax(conc * (1 + rnorm(length(obs_times), 0, sigma_prop)), 0)

  dose_row <- data.frame(ID = i, TIME = 0, DV = NA, EVID = 1,
                         AMT = DOSE, CMT = 1, MDV = 1)
  obs_rows <- data.frame(ID = i, TIME = obs_times, DV = round(DV, 4),
                         EVID = 0, AMT = NA, CMT = 2, MDV = 0)
  records[[i]] <- rbind(dose_row, obs_rows)
}

dat <- do.call(rbind, records)
dat$DV[dat$EVID == 1] <- "."
dat$AMT[dat$EVID == 0] <- "."

cat("Dataset: n =", n_subj, "subjects,", nrow(dat), "rows\n")
print(head(dat, 14))

write.csv(dat, "sequential_absorption.csv", row.names = FALSE, quote = FALSE, na = ".")
cat("Written: sequential_absorption.csv\n")
