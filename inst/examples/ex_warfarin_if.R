library(ferx)

# Conditional `if (cond) { ... } else { ... }` in [individual_parameters]:
# heavier subjects (WT > 70) get an allometric CL bump, lighter subjects use
# the typical clearance directly. Equivalent to NONMEM's IF/THEN/ENDIF or
# nlmixr2's if/else.
ex <- ferx_example("warfarin_if")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Per-subject CL values reflect the WT-driven branching:
head(fit$individual_estimates)
