# given: per-plot annual profits -> summarize mu(SR), cv(SR)
# or use your smooths for CV and mu vs sowndiv

library(dplyr)

# 1) fit your smooths (quadratic like your CV model)
fit_cv <- lm(cv ~ sowndiv + I(sowndiv^2), data = by_plot)   # by_plot: one row per plot with mean/sd/cv
fit_mu <- lm(mu ~ sowndiv + I(sowndiv^2), data = by_plot)

cv_fun  <- function(SR){ c <- coef(fit_cv); c[1] + c[2]*SR + c[3]*SR^2 }
dcv_fun <- function(SR){ c <- coef(fit_cv); c[2] + 2*c[3]*SR }

mu_fun  <- function(SR){ c <- coef(fit_mu); c[1] + c[2]*SR + c[3]*SR^2 }
dmu_fun <- function(SR){ c <- coef(fit_mu); c[2] + 2*c[3]*SR }

# 2) your dimensionless index
iv_index_fun <- function(SR){ - SR * dcv_fun(SR) }     # exactly your code’s definition

# 3) choose CARA risk aversion (sensitivity: try A in [0.001, 0.01])
A <- 0..001

# 4) converters to €/ha
ins_marg_eur <- function(SR){
  CV <- cv_fun(SR); MU <- mu_fun(SR)
  A * (CV*MU) * ( MU * iv_index_fun(SR)/SR - CV * dmu_fun(SR) )
}

rp_fun <- function(SR){
  CV <- cv_fun(SR); MU <- mu_fun(SR)
  0.5 * A * (CV*MU)^2
}

# 5) level insurance vs baseline SR0 (e.g. monocultures = 1)
SR0 <- 1
ins_level_eur <- function(SR){ rp_fun(SR0) - rp_fun(SR) }

# ----- Example evaluation on your SR values -----
SR_vals <- c(1,2,4,8,16)
out <- data.frame(
  SR = SR_vals,
  mu_hat = sapply(SR_vals, mu_fun),
  cv_hat = sapply(SR_vals, cv_fun),
  IV_index = sapply(SR_vals, iv_index_fun),
  Ins_marg_eur_per_species = sapply(SR_vals, ins_marg_eur),
  RP_eur = sapply(SR_vals, rp_fun),
  Ins_level_eur_vs_SR1 = sapply(SR_vals, ins_level_eur)
)
out
