# R/11b_prior_sensitivity.R
# Step 11b - PRIOR-SENSITIVITY CHECK on the Bayesian AHP weight posterior.
#            Companion to 11_bayesian_ahp.R. Answers the obvious examiner
#            question: "is your headline tau finding an artefact of the prior?"
#
# Why this EXISTS
#   The headline result of 11 is that the data-implied judgement noise tau is
#   larger than the Monte Carlo sigma. tau has a prior, tau ~ Half-Normal(0, s).
#   If the posterior tau (and the weight intervals) barely move when s changes
#   over a sensible range, the finding is driven by the DATA, not the prior --
#   which is exactly what you want to demonstrate. This script refits the model
#   under several priors and tabulates how little (or much) things move.
#
# What is VARIED
#   Primary : tau prior SD  in {0.5, 1, 2}   (1 = the baseline used in 11)
#   Optional: theta prior SD in {3} by default; add values to also probe it.
#   Everything else (model, data, seed, sampler depth) is held fixed, so any
#   change in the posterior is attributable to the prior alone.
#
# Engine
#   Stan only (uses R/stan/bayes_ahp_priors.stan, a parameterised-prior twin of
#   bayes_ahp.stan). If rstan is unavailable the script stops with a clear
#   message -- a prior sweep is a methodological robustness check best done with
#   the same engine as the main run, not the MH fallback.
#
# Reads : output/ahp_5factor_baseline.rds      (same matrix as 11)
#         R/stan/bayes_ahp_priors.stan
# Writes: output/tables/prior_sensitivity_tau.csv      (tau by prior)
#         output/tables/prior_sensitivity_weights.csv  (each weight by prior)
#         output/figures/prior_sensitivity_tau.png     (tau posteriors overlaid)
#
# READ THE RESULT LIKE THIS
#   If posterior tau median shifts by only a few hundredths and the weight CrIs
#   are visually indistinguishable across priors -> report "the tau finding is
#   robust to reasonable prior choices." If tau tracks the prior strongly, the
#   data are weakly informative about tau and you must say so (and lean on the
#   weakest-informative prior). Let the numbers decide; do not assume the answer.

cat("\n--- 11b_prior_sensitivity: does the tau finding survive the prior? ---\n")
set.seed(2026)

## ---- settings -------------------------------------------------------------
STAN_FILE<- "R/stan/bayes_ahp_priors.stan"
TAU_PRIORS<- c(0.5, 1, 2)     # half-normal SD on tau; 1 = baseline (matches 11)
THETA_PRIORS<- c(3)             # add e.g. c(1.5, 3, 6) to also probe theta prior
STAN_CHAINS<- 4
STAN_WARMUP<- 1000
STAN_ITER<- 2000             # 8000 draws/run, matching 11

## ---- guard: Stan required -------------------------------------------------
if (!requireNamespace("rstan", quietly = TRUE))
  stop("rstan is required for the prior-sensitivity sweep. Install rstan, or ",
       "run 11 with the MH engine and adapt this script's prior block by hand.")

## ---- load the same matrix as 11 ------------------------------------------
rds_path<- "output/ahp_5factor_baseline.rds"
if (!file.exists(rds_path))
  stop("Missing ", rds_path, " -- run R/01b_build_ahp_matrix.R first.")
if (!file.exists(STAN_FILE))
  stop("Missing ", STAN_FILE, " -- place bayes_ahp_priors.stan there.")

ahp<- readRDS(rds_path)
A<- ahp$matrix
stopifnot(is.matrix(A), nrow(A) == 5, ncol(A) == 5)
factors<- rownames(A); if (is.null(factors))
  factors<- c("Slope","LandCover","NDVI","SoilType","DistRoad")
n<- nrow(A)
idx<- which(upper.tri(A), arr.ind = TRUE)
y<- log(A[upper.tri(A)]); ii<- idx[,1]; jj<- idx[,2]; J<- length(y)

w_base<- if (!is.null(ahp$weights)) as.numeric(ahp$weights) else {
  ev<- eigen(A)$vectors[,1]; Re(ev)/sum(Re(ev)) }
names(w_base)<- factors

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
mod<- rstan::stan_model(STAN_FILE)   # compile once, reuse across priors

## ---- sweep ---------------------------------------------------------------
grid<- expand.grid(tau_sd = TAU_PRIORS, theta_sd = THETA_PRIORS)
tau_rows<- list(); w_rows<- list()
tau_draws_by_prior<- list()

for (r in seq_len(nrow(grid))) {
  ts<- grid$tau_sd[r]; th<- grid$theta_sd[r]
  cat(sprintf("\n[%d/%d] tau_prior_sd = %.2f, theta_prior_sd = %.2f\n",
              r, nrow(grid), ts, th))
  stan_data<- list(n=n, J=J, ii=ii, jj=jj, y=y,
                    tau_prior_sd = ts, theta_prior_sd = th)
  fit<- rstan::sampling(mod, data = stan_data,
    chains = STAN_CHAINS, warmup = STAN_WARMUP, iter = STAN_WARMUP + STAN_ITER,
    seed = 2026, refresh = 0)

  # convergence guard for THIS prior
  sm<- rstan::summary(fit, pars = c("tau","w"))$summary
  max_rhat<- max(sm[, "Rhat"], na.rm = TRUE)
  cat(sprintf("    max Rhat = %.4f %s\n", max_rhat,
              ifelse(max_rhat < 1.01, "(ok)", "(!! check convergence)")))

  tau<- as.numeric(rstan::extract(fit, "tau")$tau)
  W<- rstan::extract(fit, "w")$w; colnames(W)<- factors
  tau_draws_by_prior[[sprintf("tau~HN(0,%.1f)", ts)]]<- tau

  tau_rows[[r]]<- data.frame(
    tau_prior_sd = ts, theta_prior_sd = th,
    tau_med = round(median(tau), 3),
    tau_lo  = round(quantile(tau, 0.025), 3),
    tau_hi  = round(quantile(tau, 0.975), 3),
    max_rhat = round(max_rhat, 4)
  )
  w_rows[[r]]<- data.frame(
    tau_prior_sd = ts, theta_prior_sd = th,
    factor   = factors,
    baseline = round(as.numeric(w_base), 3),
    w_med    = round(apply(W, 2, median), 3),
    w_lo     = round(apply(W, 2, function(x) quantile(x, 0.025)), 3),
    w_hi     = round(apply(W, 2, function(x) quantile(x, 0.975)), 3)
  )
}

tau_tab<- do.call(rbind, tau_rows)
w_tab<- do.call(rbind, w_rows)

dir.create("output/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
write.csv(tau_tab, "output/tables/prior_sensitivity_tau.csv",     row.names = FALSE)
write.csv(w_tab,   "output/tables/prior_sensitivity_weights.csv", row.names = FALSE)

## ---- robustness verdict (computed, not assumed) --------------------------
tau_spread<- max(tau_tab$tau_med) - min(tau_tab$tau_med)
cat("\n=== tau posterior across priors ===\n"); print(tau_tab)
cat(sprintf("\nRange of posterior tau medians across priors: %.3f\n", tau_spread))
cat(sprintf("All posterior tau medians %s the Monte Carlo sigma = 0.3.\n",
            ifelse(all(tau_tab$tau_med > 0.3), "exceed", "do NOT all exceed")))
cat(ifelse(tau_spread < 0.10,
  ">> tau moves little across priors -> the finding is DATA-driven (robust).\n",
  ">> tau tracks the prior noticeably -> report this; lean on the weakest prior.\n"))

## ---- figure: overlaid tau posteriors -------------------------------------
png("output/figures/prior_sensitivity_tau.png", width = 1300, height = 900, res = 150)
dens<- lapply(tau_draws_by_prior, density)
xr<- range(sapply(dens, function(d) range(d$x)))
yr<- range(sapply(dens, function(d) range(d$y)))
cols<- c("#1f4e79", "#c0392b", "#2e7d32", "#c87f0a")[seq_along(dens)]
plot(NA, xlim = xr, ylim = yr, xlab = "tau (judgement noise)", ylab = "posterior density",
     main = "Posterior of tau under different priors")
for (i in seq_along(dens)) lines(dens[[i]], col = cols[i], lwd = 2)
abline(v = 0.3, lty = 3)   # the Monte Carlo sigma, for reference
legend("topright", c(names(dens), "MC sigma = 0.3"),
       col = c(cols, "black"), lwd = c(rep(2, length(dens)), 1),
       lty = c(rep(1, length(dens)), 3), bty = "n")
dev.off()

cat("\n11b_prior_sensitivity.R complete.\n")
cat("Tables: prior_sensitivity_tau.csv, prior_sensitivity_weights.csv\n")
cat("Figure: prior_sensitivity_tau.png\n")
