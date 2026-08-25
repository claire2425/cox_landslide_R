# R/11_bayesian_ahp.R
# Step 11: Bayesian AHP (RQ3a) - getting the actual posterior weights.
#
# Why do this? MC assumes the noise (sigma = 0.3). Bayesian actually INFERS it (tau) 
# straight from how (in)consistent the expert's matrix is.
#
# The setup: Log-linear model. theta = log-weights. 
# I'm using Stan (NUTS) as the main engine because it's fast and gives good diagnostics.
# 
# Note: Don't trust the numbers until I check Rhat < 1.01 and n_eff!

cat("\n--- 11_bayesian_ahp: posterior over AHP weights ---\n")
set.seed(2026)   # reproducibility, matching the rest of the pipeline

## 0. Choice of engine + MCMC settings
ENGINE<- "stan"      # "stan" (recommended) or "mh" (base-R fallback)
STAN_FILE<- "R/stan/bayes_ahp.stan"

# Stan depth (used if ENGINE = "stan")
STAN_CHAINS<- 4
STAN_WARMUP<- 1000
STAN_ITER<- 2000   # post-warmup per chain -> 4 chains x 2000 = 8000 draws total

# MH depth (used if ENGINE = "mh"). These do NOT need supervisor sign-off --
# run, read the diagnostics block, and adjust:
#   - acceptance far from ~0.2-0.4  -> change PROP_SD
#   - sticky / non-stationary trace -> raise N_ITER and BURN_IN
#   - thinning only controls how many draws you store, not correctness
N_ITER<- 60000
BURN_IN<- 10000
THIN<- 10         # -> 5000 stored draws
PROP_SD<- 0.25       # random-walk proposal SD on (theta_2..n, log_tau).
                      # Starting guess only -- the acceptance rate depends on
                      # YOUR matrix. Run, read the printed acceptance rate, and
                      # nudge PROP_SD down (if acceptance < ~0.2) or up (if
                      # > ~0.4) until it lands in the 0.2-0.4 band.

## 1. Load the judgements
rds_path<- "output/ahp_5factor_baseline.rds"
if (!file.exists(rds_path))
  stop("Missing ", rds_path, " -- run R/01b_build_ahp_matrix.R first.")

ahp<- readRDS(rds_path)
A<- ahp$matrix                                    # 5x5 reciprocal pairwise matrix
stopifnot(is.matrix(A), nrow(A) == 5, ncol(A) == 5)

factors<- rownames(A)
if (is.null(factors)) factors<- c("Slope","LandCover","NDVI","SoilType","DistRoad")
n<- nrow(A)

idx<- which(upper.tri(A), arr.ind = TRUE)         # each row is (i, j), i < j
y<- log(A[upper.tri(A)])                        # observed log-ratios
ii<- idx[, 1]; jj<- idx[, 2]
J<- length(y)
cat(sprintf("Loaded %d upper-triangle judgements from a %dx%d matrix.\n", J, n, n))

# baseline (point-estimate) weights for the side-by-side comparison.
# Prefer the stored baseline; else principal eigenvector (NOT a uniform fallback).
if (!is.null(ahp$weights)) {
  w_base<- as.numeric(ahp$weights)
} else {
  ev<- eigen(A)$vectors[, 1]
  w_base<- Re(ev) / sum(Re(ev))
}
names(w_base)<- factors

softmax<- function(theta) { e<- exp(theta - max(theta)); e / sum(e) }

## 2. Decide engine (auto-fallback if rstan missing)
have_rstan<- requireNamespace("rstan", quietly = TRUE)
if (ENGINE == "stan" && !have_rstan) {
  cat("rstan not available -> falling back to the base-R Metropolis sampler.\n")
  ENGINE<- "mh"
}

## 3a. ENGINE = stan
if (ENGINE == "stan") {
  if (!file.exists(STAN_FILE))
    stop("Missing ", STAN_FILE, " -- place bayes_ahp.stan there.")
  cat("Engine: Stan (NUTS).\n")
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1)   # single-core: avoids the parallel-cluster unserialize error on this setup

  stan_data<- list(n = n, J = J, ii = ii, jj = jj, y = y)
  fit<- rstan::stan(
    file    = STAN_FILE, data = stan_data,
    chains  = STAN_CHAINS, warmup = STAN_WARMUP, iter = STAN_WARMUP + STAN_ITER,
    seed    = 2026, refresh = 500,
    cores   = 1                          # run the chains sequentially (matrix is tiny, so no speed loss)
  )

  # convergence evidence (report these in the write-up)
  print(rstan::summary(fit, pars = c("theta","tau","w"))$summary[, c("mean","sd","n_eff","Rhat")])
  cat("Check: all Rhat < 1.01 and n_eff healthy before trusting the posterior.\n")

  W<- rstan::extract(fit, "w")$w                 # draws x n
  colnames(W)<- factors
  tau<- as.numeric(rstan::extract(fit, "tau")$tau)
}

## 3b. ENGINE = mh  (fallback; identical model)
if (ENGINE == "mh") {
  cat("Engine: hand-written random-walk Metropolis-Hastings (base R).\n")

  theta_from_par<- function(par) c(0, par[1:(n - 1)])
  tau_from_par<- function(par) exp(par[n])       # log_tau -> tau > 0

  log_prior<- function(par) {
    th<- par[1:(n - 1)]                            # theta_2..n
    tau<- tau_from_par(par)
    lp_theta<- sum(dnorm(th, 0, 3, log = TRUE))     # theta_k ~ N(0, 3^2)
    # tau ~ Half-Normal(0,1) on tau>0, + Jacobian for log_tau
    lp_tau<- dnorm(tau, 0, 1, log = TRUE) + log(2) + par[n]
    lp_theta + lp_tau
  }
  log_lik<- function(par) {
    th<- theta_from_par(par)
    tau<- tau_from_par(par)
    mu<- th[ii] - th[jj]
    sum(dnorm(y, mu, tau, log = TRUE))
  }
  log_post<- function(par) log_prior(par) + log_lik(par)

  p<- n
  cur<- c(rep(0, n - 1), log(0.5))              # start: flat weights, tau=0.5
  cur_lp<- log_post(cur)
  keep_n<- floor((N_ITER - BURN_IN) / THIN)
  par_draws<- matrix(NA_real_, nrow = keep_n, ncol = p)
  accepts<- 0L; k<- 0L

  for (it in 1:N_ITER) {
    prop<- cur + rnorm(p, sd = PROP_SD)
    prop_lp<- log_post(prop)
    if (is.finite(prop_lp) && log(runif(1)) < (prop_lp - cur_lp)) {
      cur<- prop; cur_lp<- prop_lp; accepts<- accepts + 1L
    }
    if (it > BURN_IN && ((it - BURN_IN) %% THIN == 0)) { k<- k + 1L; par_draws[k, ]<- cur }
  }
  acc_rate<- accepts / N_ITER
  cat(sprintf("MCMC done. Acceptance rate = %.3f (aim ~0.20-0.40).\n", acc_rate))

  W<- matrix(NA_real_, nrow = keep_n, ncol = n, dimnames = list(NULL, factors))
  tau<- numeric(keep_n)
  for (m in 1:keep_n) {
    W[m, ]<- softmax(theta_from_par(par_draws[m, ]))
    tau[m]<- exp(par_draws[m, n])
  }
}

## 4. Save weight draws (identical format regardless of engine)
dir.create("output/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
saveRDS(W, "output/bayes_weights.rds")
M<- nrow(W)

## 5. Posterior summaries (median + 95% credible interval per weight)
tab<- data.frame(
  factor   = factors,
  baseline = round(as.numeric(w_base), 3),
  post_med = round(apply(W, 2, median), 3),
  cri_lo   = round(apply(W, 2, function(x) quantile(x, 0.025)), 3),
  cri_hi   = round(apply(W, 2, function(x) quantile(x, 0.975)), 3)
)
tab$cri_width<- round(tab$cri_hi - tab$cri_lo, 3)
write.csv(tab, "output/tables/bayes_intervals.csv", row.names = FALSE)

tau_summary<- data.frame(
  param = "tau",
  med   = round(median(tau), 3),
  lo    = round(quantile(tau, 0.025), 3),
  hi    = round(quantile(tau, 0.975), 3)
)
write.csv(tau_summary, "output/tables/bayes_tau.csv", row.names = FALSE)

cat("\nPosterior weight summary:\n"); print(tab)
cat(sprintf("\nJudgement-noise tau: median %.3f, 95%% CrI [%.3f, %.3f]\n",
            tau_summary$med, tau_summary$lo, tau_summary$hi))
# Interpretation to make AFTER you see real numbers: compare this tau against
# the sigma you assumed in Monte Carlo (09). If they differ, that gap is itself
# a finding about how large the expert judgement noise really is -- do not
# pre-judge the direction; let your run decide.

## 6. Figures: posterior densities + trace
png("output/figures/bayes_weight_posterior.png", width = 1500, height = 1000, res = 150)
op<- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (f in seq_len(n)) {
  plot(density(W[, f]), main = factors[f], xlab = "weight", lwd = 2)
  abline(v = w_base[f], lty = 3)
  abline(v = tab$post_med[f], col = "steelblue", lwd = 2)
}
plot.new(); legend("center", c("posterior","post. median","baseline"),
                   lwd = c(2,2,1), lty = c(1,1,3),
                   col = c("black","steelblue","black"), bty = "n")
par(op); dev.off()

png("output/figures/bayes_trace.png", width = 1500, height = 1000, res = 150)
op<- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (f in seq_len(n)) plot(W[, f], type = "l", main = factors[f],
                           xlab = "draw", ylab = "weight")
plot(tau, type = "l", main = "tau", xlab = "draw", ylab = "tau")
par(op); dev.off()

## 7. Diagnostics to read before trusting the run
cat("\n--- diagnostics ---\n")
if (ENGINE == "stan") {
  cat("Stan: rely on Rhat (<1.01) and n_eff printed above; check for any\n")
  cat("divergent transitions in the sampler messages.\n")
} else {
  ess_ar1<- function(x) {
    ac<- acf(x, plot = FALSE, lag.max = 1)$acf[2]
    length(x) * (1 - ac) / (1 + ac)
  }
  cat(sprintf("acceptance rate : %.3f  (retune PROP_SD if outside ~0.2-0.4)\n", acc_rate))
  cat("approx ESS per weight:\n"); print(round(apply(W, 2, ess_ar1)))
  cat("If trace plots wander or ESS < ~400, raise N_ITER/BURN_IN.\n")
}

cat(sprintf("\n11_bayesian_ahp.R complete  [engine: %s, %d draws].\n", ENGINE, M))
