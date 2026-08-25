# R/11c_sample_size_check.R
# Does the number of posterior draws change the answer?
# Re-fits the SAME Stan model at 5000 / 8000 / 10000 retained draws
# (4 chains, 1000 warmup each; only the post-warmup iterations change),
# and compares tau. This answers Katerina's "would 10000 be different?".
#
# Run from the project root (cox_landslide_R), after 01b has produced
# output/ahp_5factor_baseline.rds. Needs rstan (same as 11_bayesian_ahp.R).

cat("\n--- 11c: posterior draw-count sensitivity (5000 / 8000 / 10000) ---\n")
suppressMessages(library(rstan))
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 1)

STAN_FILE <- "R/stan/bayes_ahp.stan"
CHAINS    <- 4
WARMUP    <- 1000
SEED      <- 2026          # same seed as the main script

# --- load the same judgements the main script uses ---
ahp <- readRDS("output/ahp_5factor_baseline.rds")
A   <- ahp$matrix
n   <- nrow(A)
idx <- which(upper.tri(A), arr.ind = TRUE)     # (i, j), i < j
y   <- log(A[upper.tri(A)])                     # observed log-ratios
stan_data <- list(n = n, J = length(y), ii = idx[, 1], jj = idx[, 2], y = y)

# --- target TOTAL retained draws -> per-chain post-warmup iterations ---
targets   <- c(5000, 8000, 10000)
per_chain <- targets / CHAINS                   # 1250, 2000, 2500

rows <- lapply(seq_along(targets), function(k) {
  fit <- stan(file = STAN_FILE, data = stan_data,
              chains = CHAINS, warmup = WARMUP,
              iter = WARMUP + per_chain[k],
              seed = SEED, refresh = 0, cores = 1)
  s   <- summary(fit, pars = "tau")$summary      # mean, se_mean, sd, 2.5%, 50%, 97.5%, n_eff, Rhat
  tau <- as.numeric(rstan::extract(fit, "tau")$tau)
  data.frame(
    draws      = targets[k],
    tau_median = round(median(tau), 4),
    cri_lo     = round(quantile(tau, 0.025), 4),
    cri_hi     = round(quantile(tau, 0.975), 4),
    tau_mean   = round(s[, "mean"], 4),
    mcse_mean  = round(s[, "se_mean"], 4),        # Monte Carlo standard error of the mean
    n_eff      = round(s[, "n_eff"]),
    Rhat       = round(s[, "Rhat"], 4)
  )
})
tab <- do.call(rbind, rows)

cat("\nDraw-count comparison for tau:\n")
print(tab, row.names = FALSE)

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(tab, "output/tables/tau_drawcount_check.csv", row.names = FALSE)
cat("\nSaved -> output/tables/tau_drawcount_check.csv\n")
cat("How to read it: tau_median and the 95% CrI should barely move across the\n")
cat("three rows, and mcse_mean (the Monte Carlo error) should be tiny. That is\n")
cat("the evidence that going to 10000 would not change the conclusions.\n")
