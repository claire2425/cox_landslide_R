## ============================================================
## tail-ESS check for tau  (SELF-CONTAINED -- re-fits the baseline model)
## Run from the project root (the folder that contains output/ and R/stan/),
## exactly like 11c_sample_size_check.R. No pre-existing `fit` needed.
##
## Reproduces the reported baseline fit: 4 chains x 2000 post-warmup = 8000
## draws, seed 2026, using R/stan/bayes_ahp.stan (half-normal(1.0) on tau).
## (You may see "hash mismatch so recompiling" -- harmless, as before.)
## ============================================================

STAN_FILE   <- "R/stan/bayes_ahp.stan"
STAN_CHAINS <- 4
STAN_WARMUP <- 1000
STAN_ITER   <- 2000          # post-warmup per chain -> 8000 draws total

## 1. Load the judgements (identical to 11_bayesian_ahp.R) --------------------
ahp <- readRDS("output/ahp_5factor_baseline.rds")
A   <- ahp$matrix
n   <- nrow(A)
idx <- which(upper.tri(A), arr.ind = TRUE)   # rows = (i, j), i < j
y   <- log(A[upper.tri(A)])                  # observed log-ratios
stan_data <- list(n = n, J = length(y), ii = idx[, 1], jj = idx[, 2], y = y)

## 2. Fit (same settings as the reported baseline) ---------------------------
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
fit <- stan(file = STAN_FILE, data = stan_data,
            chains = STAN_CHAINS, warmup = STAN_WARMUP,
            iter = STAN_WARMUP + STAN_ITER, seed = 2026, refresh = 0)

## 3. bulk-ESS and tail-ESS for tau ------------------------------------------
if (requireNamespace("posterior", quietly = TRUE)) {
  library(posterior)
  tab <- summarise_draws(
    subset_draws(as_draws_df(fit), variable = "tau"),
    median,
    q2.5  = ~ quantile(.x, 0.025),
    q97.5 = ~ quantile(.x, 0.975),
    rhat, ess_bulk, ess_tail)
  print(tab)
} else {
  ## fallback: rstan's monitor() also reports Bulk_ESS and Tail_ESS
  print(rstan::monitor(as.array(fit, pars = "tau"), print = FALSE)[, c("Rhat", "Bulk_ESS", "Tail_ESS")])
}
cat("\n--> report ess_tail for tau; Vehtari et al. (2021) recommend ESS > 400.\n")
