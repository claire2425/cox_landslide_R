# R/09_monte_carlo.R
# Step 9: Monte Carlo weight perturbation (RQ2).
#
# What it does: Instead of moving one thing at a time like OAT, let's shake up 
# ALL the judgements at once. I'm adding log-normal noise to the matrix 2000 times 
# per sigma level (0.2, 0.3, 0.4) to see how wide the weight distributions get.
# Also tracking how many matrices stay consistent (CR < 0.10).

cat("\n--- 09_monte_carlo: joint uncertainty of AHP weights ---\n")

if (!exists("ahp_weights")) source("R/00_functions.R")

## --- load baseline submatrix -----------------------------------------------
rds_path<- "output/ahp_5factor_baseline.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path, " -- run R/01b_build_ahp_matrix.R first.")
}
base<- readRDS(rds_path)
A0<- base$matrix
w0<- ahp_weights(A0)
n<- nrow(A0)
fac<- rownames(A0)

cat("\nBaseline weights:\n"); print(round(w0, 4))
cat(sprintf("Baseline CR = %.4f\n", ahp_cr(A0)$CR))

## --- Monte Carlo settings --------------------------------------------------
N<- 2000          # draws per sigma
sigmas<- c(0.2, 0.3, 0.4)
main_sig<- 0.3           # the sigma we report in the text
SEED<- 2026

## one Monte Carlo draw: perturb every upper-triangle entry, keep reciprocity
mc_draw<- function(A, sigma) {
  Ap<- A
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    Ap[i, j]<- A[i, j] * exp(rnorm(1, 0, sigma))
    Ap[j, i]<- 1 / Ap[i, j]
  }
  list(w = ahp_weights(Ap), cr = ahp_cr(Ap)$CR)
}

## --- run all sigmas --------------------------------------------------------
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

interval_rows<- list()
mc_all_draws<- list()       # create a list to store all sigma results

for (sigma in sigmas) {
  set.seed(SEED)               # same seed per sigma -> reproducible
  W<- matrix(NA_real_, N, n, dimnames = list(NULL, fac))
  CR<- numeric(N)
  for (d in 1:N) {
    out<- mc_draw(A0, sigma)
    W[d, ]<- out$w
    CR[d]<- out$cr
  }
  cr_ok<- mean(CR < 0.10)
  keep<- CR < 0.10                 # the consistent draws
  W_keep<- W[keep, , drop = FALSE]
  
  cat(sprintf("\nsigma = %.1f : CR<0.10 in %.1f%% of %d draws (%d kept)\n",
              sigma, 100 * cr_ok, N, sum(keep)))
  
  for (method in c("all", "cr_pass")) {
    Wm<- if (method == "all") W else W_keep
    for (f in fac) {
      q<- quantile(Wm[, f], c(0.05, 0.50, 0.95))
      interval_rows[[length(interval_rows) + 1]]<- data.frame(
        sigma     = sigma,
        method    = method,
        factor    = f,
        baseline  = round(w0[f], 4),
        median    = round(q[2], 4),
        q05       = round(q[1], 4),
        q95       = round(q[3], 4),
        width90   = round(q[3] - q[1], 4),
        cr_ok_pct = round(100 * cr_ok, 1),
        n_draws   = nrow(Wm),
        row.names = NULL, stringsAsFactors = FALSE)
    }
  }
  
  # store the complete results for each sigma in the list for subsequent SMAA use
  mc_all_draws[[as.character(sigma)]]<- list(
    all = W, 
    cr_pass = W_keep
  )
  
  if (sigma == main_sig) {
    main_weights<- W
    main_weights_keep<- W_keep
  }
}

mc_intervals<- do.call(rbind, interval_rows)
write.csv(mc_intervals, "output/tables/mc_intervals.csv", row.names = FALSE)

# save the complete list containing all sigmas
saveRDS(mc_all_draws, "output/mc_weights.rds")
cat("\nSaved: output/mc_weights.rds (contains all sigmas), output/tables/mc_intervals.csv\n")

## density plot at all sigmas: all draws (solid) vs CR<0.10 (dashed)
## Multi-panel density plot comparing all sigmas
png("output/figures/mc_weight_densities_comparison.png", width = 1000, height = 1500, res = 130)

op<- par(mfrow = c(3, 1), mar = c(3, 4, 3, 1), oma = c(2, 0, 0, 0))

cols<- c("#4A7C89", "#DFA350", "#639A67", "#B96A4B", "#8675A9")

all_x<- c()
all_y<- c()
for (s in as.character(sigmas)) {
  for (f in fac) {
    all_x<- c(all_x, density(mc_all_draws[[s]]$all[, f])$x)
    all_y<- c(all_y, density(mc_all_draws[[s]]$all[, f])$y)
  }
}
global_xlim<- range(all_x)
global_ylim<- c(0, max(all_y) * 1.15) 

for (i in seq_along(sigmas)) {
  s<- as.character(sigmas[i])
  W_all<- mc_all_draws[[s]]$all
  W_keep<- mc_all_draws[[s]]$cr_pass
  cr_pass_rate<- 100 * nrow(W_keep) / N
  
  d_all<- lapply(fac, function(f) density(W_all[, f]))
  d_keep<- lapply(fac, function(f) density(W_keep[, f]))
  
  plot(NA, xlim = global_xlim, ylim = global_ylim, 
       xlab = "", ylab = "Density", bty = "l",
       main = sprintf("sigma = %.1f  (CR<0.10 pass rate: %.1f%%)", sigmas[i], cr_pass_rate))
  grid(nx = NULL, ny = NULL, col = "#E0E0E0", lty = "dotted")
  
  for (k in seq_along(fac)) {
    polygon(d_all[[k]], col = adjustcolor(cols[k], alpha.f = 0.15), border = NA) 
    lines(d_all[[k]],  col = adjustcolor(cols[k], alpha.f = 0.6), lwd = 2, lty = 1)            
    lines(d_keep[[k]], col = cols[k], lwd = 3.5, lty = 5) 
    abline(v = w0[fac[k]], col = cols[k], lty = 3, lwd = 1.5)        
    
    if (i == 1) {
      text(x = w0[fac[k]], y = max(d_all[[k]]$y) * 1.02, 
           labels = sprintf("%.3f", w0[fac[k]]), 
           col = cols[k], pos = 4, cex = 0.9, font = 2)       
    }
  }
  
  if (i == 1) {
    legend("topright", legend = fac, col = cols, lwd = 3.5, bty = "n", cex = 1)
    legend("right", legend = c("all draws", "CR<0.10 only", "baseline"),
           lty = c(1, 5, 3), col = "grey50", bty = "n", cex = 0.9, lwd = c(2, 3.5, 1.5))
  }
}

mtext("Weight", side = 1, outer = TRUE, line = 0.5, font = 2)

par(op); dev.off()
cat("Saved: output/figures/mc_weight_densities_comparison.png\n")
