# R/05e_posterior_auc.R
# -------------------------------------------------------------------------
# Propagate the Bayesian posterior into the validation metric: recompute the
# AUC under each of the 8000 posterior weight draws, turning the single
# Ahmed-weighted point estimate (AUC 0.589) into an uncertainty-aware interval.
# This closes the loop between weight uncertainty and performance uncertainty
# (Conclusion, future work) with no new data.
#
# Method: extract the five reclassified factor values at the presence points
# and at a fixed 1:10 background sample (drawn once, seed 2026, so only the
# weights vary), then score both under every posterior draw and compute AUC.
# The background is held fixed across draws so the spread reflects weight
# uncertainty alone. AUC is the Mann-Whitney statistic (exact, fast).
#
# SAFETY: writes ONLY output/tables/posterior_auc.csv. Nothing else is touched.
# Run from the project root after step 11:  Rscript R/05e_posterior_auc.R
# -------------------------------------------------------------------------

cat("\n--- 05e_posterior_auc: AUC under the 8000 posterior weight draws ---\n")
suppressPackageStartupMessages({ library(terra); library(sf) })

W  <- readRDS("output/bayes_weights.rds")            # M x 5 posterior draws
rs <- unwrap(readRDS("output/reclass_stack.rds"))    # 5 reclassified factors
inventory <- readRDS("output/inventory.rds")
M  <- nrow(W)

# align posterior weight columns to the model's factor order
criteria <- c("Slope","LandCover","NDVI","SoilType","DistRoad")
if (is.null(colnames(W))) colnames(W) <- criteria
W <- W[, criteria, drop = FALSE]

# ---- presence: reclass values at inventory points ----
gt <- as.character(unique(st_geometry_type(inventory)))[1]
inv_v <- vect(if (gt %in% c("POLYGON","MULTIPOLYGON")) st_centroid(inventory) else inventory)
Xp <- extract(rs, inv_v, ID = FALSE)
Xp <- Xp[stats::complete.cases(Xp), criteria, drop = FALSE]

# ---- background: fixed 1:10 sample, drawn once, inventory masked out ----
set.seed(2026)
n_bg <- max(500, nrow(Xp) * 10)
rs_bg <- terra::mask(rs, terra::vect(inventory), inverse = TRUE)
Xb <- spatSample(rs_bg, n_bg, method = "random", na.rm = TRUE)[, criteria, drop = FALSE]
Xb <- Xb[stats::complete.cases(Xb), , drop = FALSE]

Xall <- as.matrix(rbind(Xp, Xb))
y    <- c(rep(1L, nrow(Xp)), rep(0L, nrow(Xb)))
n1 <- sum(y == 1); n0 <- sum(y == 0)
cat(sprintf("%d presence, %d background, %d posterior draws.\n", n1, n0, M))

fast_auc <- function(score) {
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- score under every posterior draw ----
aucs <- numeric(M)
pb <- txtProgressBar(min = 0, max = M, style = 3)
for (m in 1:M) { aucs[m] <- fast_auc(as.numeric(Xall %*% W[m, ])); setTxtProgressBar(pb, m) }
close(pb)

qs <- quantile(aucs, c(0.025, 0.5, 0.975))
cat(sprintf("\nPosterior AUC:  mean %.3f  median %.3f  95%% interval [%.3f, %.3f]\n",
            mean(aucs), qs[2], qs[1], qs[3]))
cat(sprintf("Share of draws with AUC > 0.5: %.1f%%   (Ahmed point estimate: 0.589)\n",
            100 * mean(aucs > 0.5)))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(data.frame(metric = c("mean","median","q2.5","q97.5","pct_above_0.5"),
                     value  = round(c(mean(aucs), qs[2], qs[1], qs[3], 100*mean(aucs>0.5)), 3)),
          "output/tables/posterior_auc.csv", row.names = FALSE)
cat("Saved ONLY: output/tables/posterior_auc.csv\n")
