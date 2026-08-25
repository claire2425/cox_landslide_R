# R/04_hazard_overlay.R
# Step 4: Weighted overlay -> baseline landslide susceptibility map.
#
# Hazard(x) = sum_j  weight_j * reclass_layer_j(x)
# Inputs are 1/2/3 (NA = Restricted), so Hazard lies in [1, 3].

cat("\n--- 04_hazard_overlay: compute baseline hazard map ---\n")

suppressPackageStartupMessages({ library(terra) })
if (!exists("ahp_weights")) source("R/00_functions.R")

if (!file.exists("output/reclass_stack.rds")) stop("Run R/03_reclassify.R first.")
if (!file.exists("output/ahp_weights.rds"))   stop("Run R/01_ahp_pairwise.R first.")

reclass_stack<- unwrap(readRDS("output/reclass_stack.rds"))
weights<- readRDS("output/ahp_weights.rds")

# ---- Re-order weights to the stack layer order (robust) ---------------------
if (!all(names(reclass_stack) %in% names(weights)))
  stop("Weight names do not cover all stack layers:\n  stack: ",
       paste(names(reclass_stack), collapse=", "),
       "\n  weights: ", paste(names(weights), collapse=", "))
weights<- weights[names(reclass_stack)]   # align by NAME, not position
stopifnot(abs(sum(weights) - 1) < 1e-6)
cat("Weights (aligned to layer order):\n"); print(round(weights, 3))

# ---- Weighted overlay (NA propagates = Restricted) --------------------------
cat("\nComputing weighted overlay...\n")
hazard<- terra::app(
  reclass_stack,
  fun = function(v) sum(v * weights)   # v is a vector of the 5 layer values
)
names(hazard)<- "Hazard_continuous"

v<- values(hazard, na.rm = TRUE)
cat(sprintf("Continuous hazard: range=[%.3f, %.3f]  mean=%.3f  median=%.3f  N=%d\n",
            min(v), max(v), mean(v), median(v), length(v)))

# ---- Classify: equal-width 5 bands over [1,3] -------------------------------
breaks_ew<- seq(1, 3, length.out = 6)
hazard_ew<- classify(hazard,
  cbind(breaks_ew[-6], breaks_ew[-1], 1:5), include.lowest = TRUE)
names(hazard_ew)<- "Hazard_class_equalwidth"

# ---- Classify: quantile 5 classes (equal pixel counts) ----------------------
qb<- quantile(v, probs = seq(0, 1, 0.2), na.rm = TRUE)
qb[1]<- qb[1] - 1e-9                      # ensure lowest value is included
qb<- unique(qb)                           # guard against tied quantiles
if (length(qb) >= 2) {
  hazard_q<- classify(hazard,
    cbind(qb[-length(qb)], qb[-1], seq_len(length(qb)-1)), include.lowest = TRUE)
  names(hazard_q)<- "Hazard_class_quantile"
} else {
  hazard_q<- hazard_ew
  cat("  (quantile breaks degenerate; falling back to equal-width)\n")
}

# Baseline used downstream = equal-width (keeps "near 3 = high" meaning).
hazard_class<- hazard_ew
names(hazard_class)<- "Hazard_class"

# ---- Class summary ----------------------------------------------------------
labs<- c("Very Low","Low","Moderate","High","Very High")
tab<- table(factor(values(hazard_class, na.rm = TRUE), levels = 1:5))
class_summary<- data.frame(
  class = 1:5, label = labs,
  lower = round(breaks_ew[-6], 3), upper = round(breaks_ew[-1], 3),
  pixels = as.integer(tab),
  pct = round(100 * as.integer(tab) / sum(as.integer(tab)), 2)
)
cat("\nHazard class distribution (equal-width baseline):\n"); print(class_summary)

# ---- Save -------------------------------------------------------------------
saveRDS(wrap(hazard),       "output/hazard_continuous.rds")
saveRDS(wrap(hazard_class), "output/hazard_class.rds")
saveRDS(wrap(hazard_q),     "output/hazard_class_quantile.rds")
# Save the exact breaks this baseline classification used, so the posterior (12) and SMAA (13) reuse the
# SAME equal-width cuts instead of recomputing quantiles on a possibly-different data state. That's what makes
# the uncertainty layer actually describe the uncertainty of the validated baseline map, not some other map.
saveRDS(as.numeric(breaks_ew), "output/hazard_breaks.rds")
writeRaster(hazard,       "output/hazard_continuous.tif", overwrite = TRUE)
writeRaster(hazard_class, "output/hazard_class.tif",      overwrite = TRUE)
writeRaster(hazard_q,     "output/hazard_class_quantile.tif", overwrite = TRUE)
write.csv(class_summary,  "output/tables/hazard_class_summary.csv", row.names = FALSE)

cat("\nSaved: hazard_continuous.tif, hazard_class.tif (equal-width baseline),\n")
cat("       hazard_class_quantile.tif, tables/hazard_class_summary.csv\n\n")

