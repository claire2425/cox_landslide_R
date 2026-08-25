# R/10_bootstrap.R
# Step 10: Bootstrap the landslide inventory to put a confidence interval on
#          the validation results (answers RQ2, part 2).
#
# Monte Carlo (09) handled WEIGHT uncertainty -- how much the weights move when
# the expert judgements are uncertain. Bootstrap handles a DIFFERENT, separate
# uncertainty: SAMPLING uncertainty. The baseline AUC (0.646) was computed from
# one particular set of 356 inventory points. If a slightly different set of
# landslides had been recorded, would the AUC be roughly the same? Bootstrap
# answers that by resampling the 356 points with replacement many times and
# recomputing the AUC (and the top-2 lift) each time, against the SAME fixed
# baseline susceptibility map.
#
# This is a genuine bootstrap: the 356 points are a real sample, so resampling
# them is legitimate. (We do NOT bootstrap the AHP matrix -- there is only one,
# and perturbing it is what Monte Carlo already does.)
#
# Design choice: we resample only the PRESENCE points (the 356 landslides). The
# background sample is drawn once and held fixed, so the resulting interval
# reflects uncertainty in the inventory itself, not in the background draw.
#
# Reads:  output/hazard_continuous.rds   (the baseline hazard, from 04)
#         the landslide inventory         (same object 05_validate uses)
# Writes: output/tables/bootstrap_validation.csv   (AUC + lift summary)
#         output/figures/bootstrap_auc.png         (histogram of bootstrap AUC)
#
# Needs terra + sf + pROC (GIS). Seeded, so re-running gives identical numbers.

cat("\n--- 10_bootstrap: confidence interval on validation (AUC, lift) ---\n")

if (!exists("reclass_codes")) source("R/00_functions.R")
suppressPackageStartupMessages({ library(terra); library(sf) })
stopifnot(requireNamespace("pROC", quietly = TRUE))

## --- load the baseline hazard ----------------------------------------------
haz_path<- "output/hazard_continuous.rds"
if (!file.exists(haz_path)) {
  stop("Missing ", haz_path, " -- run R/04_hazard_overlay.R first.")
}
hazard_cont<- unwrap(readRDS(haz_path))

## --- load the inventory (mirror how 05_validate gets it) -------------------
# 05_validate uses an object `inventory` (sf) and `inv_v` (terra vect). We rebuild
# them here from the saved inventory so this script can run on its own.
inv_path<- "output/inventory.rds"
if (file.exists(inv_path)) {
  inventory<- readRDS(inv_path)
} else {
  stop("Missing ", inv_path,
       " -- run R/02_load_data.R first (it saves the inventory).")
}
inv_v<- terra::vect(inventory)

## --- load the pipeline's 5-class hazard (so lift uses the SAME definition) --
# IMPORTANT: lift must use the same 5-class split as 04/05, not a fresh one cut
# from the background sample. 04 builds the classes from the WHOLE-MAP quantiles
# and saves them in hazard_class.rds; 05 computes lift against the whole-map
# background %. We reuse exactly that here so the bootstrap lift is directly
# comparable to the baseline lift table, instead of a slightly different number.
hc_path<- "output/hazard_class.rds"
if (!file.exists(hc_path)) {
  stop("Missing ", hc_path, " -- run R/04_hazard_overlay.R first.")
}
hazard_class<- unwrap(readRDS(hc_path))   # quantile 5-class raster from 04

## --- presence (fixed set) and a fixed background ---------------------------
set.seed(2026)

# presence: continuous score (for AUC) AND class label (for lift) at each point
pres_cont<- terra::extract(hazard_cont,  inv_v, ID = FALSE)[, 1]
pres_class<- terra::extract(hazard_class, inv_v, ID = FALSE)[, 1]
ok<- !is.na(pres_cont) & !is.na(pres_class)
pres_cont<- pres_cont[ok]
pres_class<- as.integer(pres_class[ok])
n_pres<- length(pres_cont)
cat(sprintf("Presence points with a hazard score: %d\n", n_pres))

# fixed background for AUC: continuous scores, landslides masked out
inv_full_v<- terra::vect(inventory)
hazard_bg_only<- terra::mask(hazard_cont, inv_full_v, inverse = TRUE)
n_bg<- max(500, n_pres * 10)
bg<- terra::spatSample(hazard_bg_only, n_bg, method = "random",
                          na.rm = TRUE)[, 1]
cat(sprintf("Fixed background sample (for AUC): %d\n", length(bg)))

# lift denominator: whole-map background % in the top-2 classes (matches 05)
bg_class<- as.integer(values(hazard_class, na.rm = TRUE))
bg_top2<- mean(bg_class >= 4)                       # fixed, whole-map denominator
cat(sprintf("Whole-map background in top-2 classes: %.1f%%\n", 100 * bg_top2))

## --- helper: AUC and top-2 lift for a resampled presence index -------------
# AUC uses continuous scores; lift uses the pipeline's class labels and the
# whole-map background denominator -- exactly the 05_validate definition.
auc_lift<- function(idx) {
  y<- c(rep(1, length(idx)), rep(0, length(bg)))
  x<- c(pres_cont[idx], bg)
  auc<- as.numeric(pROC::auc(pROC::roc(y, x, quiet = TRUE)))
  pres_top2<- mean(pres_class[idx] >= 4)
  lift<- if (bg_top2 > 0) pres_top2 / bg_top2 else NA
  c(auc = auc, lift = lift)
}

## --- point estimate (no resampling) ----------------------------------------
obs<- auc_lift(seq_len(n_pres))
cat(sprintf("\nPoint estimate: AUC = %.3f, top-2 lift = %.3f\n", obs["auc"], obs["lift"]))

## --- bootstrap -------------------------------------------------------------
B<- 2000
boot<- matrix(NA_real_, B, 2, dimnames = list(NULL, c("auc", "lift")))
set.seed(2026)
for (b in 1:B) {
  idx<- sample.int(n_pres, n_pres, replace = TRUE)   # resample presence w/ replacement
  boot[b, ]<- auc_lift(idx)
}

## --- summarise: 95% percentile intervals -----------------------------------
summ<- function(v) c(mean = mean(v),
                      q025 = quantile(v, 0.025),
                      q975 = quantile(v, 0.975))
auc_ci<- summ(boot[, "auc"])
lift_ci<- summ(boot[, "lift"])

cat(sprintf("\nBootstrap (%d resamples):\n", B))
cat(sprintf("  AUC      mean %.3f, 95%% CI [%.3f, %.3f]\n",
            auc_ci[1], auc_ci[2], auc_ci[3]))
cat(sprintf("  top-2 lift mean %.3f, 95%% CI [%.3f, %.3f]\n",
            lift_ci[1], lift_ci[2], lift_ci[3]))
cat(sprintf("  Fraction of bootstrap AUCs above 0.5: %.1f%%\n",
            100 * mean(boot[, "auc"] > 0.5)))

## --- save table ------------------------------------------------------------
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
res<- data.frame(
  metric   = c("AUC", "top2_lift"),
  point    = round(c(obs["auc"], obs["lift"]), 4),
  boot_mean= round(c(auc_ci[1], lift_ci[1]), 4),
  ci_lo    = round(c(auc_ci[2], lift_ci[2]), 4),
  ci_hi    = round(c(auc_ci[3], lift_ci[3]), 4),
  row.names = NULL)
write.csv(res, "output/tables/bootstrap_validation.csv", row.names = FALSE)
cat("\nSaved: output/tables/bootstrap_validation.csv\n")

## --- histogram of bootstrap AUC --------------------------------------------
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
png("output/figures/bootstrap_auc.png", width = 1000, height = 650, res = 130)
op<- par(mar = c(4, 4, 3, 1))

col_bars<- adjustcolor("cadetblue3", alpha.f = 0.75)
col_pt<- "darkorange"
col_ci<- "gray40"

x_min<- 0.49
x_max<- max(boot[, "auc"])*1.025

hist(boot[, "auc"], breaks = 40, col = col_bars, border = "white",
     main = sprintf("Bootstrap AUC (B = %d)", B),
     xlab = "AUC", ylab = "Frequency",
     xlim = c(x_min, x_max))

abline(v = 0.5, col = "gray70", lwd = 1.5, lty = 3)         # Random baseline
abline(v = obs["auc"], col = col_pt, lwd = 3)               # Point estimate
abline(v = auc_ci[2], col = col_ci, lwd = 2.5, lty = 2)     # CI lower
abline(v = auc_ci[3], col = col_ci, lwd = 2.5, lty = 2)     # CI upper

legend("topleft", 
       legend = c(sprintf("Point Estimate (AUC = %.3f)", obs["auc"]), 
                  sprintf("95%% CI [%.3f, %.3f]", auc_ci[2], auc_ci[3]),
                  "Random Baseline (AUC = 0.5)"),
       col = c(col_pt, col_ci, "gray70"),
       lwd = c(3, 2.5, 1.5), lty = c(1, 2, 3), bty = "n", cex = 0.95)

par(op); dev.off()
cat("Saved: output/figures/bootstrap_auc.png\n")

