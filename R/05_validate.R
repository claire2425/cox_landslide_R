# R/05_validate.R
# Step 5: Validate the hazard map against the landslide inventory.
#
# Two checks:
#   (1) Lift: do inventory points concentrate in High/Very-High classes,
#       relative to the background pixel distribution?
#   (2) AUC: presence (inventory) vs background-absence ROC on the continuous
#       hazard score -- a standard, single-number accuracy statistic.


cat("\n--- 05_validate: cross-check inventory against hazard ---\n")

suppressPackageStartupMessages({ library(terra); library(sf) })
if (!exists("ahp_weights")) source("R/00_functions.R")

if (!file.exists("output/hazard_class.rds")) stop("Run R/04_hazard_overlay.R first.")
if (!file.exists("output/inventory.rds"))    stop("Run R/02_load_data.R first.")

hazard_class<- unwrap(readRDS("output/hazard_class.rds"))
hazard_cont<- unwrap(readRDS("output/hazard_continuous.rds"))
inventory<- readRDS("output/inventory.rds")

# ---- Inventory -> points -----------------------------------------------------
gt<- as.character(unique(st_geometry_type(inventory)))[1]
inv_pts<- if (gt %in% c("POLYGON","MULTIPOLYGON")) st_centroid(inventory) else inventory
inv_v<- vect(inv_pts)

# ---- (1) Lift analysis ------------------------------------------------------
sampled<- extract(hazard_class, inv_v, ID = FALSE)[, 1]
n_na<- sum(is.na(sampled))
cat(sprintf("Sampled %d inventory points (%d in NA/Restricted cells).\n",
            length(sampled), n_na))

labs<- c("Very Low","Low","Moderate","High","Very High")
inv_dist<- table(factor(sampled, levels = 1:5, labels = labs))
bg_vals<- values(hazard_class, na.rm = TRUE)
bg_dist<- table(factor(bg_vals, levels = 1:5, labels = labs))

inv_pct<- 100 * as.numeric(inv_dist) / sum(as.numeric(inv_dist))
bg_pct<- 100 * as.numeric(bg_dist)  / sum(as.numeric(bg_dist))
# Guard divide-by-zero: lift = NA where background class is empty.
lift<- ifelse(bg_pct > 0, round(inv_pct / bg_pct, 2), NA_real_)

validation<- data.frame(
  class = labs,
  inventory_n   = as.integer(inv_dist),
  inventory_pct = round(inv_pct, 2),
  background_pct = round(bg_pct, 2),
  lift = lift
)
cat("\nValidation (lift) summary:\n"); print(validation)

top2_inv<- sum(validation$inventory_pct[4:5])
top2_bg<- sum(validation$background_pct[4:5])
lift_top2<- if (top2_bg > 0) top2_inv / top2_bg else NA
cat(sprintf("\n%.1f%% of inventory in High/VeryHigh vs %.1f%% of background. Lift = %.2fx\n",
            top2_inv, top2_bg, lift_top2))
cat(if (!is.na(lift_top2) && lift_top2 > 1.5)
  "Lift check: PASSES.\n" else "Lift check: WEAK -- review reclass mappings.\n")

# ---- (2) AUC (presence vs background) ----
# I use AUC partly as an accuracy number and partly as a direction check. Each factor is scored 1-5, and it's
# easy to code the causal direction backwards (e.g. steeper slope should be higher hazard; sparse vegetation /
# low NDVI should be higher hazard). If I accidentally invert a rule, nothing errors -- the logic is just wrong.
# So AUC > 0.5 means the map gives higher scores where landslides actually happened (direction is fine)
# AUC = 0.5 is no better than chance; AUC < 0.5 means at least one rule is applied backwards and I need to fix it.

auc_val<- NA
if (requireNamespace("pROC", quietly = TRUE)) {
  
  set.seed(2026)
  
  #1. Extract hazard scores for presence (landslide) locations
  pres <- extract(hazard_cont, inv_v, ID = FALSE)[, 1]
  pres <- pres[!is.na(pres)]
  
  #2. Define the number of background samples
  n_bg <- max(500, length(pres) * 10)
  
  #3. Create pure absence areas by masking out landslide locations
  inv_full_v <- terra::vect(inventory) 
  hazard_bg_only <- terra::mask(hazard_cont, inv_full_v, inverse = TRUE)
  
  #4. Draw background samples from the clean map (excluding landslide areas)
  bg_s <- spatSample(hazard_bg_only, n_bg, method = "random", na.rm = TRUE)[, 1]
  
  y<- c(rep(1, length(pres)), rep(0, length(bg_s)))
  x<- c(pres, bg_s)
  roc_obj<- pROC::roc(y, x, quiet = TRUE)
  auc_val<- as.numeric(pROC::auc(roc_obj))
  
  cat(sprintf("\nAUC (presence vs background) = %.3f\n", auc_val))
  
  # Automated sanity check and diagnostic feedback
  if (auc_val >= 0.8) {
    cat("  ✅ [PASS] Excellent predictive capacity. Hazard direction is correctly aligned.\n")
  } else if (auc_val >= 0.7) {
    cat("  ✅ [PASS] Acceptable predictive capacity. Hazard direction is correctly aligned.\n")
  } else if (auc_val >= 0.5) {
    cat("  ⚠️ [WARN] Weak predictive capacity, but the overall direction is not inverted.\n")
  } else {
    cat("  ❌ [FATAL ERROR] AUC < 0.5! Severe warning: your reclassification direction is inverted!\n")
    cat("     (e.g., categorising dense NDVI or flat terrain as high risk)\n")
    cat("     Please return to 03_reclassify.R immediately to check and invert the 1-5 mappings of the erroneous factor.\n")
  }
  
  if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
  png("output/figures/roc_baseline.png", width = 700, height = 700)
  plot(roc_obj, main = sprintf("Baseline ROC  (AUC = %.3f)", auc_val))
  dev.off()
  cat("Saved: output/figures/roc_baseline.png\n")
} else {
  cat("\n(pROC not installed -> skipping AUC. install.packages('pROC') to enable.)\n")
}

# ---- Save -------------------------------------------------------------------
write.csv(validation, "output/tables/validation_summary.csv", row.names = FALSE)
saveRDS(list(validation = validation, lift_top2 = lift_top2, auc = auc_val),
        "output/validation.rds")
cat("\nSaved: output/tables/validation_summary.csv, output/validation.rds\n\n")

