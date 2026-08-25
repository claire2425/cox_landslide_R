# R/05d_ahmed_auc.R
# -------------------------------------------------------------------------
# Non-destructive check: AUC and lift for the ALTERNATIVE Ahmed-weighted
# baseline map. This lets the confidence layer (which perturbs the Ahmed
# submatrix) be reported against the accuracy of the map it actually certifies,
# instead of only against the WO1 map.
#
# SAFETY: this script writes ONLY to output/tables/ahmed_weighted_auc.csv.
# It does NOT overwrite ahp_weights.rds, hazard_*.rds/.tif, validation.rds,
# or any other WO1 output. Your existing WO1 results are left untouched.
#
# It reproduces 05_validate.R's AUC protocol exactly (same seed, same 1:10
# presence-to-background sampling, same pROC call), so the number is directly
# comparable to your WO1 baseline AUC of 0.646.
#
# Run from the project root, after step 03 has been run at least once:
#     Rscript R/05d_ahmed_auc.R
# or in RStudio:  source("R/05d_ahmed_auc.R")
# -------------------------------------------------------------------------

cat("\n--- 05d_ahmed_auc: Ahmed-weighted baseline AUC (non-destructive) ---\n")
suppressPackageStartupMessages({ library(terra); library(sf) })

if (!file.exists("output/reclass_stack.rds"))
  stop("Missing output/reclass_stack.rds -- run R/03_reclassify.R first.")
if (!file.exists("output/inventory.rds"))
  stop("Missing output/inventory.rds -- run R/02_load_data.R first.")

reclass_stack <- unwrap(readRDS("output/reclass_stack.rds"))
inventory     <- readRDS("output/inventory.rds")

## ---- Ahmed 5x5 submatrix: identical to 01_ahp_pairwise.R --------------------
criteria <- c("Slope", "LandCover", "NDVI", "SoilType", "DistRoad")
A_ahmed <- matrix(c(
  1,    2,    6,    1/5,  1/5,
  1/2,  1,    5,    1/4,  1/6,
  1/6,  1/5,  1,    1/6,  1/9,
  5,    4,    6,    1,    1/2,
  5,    6,    9,    2,    1
), nrow = 5, byrow = TRUE,
   dimnames = list(c("LandCover", "NDVI", "DistRoad", "Slope", "SoilType"),
                   c("LandCover", "NDVI", "DistRoad", "Slope", "SoilType")))
A_ahmed <- A_ahmed[criteria, criteria]

## principal-eigenvector weights (base R, matches ahp_weights())
ev  <- eigen(A_ahmed)
idx <- which.max(Re(ev$values))
w_a <- abs(Re(ev$vectors[, idx])); w_a <- w_a / sum(w_a)
names(w_a) <- criteria
cat("Ahmed-weighted overlay weights:\n"); print(round(w_a, 4))

## align to the stack's layer order BY NAME (exactly as 04_hazard_overlay.R does)
if (!all(names(reclass_stack) %in% names(w_a)))
  stop("Weight names do not cover all stack layers.")
w_a <- w_a[names(reclass_stack)]

## ---- Weighted overlay, in memory only; equal-width 5 classes over [1,3] -----
hazard_cont <- terra::app(reclass_stack, fun = function(v) sum(v * w_a))
names(hazard_cont) <- "Hazard_continuous_ahmed"
breaks_ew   <- seq(1, 3, length.out = 6)       # same fixed cuts as the baseline
hazard_class <- classify(hazard_cont,
  cbind(breaks_ew[-6], breaks_ew[-1], 1:5), include.lowest = TRUE)

## ---- Inventory -> points ---------------------------------------------------
gt <- as.character(unique(st_geometry_type(inventory)))[1]
inv_pts <- if (gt %in% c("POLYGON", "MULTIPOLYGON")) st_centroid(inventory) else inventory
inv_v   <- vect(inv_pts)

## ---- (1) Lift (same construction as 05_validate.R) -------------------------
labs    <- c("Very Low", "Low", "Moderate", "High", "Very High")
sampled <- extract(hazard_class, inv_v, ID = FALSE)[, 1]
inv_ct  <- as.numeric(table(factor(sampled, levels = 1:5)))
bg_ct   <- as.numeric(table(factor(values(hazard_class, na.rm = TRUE), levels = 1:5)))
inv_pct <- 100 * inv_ct / sum(inv_ct)
bg_pct  <- 100 * bg_ct  / sum(bg_ct)
lift    <- ifelse(bg_pct > 0, round(inv_pct / bg_pct, 2), NA_real_)
top2_lift <- sum(inv_pct[4:5]) / sum(bg_pct[4:5])

## ---- (2) AUC (presence vs background) -- identical protocol to 05_validate.R
auc_val <- NA
if (requireNamespace("pROC", quietly = TRUE)) {
  set.seed(2026)
  pres <- extract(hazard_cont, inv_v, ID = FALSE)[, 1]; pres <- pres[!is.na(pres)]
  n_bg <- max(500, length(pres) * 10)
  hazard_bg_only <- terra::mask(hazard_cont, terra::vect(inventory), inverse = TRUE)
  bg_s <- spatSample(hazard_bg_only, n_bg, method = "random", na.rm = TRUE)[, 1]
  y <- c(rep(1, length(pres)), rep(0, length(bg_s)))
  x <- c(pres, bg_s)
  auc_val <- as.numeric(pROC::auc(pROC::roc(y, x, quiet = TRUE)))
} else {
  cat("(pROC not installed -> install.packages('pROC') to get the AUC.)\n")
}

## ---- Report + save ONLY the new file ---------------------------------------
cat(sprintf("\nAhmed-weighted map:  AUC = %.3f   |   top-2 (High/VeryHigh) lift = %.2f\n",
            auc_val, top2_lift))
cat("Compare with WO1 baseline:  AUC = 0.646,  top-2 lift = 1.68\n")
cat("\nPer-class lift (Ahmed-weighted):\n")
print(setNames(round(lift, 2), labs))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
out <- data.frame(class = labs,
                  inventory_pct = round(inv_pct, 2),
                  background_pct = round(bg_pct, 2),
                  lift = lift)
out <- rbind(out, data.frame(class = "AUC", inventory_pct = NA,
                             background_pct = NA, lift = round(auc_val, 3)))
write.csv(out, "output/tables/ahmed_weighted_auc.csv", row.names = FALSE)
cat("\nSaved ONLY: output/tables/ahmed_weighted_auc.csv  (no other file touched).\n\n")
