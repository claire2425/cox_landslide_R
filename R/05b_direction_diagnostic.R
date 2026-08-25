# R/05b_direction_diagnostic.R
# PURPOSE: diagnose the low AUC (0.517) by testing whether any factor's
# reclass DIRECTION is reversed, then quantify the fix.
#
# How to use:
#   1. Run the normal pipeline first (01 -> 06) so output/ exists.
#   2. source("R/05b_direction_diagnostic.R")
#   3. Read the printed table. It re-runs the overlay+AUC with each factor
#      flipped one at a time and tells me which flip (if any) raises the AUC.
#
# Background: validation showed AUC = 0.517 and a U-shaped lift (landslide
# points appearing in BOTH Very-Low and Very-High zones). That pattern means
# a factor is scored backwards. WO1.txt scores high-NDVI as high-susceptibility
# (bin5 -> 3), but denser vegetation generally stabilises slopes -> should be
# LOW susceptibility. This script tests that hypothesis quantitatively.

cat("\n--- 05b_direction_diagnostic: find the reversed factor ---\n")

suppressPackageStartupMessages({ library(terra); library(sf) })
if (!exists("ahp_weights")) source("R/00_functions.R")
stopifnot(file.exists("output/reclass_stack.rds"),
          file.exists("output/ahp_weights.rds"),
          file.exists("output/inventory.rds"))

stack0<- unwrap(readRDS("output/reclass_stack.rds"))
weights<- readRDS("output/ahp_weights.rds")
weights<- weights[names(stack0)]
inv<- readRDS("output/inventory.rds")
gt<- as.character(unique(st_geometry_type(inv)))[1]
inv_pts<- if (gt %in% c("POLYGON","MULTIPOLYGON")) st_centroid(inv) else inv
inv_v<- vect(inv_pts)

# Convert the sf object to a terra vector outside the function for subsequent masking
inv_full_v<- terra::vect(inv)

# Helper: build hazard from a stack, compute AUC against pure background
auc_of<- function(stk) {
  haz<- terra::app(stk, fun = function(v) sum(v * weights))
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("pROC needed: install.packages('pROC')")
  }
  set.seed(2026)
  
  #1. Extract presence locations
  pres<- extract(haz, inv_v, ID = FALSE)[, 1]
  pres<- pres[!is.na(pres)]
  
  #2. Mask out actual landslide boundaries to ensure the purity of background sampling (Pure Absence)
  haz_bg_only<- terra::mask(haz, inv_full_v, inverse = TRUE)
  
  #3. Draw background samples
  bg<- spatSample(haz_bg_only, max(500, length(pres) * 10), method = "random", na.rm = TRUE)[, 1]
  
  y<- c(rep(1, length(pres)), rep(0, length(bg)))
  x<- c(pres, bg)
  as.numeric(pROC::auc(pROC::roc(y, x, quiet = TRUE)))
}

# Flip one layer's 1/2/3 scores (1<->3, 2 stays) — simulates reversing direction
flip_layer<- function(stk, layer_name) {
  s<- stk
  r<- s[[layer_name]]
  r_flip<- classify(r, cbind(c(1,2,3), c(3,2,1)))
  s[[layer_name]]<- r_flip
  s
}

cat("\nBaseline AUC (as submitted):\n")
auc_base<- auc_of(stack0)
cat(sprintf("  AUC = %.3f\n", auc_base))

cat("\nTesting each factor flipped one at a time:\n")
results<- data.frame(flipped = "none", AUC = round(auc_base, 3),
                      stringsAsFactors = FALSE)
for (nm in names(stack0)) {
  a<- auc_of(flip_layer(stack0, nm))
  results<- rbind(results,
                   data.frame(flipped = nm, AUC = round(a, 3)))
  cat(sprintf("  flip %-10s -> AUC = %.3f  %s\n", nm, a,
              if (a > auc_base + 0.03) "<-- IMPROVES" else ""))
}

cat("\nSummary (sorted):\n")
print(results[order(-results$AUC), ], row.names = FALSE)

best<- results[which.max(results$AUC), ]
cat(sprintf("\nBest single flip: %s (AUC %.3f -> %.3f)\n",
            best$flipped, auc_base, best$AUC))
if (best$flipped != "none" && best$AUC > auc_base + 0.03) {
  cat(sprintf(paste0(
    "\n>>> ACTION: the '%s' factor is scored in the wrong direction.\n",
    "    Fix it in R/03_reclassify.R (reverse its 1/2/3 mapping), then\n",
    "    re-run 03 -> 06. Re-validate; AUC should rise toward the value above.\n"),
    best$flipped))
} else {
  cat("\nNo single flip clearly improves AUC. The issue may be elsewhere\n")
  cat("(e.g. binning thresholds, or genuinely weak 5-factor signal). See notes.\n")
}

inv<- readRDS("output/inventory.rds")
sf::st_crs(inv)$epsg   #correct: 32646(UTM)

cat("\n--- end diagnostic ---\n\n")

# Best single flip: SoilType (AUC 0.646 -> 0.688)

#  ACTION: the 'SoilType' factor is scored in the wrong direction.
# Fix it in R/03_reclassify.R (reverse its 1/2/3 mapping), then
# re-run 03 -> 06. Re-validate; AUC should rise toward the value above.

# Note for the write-up:
# The 0.646 -> 0.688 finding is a nice hook for why I need the uncertainty
# analysis + confidence maps. If I just trust the expert weights, I'd overrate
# the hazard of some soil types. The Bayes + SMAA layer flags those areas as low
# decision confidence instead. For a disaster-management agency, knowing where the
# map is shaky is worth more than nudging the AUC up by 0.04.
