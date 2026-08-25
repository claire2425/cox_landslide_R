# R/diagnostic_test.R
# Diagnostic test (NOT a setting change)
# Purpose: To quantify the AUC penalty of maintaining theoretical/objective-driven settings compared to alternative empirical settings.
#        This is not to select a higher-scoring configuration, but to obtain empirical evidence for the dissertation.
#
# Run three configurations and record the AUC for each:
#   (A) BASELINE     : NDVI=ecological, SLOPE_FLAT=include <- Main analysis (default)
#   (B) NDVI_WO1     : NDVI=wo1,        SLOPE_FLAT=include <- Quantifies the cost of adhering to ecological theory
#   (C) SLOPE_EXCL   : NDVI=ecological, SLOPE_FLAT=exclude <- Quantifies the cost of retaining flat terrain
#
# Key principles:
#   * Alter only one variable at a time (A->B modifies only NDVI; A->C modifies only SLOPE).
#   * Clear the override using rm() after each run to ensure it does not contaminate the subsequent main analysis.
#   * The default values in 03_reclassify.R (ecological / include) remain unchanged throughout.

## Before running diagnostic_test.R, back up the main results in output/tables/, or have the diagnostic script write to different filenames to prevent overwriting the main files and causing confusion.

cat("\n========== DIAGNOSTIC TEST: quantify the cost of principled choices ==========\n")

suppressPackageStartupMessages({ library(terra) })
source("R/00_functions.R")

# Helper tool: Run the complete pipeline once (03->04->05), returning the AUC and validation table
run_once<- function(label) {
  cat(sprintf("\n----- [%s] running 03 -> 04 -> 05 -----\n", label))
  source("R/03_reclassify.R")
  source("R/04_hazard_overlay.R")
  source("R/05_validate.R")
  v<- readRDS("output/validation.rds")
  # save a separate copy of these outputs to prevent them from being overwritten by the next run
  file.copy("output/figures/roc_baseline.png",
            sprintf("output/figures/roc_%s.png", label), overwrite = TRUE)
  write.csv(v$validation,
            sprintf("output/tables/validation_%s.csv", label), row.names = FALSE)
  data.frame(setting = label, auc = round(v$auc, 4),
             lift_top2 = round(v$lift_top2, 2))
}

results<- list()

# ---- (A) BASELINE: Default values, no overrides set ----
# Ensure a clean environment: clear any residual overrides first
if (exists("NDVI_DIRECTION_OVERRIDE")) rm(NDVI_DIRECTION_OVERRIDE)
if (exists("SLOPE_FLAT_OVERRIDE"))     rm(SLOPE_FLAT_OVERRIDE)
results[["A_baseline"]] <- run_once("A_baseline")

# ---- (B) NDVI = wo1 (Modifies only NDVI, SLOPE remains 'include') ----
NDVI_DIRECTION_OVERRIDE <- "wo1"
results[["B_ndvi_wo1"]] <- run_once("B_ndvi_wo1")
rm(NDVI_DIRECTION_OVERRIDE)   # Clear immediately to revert to defaults

# ---- (C) SLOPE = exclude (Modifies only SLOPE, NDVI remains 'ecological') ----
SLOPE_FLAT_OVERRIDE <- "exclude"
results[["C_slope_excl"]] <- run_once("C_slope_excl")
rm(SLOPE_FLAT_OVERRIDE)       # Clear immediately to revert to defaults

# ---- Summary: Collate the AUCs of the three configurations side-by-side ----
summary_tbl <- do.call(rbind, results)
rownames(summary_tbl) <- NULL

cat("\n\n========== DIAGNOSTIC SUMMARY ==========\n")
print(summary_tbl)

# Calculate the 'cost': the change in AUC for the alternative settings relative to the baseline
base_auc<- summary_tbl$auc[summary_tbl$setting == "A_baseline"]
cat(sprintf("\nBaseline AUC = %.4f\n", base_auc))
for (i in seq_len(nrow(summary_tbl))) {
  s<- summary_tbl$setting[i]
  if (s != "A_baseline") {
    d<- summary_tbl$auc[i] - base_auc
    cat(sprintf("  %-14s AUC = %.4f  (%+.4f vs baseline)\n",
                s, summary_tbl$auc[i], d))
  }
}

write.csv(summary_tbl, "output/tables/diagnostic_summary.csv", row.names = FALSE)
cat("\nSaved: output/tables/diagnostic_summary.csv\n")
cat("Saved per-setting: roc_{A,B,C}.png + validation_{A,B,C}.csv\n")

# ---- Final safeguard: Confirm the default values for the main analysis have been restored ----
cat("\n--- safety check: overrides cleared? ---\n")
cat(sprintf("  NDVI_DIRECTION_OVERRIDE exists: %s\n", exists("NDVI_DIRECTION_OVERRIDE")))
cat(sprintf("  SLOPE_FLAT_OVERRIDE exists:     %s\n", exists("SLOPE_FLAT_OVERRIDE")))
cat("  (both should be FALSE -> your formal analysis defaults are intact)\n")
cat("\nNOTE: 03_reclassify.R defaults remain ecological / include. This script only\n")
cat("      measured the cost of alternatives; it changed nothing permanent.\n")
cat("=========================================\n\n")

