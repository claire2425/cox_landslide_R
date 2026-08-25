# R/03b_slope_flat_comparison.R
# OPTION C SUPPORT: compare the two ways of treating gentle slopes.
#   include  = flat land kept as lowest susceptibility (MAIN analysis, B)
#   exclude  = flat land set to Restricted/NA per WO1   (APPENDIX, A)

# Produces a small table of coverage + AUC for each, so you can justify the
# main choice and report the comparison in the appendix.

# Run AFTER the pipeline works end-to-end. It re-runs 03->04->05 internally
# for each slope mode by sourcing them with SLOPE_FLAT overridden.

cat("\n--- 03b: slope-flat include vs exclude comparison ---\n")
suppressPackageStartupMessages({ library(terra); library(sf) })
if (!exists("ahp_weights")) source("R/00_functions.R")

weights<- readRDS("output/ahp_weights.rds")
inv<- readRDS("output/inventory.rds")
gt<- as.character(unique(st_geometry_type(inv)))[1]
inv_v<- vect(if (gt %in% c("POLYGON","MULTIPOLYGON")) st_centroid(inv) else inv)

auc_of<- function(haz) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NA)
  set.seed(2026)
  pres<- extract(haz, inv_v, ID = FALSE)[, 1]; pres<- pres[!is.na(pres)]
  bg<- spatSample(haz, max(500, length(pres)*10), "random", na.rm = TRUE)[, 1]
  as.numeric(pROC::auc(pROC::roc(c(rep(1,length(pres)),rep(0,length(bg))),
                                 c(pres,bg), quiet = TRUE)))
}

run_mode<- function(mode) {
  
  # include-mode should match the MAIN pipeline exactly, so reuse the hazard
  # that 04_hazard_overlay already produced and 05_validate scores (the
  # name-aligned, sum-to-one-checked hazard). Only exclude-mode needs a fresh
  # re-run, since that scenario is not what the main pipeline saved.
  
  if (mode == "include" && file.exists("output/hazard_continuous.rds")) {
    haz<- unwrap(readRDS("output/hazard_continuous.rds"))
  } else {
    SLOPE_FLAT_OVERRIDE<- mode
    
    # with local = TRUE, script 03 can see the SLOPE_FLAT_OVERRIDE set just above
    source("R/03_reclassify.R", local = TRUE)
    stk<- unwrap(readRDS("output/reclass_stack.rds"))
    w<- weights[names(stk)]
    haz<- app(stk, fun = function(v) sum(v * w))
  }
  total<- ncell(haz)
  na_pct<- 100 * sum(is.na(values(haz))) / total
  list(mode = mode, na_pct = round(na_pct, 1), auc = round(auc_of(haz), 3))
}

cat("\nNOTE: For this to switch modes, add ONE line near the top of\n")
cat("R/03_reclassify.R:\n")
cat('  if (exists("SLOPE_FLAT_OVERRIDE")) SLOPE_FLAT<- SLOPE_FLAT_OVERRIDE\n\n')

res<- rbind(
  as.data.frame(run_mode("include")),
  as.data.frame(run_mode("exclude"))
)
res$interpretation<- c("MAIN (susceptibility, full coverage)",
                        "APPENDIX (hazard zonation, flat excluded)")
cat("\n=== Slope-flat comparison ===\n")
print(res, row.names = FALSE)
write.csv(res, "output/tables/slope_flat_comparison.csv", row.names = FALSE)
cat("\nSaved: output/tables/slope_flat_comparison.csv\n")

# added here: automatically clean up and restore the main-analysis state
cat("\n--- Auto-restoring the MAIN analysis state on disk ---\n")
# just to be safe, make sure no override variable is left in the global environment
if (exists("SLOPE_FLAT_OVERRIDE", envir = .GlobalEnv)) {
  rm(SLOPE_FLAT_OVERRIDE, envir = .GlobalEnv)
}
# re-run 03 once more; with no override it falls back to the default "include",
# overwriting reclass_stack.rds on disk with the correct MAIN version.
source("R/03_reclassify.R") 
cat("\nDone! Environment and output files safely restored to MAIN configuration.\n\n")

