# R/03_reclassify.R
# Step 3: Reclassify each layer to a 1/2/3 suitability score (NA = Restricted),
# following the WO1.txt rules.

cat("\n--- 03_reclassify: apply WO1.txt reclassification rules ---\n")

suppressPackageStartupMessages({ library(terra) })
if (!exists("reclass_codes")) source("R/00_functions.R")

if (!file.exists("output/rasters_raw.rds")) stop("Run R/02_load_data.R first.")
rasters<- lapply(readRDS("output/rasters_raw.rds"), unwrap)

# Where pre-binned 1-5 rasters live (from the GIS practical's Reclass step).
reclass_dir<- "data/LSM Tutorial/Reclass"
slope_grid<- "data/LSM Tutorial/reclass_slope"

# Helper that prints the codes present, warns about unmapped ones.
report_codes<- function(r, layer_name, map) {
  u<- sort(unique(values(r, na.rm = TRUE)))
  cat(sprintf("  %-10s codes present: {%s}\n", layer_name, paste(u, collapse = ", ")))
  unmapped<- setdiff(u, as.numeric(names(map)))
  if (length(unmapped) > 0) {
    cat(sprintf("    NOTE: codes {%s} are not in your mapping -> will become NA.\n",
                paste(unmapped, collapse = ", ")))
    cat("    If that is wrong, fix the mapping after checking the .vat.dbf file.\n")
  }
}

# Continuous factors: use pre-binned 1-5 rasters if present, else self-bin

# ---- Slope (direction: steeper = more hazardous) ----------------------------
if (dir.exists(slope_grid)) {
  cat("Slope: using pre-binned reclass_slope/ grid\n")
  slope_5<- rast(slope_grid)
} else {
  cat("Slope: pre-binned grid not found -> binning raw slope (degrees)\n")
  slope_5<- classify(rasters$Slope, matrix(c(
    -Inf, 5, 1,   5, 10, 2,   10, 15, 3,   15, 25, 4,   25, Inf, 5
  ), ncol = 3, byrow = TRUE))
}

# SLOPE_FLAT decides how I treat the gentlest slope bin (bin1, 0-2.8 degrees). Cox's Bazar has a lot of
# flat land (~47% of the study area), so this one choice changes the coverage a lot:
#   "include" -> bin1 becomes 1 (lowest susceptibility), so the whole area gets a relative class. This is the
#                main analysis -- for a susceptibility map I want full spatial coverage.
#   "exclude" -> bin1 becomes NA (Restricted), following the original WO1 rule. That's strict hazard zonation
#                with flat land dropped; keep it for the appendix comparison.

SLOPE_FLAT<- "include"
if (exists("SLOPE_FLAT_OVERRIDE")) SLOPE_FLAT<- SLOPE_FLAT_OVERRIDE  # for 07b comparison
slope_map<- if (SLOPE_FLAT == "exclude") {
  c("1"=NA, "2"=1, "3"=2, "4"=3, "5"=3)    # WO1 original: flat excluded
} else {
  c("1"=1,  "2"=1, "3"=2, "4"=3, "5"=3)    # flat = lowest susceptibility
}
slope_reclass<- reclass_codes(slope_5, slope_map)

# ---- NDVI (direction: LOWER NDVI = less vegetation = more hazardous) --------
ndvi_tif<- file.path(reclass_dir, "Reclass_NDVI.tif")
if (file.exists(ndvi_tif)) {
  cat("NDVI:  using pre-binned Reclass_NDVI.tif\n")
  ndvi_5<- rast(ndvi_tif)
} else {
  cat("NDVI:  pre-binned tif not found -> binning raw NDVI\n")
  # NOTE: bin so that bin 1 = LOWEST NDVI. Verify with the data range.
  ndvi_5<- classify(rasters$NDVI, matrix(c(
    -Inf, 0.0, 1,   0.0, 0.1, 2,   0.1, 0.2, 3,   0.2, 0.35, 4,   0.35, Inf, 5
  ), ncol = 3, byrow = TRUE))
}

# WO1 as published scores high NDVI as high hazard (bin5 -> 3). But denser vegetation usually stabilises
# slopes (so it should be lower hazard), and the baseline validation (AUC 0.517, a U-shaped lift) tells me
# this factor is reversed. So NDVI_DIRECTION lets me switch:
#   "wo1"        -> keep WO1.txt as-is (high NDVI -> hazard 3)
#   "ecological" -> the corrected version I use: low NDVI (sparse veg) -> high hazard

NDVI_DIRECTION<- "ecological"
if (exists("NDVI_DIRECTION_OVERRIDE")) NDVI_DIRECTION<- NDVI_DIRECTION_OVERRIDE  # for diagnostic test
ndvi_map<- if (NDVI_DIRECTION == "wo1") {
  c("1"=NA, "2"=1, "3"=2, "4"=2, "5"=3)         # original WO1
} else {
  c("1"=NA, "2"=3, "3"=2, "4"=2, "5"=1)         # corrected: low NDVI = hazardous
}
ndvi_reclass<- reclass_codes(ndvi_5, ndvi_map)

# ---- Distance to road (direction: CLOSER = more hazardous) ------------------
road_tif<- file.path(reclass_dir, "Reclass_road.tif")
if (file.exists(road_tif)) {
  cat("Road:  using pre-binned Reclass_road.tif\n")
  road_5<- rast(road_tif)
} else {
  cat("Road:  pre-binned tif not found -> binning raw distance (m)\n")
  road_5<- classify(rasters$DistRoad, matrix(c(
    -Inf, 100, 1,  100, 500, 2,  500, 1000, 3,  1000, 2000, 4,  2000, Inf, 5
  ), ncol = 3, byrow = TRUE))
}
# WO1: bin 1->3, 2->3, 3->2, 4->1, 5->NA  (inverse: near road scores high)
road_reclass<- reclass_codes(road_5,
                              c("1"=3, "2"=3, "3"=2, "4"=1, "5"=NA))

# Categorical factors: map by raster code. VERIFY codes against .vat.dbf!

# === [VAT Code-to-Class frequency check] ===
cat("\n--- Verify: Land Cover and Soil Type original pixel distribution ---\n")
cat("Land Cover frequency:\n")
print(terra::freq(rasters$LandCover))

cat("\nSoil Type frequency:\n")
print(terra::freq(rasters$SoilType))

# ---- Land Cover -------------------------------------------------------------
# VERIFIED against Land Cover.tif.vat.dbf by R/07_verify_inputs.R (CHECK 1: PASS).
# Codes are 0-based:
#   0 = Bare Land      -> 1
#   1 = Built Up Area  -> 3  (most hazardous land cover)
#   2 = Crop Land      -> 2
#   3 = Vegetation     -> 3
#   4 = Water          -> NA (Restricted; excluded from analysis)
landcover_map<- c("0"=1, "1"=3, "2"=2, "3"=3, "4"=NA)
report_codes(rasters$LandCover, "LandCover", landcover_map)
landcover_reclass<- reclass_codes(rasters$LandCover, landcover_map)

# ---- Soil Type --------------------------------------------------------------
# VERIFIED against Soil Type.tif.vat.dbf by R/07_verify_inputs.R (CHECK 1: PASS).
# Codes are 0-based:
#   0 = Noncalcareous aluvium         -> 2
#   1 = Noncalcareous grey floodplain -> 2
#   2 = Noncalcareous brown floodplain-> NA (Restricted; only 150 px, negligible)
#   3 = Acid sulphate                 -> 1
#   4 = Grey piedmont                 -> 2
#   5 = Brown hill                    -> 3  (most hazardous soil)
soil_map<- c("0"=2, "1"=2, "2"=NA, "3"=1, "4"=2, "5"=3)
report_codes(rasters$SoilType, "SoilType", soil_map)
soil_reclass<- reclass_codes(rasters$SoilType, soil_map)

# Assemble stack (aligned, canonical order)
ref<- slope_reclass
ndvi_reclass<- align_to(ndvi_reclass,      ref, categorical = TRUE)
road_reclass<- align_to(road_reclass,      ref, categorical = TRUE)
landcover_reclass<- align_to(landcover_reclass, ref, categorical = TRUE)
soil_reclass<- align_to(soil_reclass,      ref, categorical = TRUE)

reclass_stack<- c(slope_reclass, landcover_reclass, ndvi_reclass,
                   soil_reclass, road_reclass)
names(reclass_stack)<- c("Slope", "LandCover", "NDVI", "SoilType", "DistRoad")

cat("\nReclassified layers (should each be a subset of {1,2,3}):\n")
total_cells<- ncell(reclass_stack[[1]])
for (nm in names(reclass_stack)) {
  v<- values(reclass_stack[[nm]], na.rm = TRUE)
  na<- sum(is.na(values(reclass_stack[[nm]])))
  bad<- setdiff(unique(v), c(1, 2, 3))
  cat(sprintf("  %-10s values={%s}  NA=%d (%.1f%%)%s\n",
              nm, paste(sort(unique(v)), collapse = ", "),
              na, 100 * na / total_cells,
              if (length(bad) > 0) sprintf("  <-- UNEXPECTED: %s",
                                           paste(bad, collapse=",")) else ""))
}

# NA coverage. The rasters are rectangular but Cox's Bazar is irregular, so a
# large share of NA is EXPECTED (area outside the study boundary), not a fault.
# Only a genuine flag if ONE layer is far more NA than the others.
na_any<- sum(is.na(values(sum(reclass_stack))))
per_layer_na<- sapply(names(reclass_stack), function(nm) mean(is.na(values(reclass_stack[[nm]]))))
cat(sprintf("\nPixels NA in at least one layer (-> Restricted): %.1f%%\n", 100 * na_any / total_cells))
cat(sprintf("Per-layer NA %.1f%%-%.1f%% (similar => out-of-boundary frame, expected).\n",
            100*min(per_layer_na), 100*max(per_layer_na)))
if (max(per_layer_na) - min(per_layer_na) > 0.25)
  cat("  CHECK: one layer far more NA than others -> inspect its reclass rule.\n")

saveRDS(wrap(reclass_stack), "output/reclass_stack.rds")
cat("\nSaved: output/reclass_stack.rds\n\n")

