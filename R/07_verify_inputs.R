# R/07_verify_inputs.R
# Step 7: VERIFY the two inputs the README flags as "you MUST confirm".
#
# Until now these were left as manual TODOs in the comments of 03_reclassify.R:
#   (1) Land Cover / Soil Type numeric code -> class mapping
#   (2) Reclass DIRECTION (are landslide points concentrated in High zones?)
#
# This script performs BOTH checks programmatically, against the actual data
# files, and writes a pass/fail report to output/tables/input_verification.csv.
# It is read-only w.r.t. the analysis: it changes no weights and no rasters,
# it only audits them. Run it AFTER 01->05 (it needs reclass_stack + inventory),
# or standalone for the mapping check alone (it only needs the .vat.dbf files).
#
# Design: every check prints PASS / FAIL / WARN and the evidence behind it, so
# the verification itself is reproducible and examinable — not a claim in a
# comment, but a re-runnable test.

cat("\n--- 07_verify_inputs: audit code->class mapping and reclass direction ---\n")

# Note: terra/sf are only needed for CHECK 2 (the AUC direction test); they are
# loaded lazily there. CHECK 1 (the .vat.dbf mapping audit) is pure base R, so
# it runs even on a machine without the geospatial stack installed.
if (!exists("reclass_codes")) source("R/00_functions.R")

verification<- list()   # collect rows for the report
add_row<- function(check, status, detail) {
  verification[[length(verification) + 1]] <<-
    data.frame(check = check, status = status, detail = detail,
               stringsAsFactors = FALSE)
}

# CHECK 1 — Land Cover / Soil Type code -> class mapping vs the .vat.dbf
# The mapping used in 03_reclassify.R assumes which integer raster code means
# which class. The authoritative source is each layer's ESRI VAT (.vat.dbf),
# which stores Value -> Class_name/Type. We read it directly and compare.
#
# We rebuild the EXPECTED suitability map from two authoritative sources jointly:
#   - .vat.dbf      : Value(code) -> class name
#   - WO1.txt rules : class name  -> suitability score (1/2/3, Restricted=NA)
# and check that 03_reclassify.R's hard-coded *_map matches that join exactly.
cat("\n========== CHECK 1: code->class mapping vs .vat.dbf ==========\n")

read_vat<- function(path) {
  # Minimal dBASE III reader (no external pkg): enough for ESRI VAT files.
  con<- file(path, "rb"); on.exit(close(con))
  hdr<- readBin(con, "raw", 32)
  n_rec<- readBin(hdr[5:8],   "integer", size = 4, endian = "little")
  hdr_len<- readBin(hdr[9:10],  "integer", size = 2, endian = "little")
  rec_len<- readBin(hdr[11:12], "integer", size = 2, endian = "little")
  n_fields<- (hdr_len - 33L) %/% 32L
  fields<- vector("list", n_fields)
  for (i in seq_len(n_fields)) {
    fd<- readBin(con, "raw", 32)
    nm<- rawToChar(fd[which(fd[1:11] != as.raw(0))[seq_len(sum(fd[1:11] != as.raw(0)))]])
    fields[[i]]<- list(
      name = nm,
      type = rawToChar(fd[12]),
      len  = as.integer(fd[17]))
  }
  readBin(con, "raw", 1)  # header terminator 0x0D
  recs<- vector("list", n_rec)
  for (r in seq_len(n_rec)) {
    raw_rec<- readBin(con, "raw", rec_len)
    pos<- 2L  # skip deletion flag
    row<- list()
    for (f in fields) {
      val<- trimws(rawToChar(raw_rec[pos:(pos + f$len - 1L)]))
      row[[f$name]]<- val
      pos<- pos + f$len
    }
    recs[[r]]<- row
  }
  do.call(rbind, lapply(recs, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
}

# WO1.txt suitability rules, keyed by CLASS NAME (Restricted -> NA)
wo1_landcover<- c("Bare Land" = 1, "Built Up Area" = 3, "Crop Land" = 2,
                   "Vegetation" = 3, "Water" = NA)
wo1_soil<- c("Noncalcareous aluvium" = 2, "Noncalcareous grey floodplain" = 2,
              "Noncalcareous brown floodplain" = NA, "Acid sulphate" = 1,
              "Grey piedmont" = 2, "Brown hill" = 3)

# The mappings 03_reclassify.R actually uses (keep in sync with that script):
landcover_map_used<- c("0"=1, "1"=3, "2"=2, "3"=3, "4"=NA)
soil_map_used<- c("0"=2, "1"=2, "2"=NA, "3"=1, "4"=2, "5"=3)

verify_categorical<- function(vat_path, name_field, wo1_by_name, map_used,
                               label) {
  if (!file.exists(vat_path)) {
    add_row(sprintf("Mapping: %s", label), "SKIP",
            sprintf(".vat.dbf not found at %s", vat_path))
    cat(sprintf("  %-10s SKIP (no .vat.dbf)\n", label)); return(invisible())
  }
  vat<- read_vat(vat_path)
  vat$Value<- as.integer(vat$Value)
  cls<- vat[[name_field]]
  # expected score per code = WO1 score of that code's class name
  expected<- wo1_by_name[cls]
  names(expected)<- as.character(vat$Value)
  expected<- expected[order(as.integer(names(expected)))]
  used<- map_used[names(expected)]

  cat(sprintf("\n  %s — joining .vat.dbf (code->class) with WO1 (class->score):\n", label))
  cmp<- data.frame(code = names(expected),
                    class = cls[order(vat$Value)],
                    expected_score = unname(expected),
                    script_uses    = unname(used),
                    stringsAsFactors = FALSE)
  print(cmp, row.names = FALSE)

  same<- mapply(function(a, b) (is.na(a) && is.na(b)) || (!is.na(a) && !is.na(b) && a == b),
                 expected, used)
  if (all(same)) {
    add_row(sprintf("Mapping: %s", label), "PASS",
            "03_reclassify *_map matches .vat.dbf joined with WO1 rules")
    cat(sprintf("  -> %s: PASS (mapping verified against .vat.dbf)\n", label))
  } else {
    bad<- names(expected)[!same]
    add_row(sprintf("Mapping: %s", label), "FAIL",
            sprintf("codes {%s} disagree with .vat.dbf+WO1", paste(bad, collapse=",")))
    cat(sprintf("  -> %s: FAIL at codes {%s} — fix *_map in 03_reclassify.R\n",
                label, paste(bad, collapse=",")))
  }
}

data_dir<- "data/LSM Tutorial"
verify_categorical(file.path(data_dir, "Land Cover.tif.vat.dbf"),
                   "Class_name", wo1_landcover, landcover_map_used, "LandCover")
verify_categorical(file.path(data_dir, "Soil Type.tif.vat.dbf"),
                   "Type", wo1_soil, soil_map_used, "SoilType")

# CHECK 2 — Reclass DIRECTION via the inventory overlay (AUC test)
# A correct susceptibility map puts known landslides in High/Very-High zones.
# Quantitatively: AUC of (inventory presence vs background) should be > 0.5;
# AUC < 0.5 means at least one factor is scored backwards (direction reversed).

cat("\n========== CHECK 2: reclass direction (AUC) ==========\n")

have_stack<- file.exists("output/reclass_stack.rds")
have_w<- file.exists("output/ahp_weights.rds")
have_inv<- file.exists("output/inventory.rds")

if (!(have_stack && have_w && have_inv)) {
  add_row("Direction: AUC > 0.5", "SKIP",
          "needs reclass_stack.rds + ahp_weights.rds + inventory.rds (run 01->04 first)")
  cat("  SKIP: run R/01..R/04 first so the stack + inventory exist.\n")
} else if (!requireNamespace("terra", quietly = TRUE) ||
           !requireNamespace("pROC", quietly = TRUE) ||
           !requireNamespace("sf", quietly = TRUE)) {
  add_row("Direction: AUC > 0.5", "SKIP", "terra/sf/pROC not all installed")
  cat("  SKIP: install terra, sf, pROC to run the AUC direction test.\n")
} else {
  suppressPackageStartupMessages({ library(terra); library(sf) })
  stk<- unwrap(readRDS("output/reclass_stack.rds"))
  w<- readRDS("output/ahp_weights.rds"); w<- w[names(stk)]
  inv<- readRDS("output/inventory.rds")
  gt<- as.character(unique(st_geometry_type(inv)))[1]
  pts<- if (gt %in% c("POLYGON","MULTIPOLYGON")) st_centroid(inv) else inv
  pv<- vect(pts)

  haz<- terra::app(stk, fun = function(v) sum(v * w))
  set.seed(2026)
  pres<- terra::extract(haz, pv, ID = FALSE)[, 1]; pres<- pres[!is.na(pres)]
  bg<- terra::spatSample(haz, max(500, length(pres) * 10),
                            method = "random", na.rm = TRUE)[, 1]
  y<- c(rep(1, length(pres)), rep(0, length(bg)))
  x<- c(pres, bg)
  auc<- as.numeric(pROC::auc(pROC::roc(y, x, quiet = TRUE)))

  cat(sprintf("\n  Inventory points scored: %d (NA-dropped)\n", length(pres)))
  cat(sprintf("  Background sample: %d\n", length(bg)))
  cat(sprintf("  AUC = %.3f\n", auc))

  if (auc >= 0.55) {
    add_row("Direction: AUC > 0.5", "PASS",
            sprintf("AUC=%.3f (>=0.55): landslides concentrate in high-hazard zones", auc))
    cat("  -> PASS: direction is correct (landslides fall in High zones).\n")
  } else if (auc >= 0.50) {
    add_row("Direction: AUC > 0.5", "WARN",
            sprintf("AUC=%.3f (0.50-0.55): weak but not reversed; inspect lift/binning", auc))
    cat("  -> WARN: weak signal; not reversed, but check binning/thresholds.\n")
  } else {
    add_row("Direction: AUC > 0.5", "FAIL",
            sprintf("AUC=%.3f (<0.50): a factor is reversed — run R/05b to find it", auc))
    cat("  -> FAIL: AUC<0.5 => a factor reversed. Run R/05b_direction_diagnostic.R.\n")
  }
}

# Write report
report<- do.call(rbind, verification)
cat("\n========== VERIFICATION REPORT ==========\n")
print(report, row.names = FALSE)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(report, "output/tables/input_verification.csv", row.names = FALSE)
cat("\nSaved: output/tables/input_verification.csv\n")

n_fail<- sum(report$status == "FAIL")
n_warn<- sum(report$status == "WARN")
n_skip<- sum(report$status == "SKIP")
cat(sprintf("\nSummary: %d PASS, %d WARN, %d FAIL, %d SKIP\n",
            sum(report$status == "PASS"), n_warn, n_fail, n_skip))
if (n_fail > 0)
  cat(">> Action needed: resolve FAIL rows before trusting the hazard map.\n")
if (n_fail == 0 && n_warn == 0 && n_skip == 0)
  cat(">> All inputs verified. The two README 'MUST confirm' items are closed.\n")

cat("\n--- end 07_verify_inputs ---\n\n")

