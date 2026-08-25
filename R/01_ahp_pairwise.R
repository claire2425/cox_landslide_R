# 01_ahp_pairwise.R
# Step 1: set up the baseline criterion weights, plus a comparison pairwise matrix for the CR discussion.
#
# A note on where the weights come from, because it matters for the whole project:
# the dataset (Alam & Ahmed's LSM tutorial for Cox's Bazar) gives the factor weights directly
# as integers in WO1.txt -- Slope 25 / Distance-to-road 15 / NDVI 20 / Land Cover 20 / Soil Type 20,
# summing to 100. So strictly this is a Weighted Linear Combination, not an AHP pairwise elicitation:
# there's no comparison matrix and hence no consistency ratio attached. That's actually the whole
# motivation here -- the weights are handed to us as fixed, certain numbers with no sense of how
# uncertain the expert judgement behind them is. What I do is take this exact weight vector and push
# uncertainty through it (OAT -> Monte Carlo -> Bayesian) to get a per-pixel confidence map.
#
# WEIGHT_SOURCE:
#   "wo1"   -> the actual WO1.txt weights. This is the honest baseline; no CR (no matrix behind it).
#   "ahmed" -> a 5x5 pairwise matrix pulled straight out of Ahmed (2015) Table 3's full 9x9, no manual
#              tweaking. Comparison scenario only -- it's a genuine matrix with a real CR = 0.0824 (< 0.10).
#              The full reproducible build (rebuild the 9x9 -> check against the paper -> extract the
#              5-factor submatrix) lives in 01b_build_ahp_matrix.R.

cat("\n--- 01_ahp_pairwise: establish baseline weights ---\n")

if (!exists("ahp_weights")) source("R/00_functions.R")

criteria<- c("Slope", "LandCover", "NDVI", "SoilType", "DistRoad")

# ---- choose the baseline ----
WEIGHT_SOURCE<- "wo1"   # "wo1" = honest baseline (recommended); "ahmed" = comparison only
if (exists("WEIGHT_SOURCE_OVERRIDE")) WEIGHT_SOURCE<- WEIGHT_SOURCE_OVERRIDE

# (1) WO1.txt weights -- the honest baseline, straight from the supplied data
w_wo1<- c(Slope = 0.25, LandCover = 0.20, NDVI = 0.20,
          SoilType = 0.20, DistRoad = 0.15)[criteria]
stopifnot(abs(sum(w_wo1) - 1) < 1e-9)   # sanity: the weights must sum to 1

# (2) Ahmed (2015) comparison matrix -- the 5 relevant factors' pairwise values taken directly from
#     the published 9x9 (verified identical to 01b's from-scratch rebuild). No manual adjustment.
A_ahmed<- matrix(c(
  1,    2,    6,    1/5,  1/5,
  1/2,  1,    5,    1/4,  1/6,
  1/6,  1/5,  1,    1/6,  1/9,
  5,    4,    6,    1,    1/2,
  5,    6,    9,    2,    1
), nrow = 5, byrow = TRUE,
   dimnames = list(c("LandCover","NDVI","DistRoad","Slope","SoilType"),
                   c("LandCover","NDVI","DistRoad","Slope","SoilType")))
A_ahmed<- A_ahmed[criteria, criteria]   # reorder rows/cols to match the criteria vector
w_ahmed<- ahp_weights(A_ahmed)
cr_ahmed<- ahp_cr(A_ahmed)

# (3) Ahmed's 9-factor subset, just renormalised -- kept only as a transparency reference
ahmed_full<- c(Slope = 0.1989, LandCover = 0.0975, NDVI = 0.0706,
               SoilType = 0.3074, DistRoad = 0.0243)
w_ahmed_subset<- (ahmed_full / sum(ahmed_full))[criteria]

# ---- pick the active weight vector, and matrix if there is one ----
if (WEIGHT_SOURCE == "wo1") {
  w<- w_wo1; A<- NULL; cr<- NULL
  cat("\nUsing WEIGHT_SOURCE = 'wo1' (honest baseline, directly from WO1.txt).\n")
  cat("Pre-assigned WLC weights; no pairwise matrix and no CR.\n")
} else if (WEIGHT_SOURCE == "ahmed") {
  w<- w_ahmed; A<- A_ahmed; cr<- cr_ahmed
  cat("\nUsing WEIGHT_SOURCE = 'ahmed' (comparison scenario, Ahmed Table 3 submatrix).\n")
} else stop("WEIGHT_SOURCE must be 'wo1' or 'ahmed'.")

cat("\nActive baseline weights:\n"); print(round(w[criteria], 4))
if (!is.null(A)) {
  check_ahp_matrix(A)
  cat(sprintf("\nComparison matrix CR = %.4f %s\n", cr$CR,
              if (cr$CR < 0.10) "[< 0.10 acceptable]" else "[>= 0.10 review]"))
}

# put the three weight versions side by side, for the methodology/appendix
provenance<- data.frame(
  factor               = criteria,
  WO1_baseline         = round(as.numeric(w_wo1[criteria]), 4),
  Ahmed_submatrix      = round(as.numeric(w_ahmed[criteria]), 4),
  Ahmed_subset_renorm  = round(as.numeric(w_ahmed_subset[criteria]), 4)
)
cat("\nProvenance comparison (for the methodology/appendix):\n"); print(provenance)
cat("\nNote: WO1 is the baseline (data-faithful). Ahmed columns are context only,\n")
cat("NOT used in the main analysis.\n")

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(provenance, "output/tables/ahp_weights.csv", row.names = FALSE)
saveRDS(w,  "output/ahp_weights.rds")
saveRDS(A,  "output/ahp_matrix.rds")
saveRDS(cr, "output/ahp_diagnostics.rds")
cat("\nSaved: ahp_weights.csv/.rds, ahp_matrix.rds, ahp_diagnostics.rds\n\n")

