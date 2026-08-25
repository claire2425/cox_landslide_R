# R/13_smaa.R
# Step 13: SMAA - Formalising the confidence map (RQ3c).
#
# What it does: Takes the posterior draws and calculates the SMAA stats.
# The most useful part here is the central weights: "What kind of weighting 
# actually pushes a pixel into the 'High' class?"
#
# Sanity check: The rank-1 acceptability must perfectly match the P(modal) from Step 12.

cat("\n--- 13_smaa: rank acceptability + central weights per pixel ---\n")
suppressMessages({ library(terra) })
set.seed(2026)

CENTRAL_WEIGHTS_PER_PIXEL<- TRUE   # FALSE -> map-wide central weights only (low memory)

## 1. Load posterior weights and the reclassified factor stack (as in 12)
W<- readRDS("output/bayes_weights.rds")            # M x 5
factors<- colnames(W); M<- nrow(W); nf<- ncol(W)
cat(sprintf("Loaded %d posterior weight draws over %d factors.\n", M, nf))

stack_path<- "output/reclass_stack.rds"
if (file.exists(stack_path)) {
  rs<- readRDS(stack_path); if (!inherits(rs, "SpatRaster")) rs<- rast(rs)
} else {
  rs<- rast(list(
    Slope     = rast("output/reclass_Slope.tif"),
    LandCover = rast("output/reclass_LandCover.tif"),
    NDVI      = rast("output/reclass_NDVI.tif"),
    SoilType  = rast("output/reclass_SoilType.tif"),
    DistRoad  = rast("output/reclass_DistRoad.tif")))
}
rs<- rs[[ match(factors, names(rs)) ]]
stopifnot(all(names(rs) == factors))

## 2. SAME class breaks as 12 (do not re-derive a different rule)
get_breaks<- function() {
  if (file.exists("output/hazard_breaks.rds")) {
    b<- readRDS("output/hazard_breaks.rds")
    if (is.numeric(b) && length(b) == 6) return(b)
  }
  if (file.exists("output/hazard_class.rds")) {
    hc<- readRDS("output/hazard_class.rds")
    if (is.list(hc) && !is.null(hc$breaks) && length(hc$breaks) == 6) return(hc$breaks)
  }
  if (file.exists("output/hazard_continuous.rds")) {
    hcont<- readRDS("output/hazard_continuous.rds")
    sc<- values(if (inherits(hcont, "SpatRaster")) hcont else rast(hcont))
    sc<- sc[is.finite(sc)]
    return(as.numeric(quantile(sc, probs = seq(0, 1, by = 0.2))))
  }
  stop("Could not obtain 6 class breaks -- must match 12_hazard_posterior.R.")
}
brk<- get_breaks(); stopifnot(is.numeric(brk), length(brk) == 6)
cat("Using fixed baseline class breaks (same as 12):\n"); print(round(brk, 4))

vals<- values(rs); ok<- stats::complete.cases(vals)
Xok<- vals[ok, , drop = FALSE]; n_ok<- nrow(Xok)
cat(sprintf("%d valid pixels (of %d).\n", n_ok, nrow(vals)))

cut_to_class<- function(score) {
  cl<- findInterval(score, brk, rightmost.closed = TRUE)
  cl[cl < 1]<- 1L; cl[cl > 5]<- 5L; as.integer(cl)
}

## 3. Stream the posterior: class counts + per-class weight sums
class_count<- matrix(0, nrow = n_ok, ncol = 5)      # double, not 0L: 8000 draws x ~1M pixels overflows int -> acceptability indices
# per-class running sum of weight vectors (for central weights)
if (CENTRAL_WEIGHTS_PER_PIXEL) {
  wsum_by_class<- vector("list", 5)
  for (c in 1:5) wsum_by_class[[c]]<- matrix(0, nrow = n_ok, ncol = nf)
}
# map-wide weight sums per class (always cheap)
wsum_mapwide<- matrix(0, nrow = 5, ncol = nf, dimnames = list(paste0("class",1:5), factors))
count_mapwide<- numeric(5)   # double, not integer(5): sum over 8000 draws x ~1M pixels exceeds .Machine$integer.max

cat("Streaming draws"); pb<- txtProgressBar(min = 0, max = M, style = 3)
for (m in 1:M) {
  wv<- W[m, ]
  score<- as.numeric(Xok %*% wv)
  cl<- cut_to_class(score)
  # acceptability counts
  class_count[cbind(seq_len(n_ok), cl)]<- class_count[cbind(seq_len(n_ok), cl)] + 1
  # central-weight accumulation
  if (CENTRAL_WEIGHTS_PER_PIXEL) {
    for (c in 1:5) {
      sel<- which(cl == c)
      if (length(sel)) wsum_by_class[[c]][sel, ]<- sweep(wsum_by_class[[c]][sel, , drop=FALSE], 2, wv, "+")
    }
  }
  # map-wide central weights: how many pixels in each class this draw, weight wv
  tb<- tabulate(cl, nbins = 5)
  wsum_mapwide<- wsum_mapwide + outer(tb, wv)
  count_mapwide<- count_mapwide + tb
  setTxtProgressBar(pb, m)
}
close(pb)

## 4. SMAA quantities
accept<- class_count / M                       # n_ok x 5  rank acceptability
modal_class<- max.col(accept, ties.method = "first")
confidence<- accept[cbind(seq_len(n_ok), modal_class)]   # confidence factor
p_high<- accept[, 4] + accept[, 5]             # P(class >= High)

# CONSISTENCY CHECK vs 12: rank-1 (modal) acceptability == 12's P(modal)
if (file.exists("output/conf_pmodal.rds")) {
  pm12<- values(readRDS("output/conf_pmodal.rds"))[ok]
  d<- max(abs(pm12 - confidence), na.rm = TRUE)
  cat(sprintf("Consistency check vs 12 P(modal): max abs diff = %.4g %s\n",
              d, ifelse(d < 1e-6, "(identical, as expected)",
                        "(NONZERO -- breaks/overlay/mask drifted; investigate)")))
}

# per-pixel central weight vector for the MODAL class
if (CENTRAL_WEIGHTS_PER_PIXEL) {
  central_modal<- matrix(NA_real_, nrow = n_ok, ncol = nf, dimnames = list(NULL, factors))
  for (c in 1:5) {
    sel<- which(modal_class == c)
    if (length(sel)) {
      cnt<- class_count[cbind(sel, rep(c, length(sel)))]
      central_modal[sel, ]<- wsum_by_class[[c]][sel, ] / pmax(cnt, 1)
    }
  }
}

# map-wide central weight per class (always)
central_mapwide<- wsum_mapwide / pmax(count_mapwide, 1)

## 5. Write rasters + tables
mk_raster<- function(vec) { r<- rs[[1]]; r[]<- NA_real_; r[ok]<- vec; r }
saveRDS(accept,                       "output/smaa_acceptability.rds")
saveRDS(mk_raster(confidence),        "output/smaa_confidence.rds")
if (CENTRAL_WEIGHTS_PER_PIXEL) saveRDS(central_modal, "output/smaa_central_weights.rds")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
cw_tab<- data.frame(class = paste0("class", 1:5),
                     label = c("Very Low","Low","Moderate","High","Very High"),
                     round(as.data.frame(central_mapwide), 4),
                     n_pixel_draws = count_mapwide)
write.csv(cw_tab, "output/tables/smaa_central_weights_by_class.csv", row.names = FALSE)

smaa_summary<- data.frame(
  metric = c("mean confidence factor", "median confidence factor",
             "% pixels confidence > 0.8", "% pixels confidence < 0.5",
             "mean P(class >= High)", "% pixels P(>=High) > 0.8"),
  value  = c(round(mean(confidence), 3), round(median(confidence), 3),
             round(100 * mean(confidence > 0.8), 1),
             round(100 * mean(confidence < 0.5), 1),
             round(mean(p_high), 3),
             round(100 * mean(p_high > 0.8), 1))
)
write.csv(smaa_summary, "output/tables/smaa_summary.csv", row.names = FALSE)
cat("\nSMAA summary:\n"); print(smaa_summary)
cat("\nMap-wide central weights per class (what weighting lands a pixel there):\n")
print(cw_tab[, c("label", factors)])

## 6. Figures
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# (a) the five class-acceptability maps
png("output/figures/smaa_acceptability.png", width = 1900, height = 850, res = 150)
op<- par(mfrow = c(1, 5), mar = c(2, 2, 3, 3))
labs<- c("Very Low","Low","Moderate","High","Very High")
for (c in 1:5) {
  plot(mk_raster(accept[, c]), main = paste0("Acceptability: ", labs[c]),
       col = hcl.colors(20, "Blues", rev = TRUE), range = c(0, 1))
}
par(op); dev.off()

# (b) confidence factor + policy category
png("output/figures/smaa_confidence.png", width = 1800, height = 900, res = 150)
op<- par(mfrow = c(1, 2), mar = c(2, 2, 3, 4))
plot(mk_raster(confidence), main = "SMAA confidence factor (modal class)",
     col = hcl.colors(20, "Blues", rev = TRUE), range = c(0, 1))
plot(mk_raster(p_high), main = "P(class >= High)",
     col = hcl.colors(20, "YlOrRd", rev = TRUE), range = c(0, 1))
par(op); dev.off()

cat("\n13_smaa.R complete.\n")
cat("Tables : smaa_summary.csv, smaa_central_weights_by_class.csv\n")
cat("Figures: smaa_acceptability.png, smaa_confidence.png\n")
cat("Note   : rank-1 acceptability should match 12's P(modal) exactly (checked above).\n")
