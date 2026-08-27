# R/12_hazard_posterior.R
# Step 12: Per-pixel hazard posterior (RQ3b).
#
# What it does: Takes the 8000 posterior weight draws from Stan and pushes them 
# through the spatial overlay to see how confident we are at EACH pixel.
#
# Memory trick: Don't hold 8000 maps in RAM at once. I'm streaming them draw by draw 
# and just updating the class counts per pixel.
#
# Outputs: Modal class (best guess) and P(modal) (how confident I am in it).

cat("\n--- 12_hazard_posterior: per-pixel hazard posterior + confidence ---\n")
suppressMessages({ library(terra) })
set.seed(2026)

## 1. Load the posterior weights and the reclassified factor stack
W<- readRDS("output/bayes_weights.rds")          # M x 5
factors<- colnames(W); M<- nrow(W)
cat(sprintf("Loaded %d posterior weight draws over %d factors.\n", M, ncol(W)))

# The reclassified 1-5 factor layers the baseline overlay multiplies by weights.
# Load them in the SAME ORDER as the columns of W. Adjust paths/names to match
# your 04_hazard_overlay inputs.
stack_path<- "output/reclass_stack.rds"
if (file.exists(stack_path)) {
  rs<- readRDS(stack_path)
  if (!inherits(rs, "SpatRaster")) rs<- rast(rs)
} else {
  # fallback: build from individual reclassified rasters
  rs<- rast(list(
    Slope     = rast("output/reclass_Slope.tif"),
    LandCover = rast("output/reclass_LandCover.tif"),
    NDVI      = rast("output/reclass_NDVI.tif"),
    SoilType  = rast("output/reclass_SoilType.tif"),
    DistRoad  = rast("output/reclass_DistRoad.tif")
  ))
}
# reorder layers to match W's columns exactly
rs<- rs[[ match(factors, names(rs)) ]]
stopifnot(all(names(rs) == factors))

## 2. The baseline class breaks (so the posterior uses the SAME 5-class cuts)
# IMPORTANT: use the baseline's fixed class boundaries, NOT per-draw quantiles.
# Otherwise each draw would be re-binned on its own scale and the confidence
# map would conflate "weight uncertainty" with "moving the goalposts".
#
# We need the 6 numeric break points that 04/05 used to cut the continuous
# baseline hazard into 5 classes. Depending on how your pipeline stored things,
# get them from ONE of these (the script tries them in order):
#   (a) output/hazard_breaks.rds        -- if you saved the breaks explicitly
#   (b) output/hazard_class.rds$breaks  -- if the class object carries them
#   (c) re-derive from the baseline CONTINUOUS hazard (hazard_continuous.rds)
#       using the SAME quantile rule 04_hazard_overlay used (edit probs= to match)
get_breaks<- function() {
  if (file.exists("output/hazard_breaks.rds")) {
    b<- readRDS("output/hazard_breaks.rds")
    if (is.numeric(b) && length(b) == 6) return(b)
  }
  if (file.exists("output/hazard_class.rds")) {
    hc<- readRDS("output/hazard_class.rds")
    if (is.list(hc) && !is.null(hc$breaks) && length(hc$breaks) == 6)
      return(hc$breaks)
  }
  if (file.exists("output/hazard_continuous.rds")) {
    hcont<- readRDS("output/hazard_continuous.rds")
    sc<- values(if (inherits(hcont, "SpatRaster")) hcont else rast(hcont))
    sc<- sc[is.finite(sc)]
    # EDIT probs to match 04_hazard_overlay's class rule (here: equal quintiles)
    return(as.numeric(quantile(sc, probs = seq(0, 1, by = 0.2))))
  }
  stop("Could not obtain 6 class breaks. Save them to output/hazard_breaks.rds ",
       "in 04_hazard_overlay, or point this function at the right object.")
}
brk<- get_breaks()
stopifnot(is.numeric(brk), length(brk) == 6)
cat("Using fixed baseline class breaks:\n"); print(round(brk, 4))

# pull factor layers into a matrix of valid (non-NA) pixels for speed
vals<- values(rs)                                 # pixels x 5
ok<- stats::complete.cases(vals)                # pixels inside study area
Xok<- vals[ok, , drop = FALSE]                   # n_ok x 5
n_ok<- nrow(Xok)
cat(sprintf("%d valid pixels (of %d).\n", n_ok, nrow(vals)))

## 3. Stream the posterior: accumulate per-pixel class counts
# class_count[pixel, class] -- integer counts over the M draws
class_count<- matrix(0L, nrow = n_ok, ncol = 5)

cut_to_class<- function(score) {
  cl<- findInterval(score, brk, rightmost.closed = TRUE)
  cl[cl < 1]<- 1L; cl[cl > 5]<- 5L
  as.integer(cl)
}

cat("Streaming draws"); pb<- txtProgressBar(min = 0, max = M, style = 3)
# IMPORTANT consistency requirement (verify once against 04_hazard_overlay):
#   score = Xok %*% W[m, ]  is a PLAIN weighted SUM of the 1-5 reclassified
#   layers. This matches the baseline ONLY IF 04_hazard_overlay also forms a
#   plain weighted sum (sum_k w_k * layer_k) before cutting into classes, with
#   weights summing to 1. Softmax weights sum to 1, so this is consistent with a
#   weighted-sum/weighted-mean baseline. If your baseline rescales the score
#   (e.g. divides by something, or min-max stretches) BEFORE applying the class
#   breaks, replicate that same transform on `score` here, or the breaks will
#   mis-bin every pixel. One-line check: confirm overlay == sum(w*layer).
for (m in 1:M) {
  score<- as.numeric(Xok %*% W[m, ])               # weighted overlay, this draw
  cl<- cut_to_class(score)                       # 1..5 per pixel
  # bump the count for each pixel's class (vectorised)
  class_count[cbind(seq_len(n_ok), cl)] <-
    class_count[cbind(seq_len(n_ok), cl)] + 1L
  setTxtProgressBar(pb, m)
}
close(pb)

## 4. Per-pixel confidence summaries
prob<- class_count / M                       # n_ok x 5 class probabilities
modal_class<- max.col(prob, ties.method = "first")  # 1..5
p_modal<- prob[cbind(seq_len(n_ok), modal_class)]
p_high<- prob[, 4] + prob[, 5]                 # P(High or Very-High)
ent<- -rowSums(ifelse(prob > 0, prob * log(prob), 0)) / log(5)  # 0..1

# write each summary back into a raster shaped like the input
mk_raster<- function(vec) {
  r<- rs[[1]]; r[]<- NA_real_; r[ok]<- vec; r
}
r_modal<- mk_raster(modal_class)
# label the five classes by name (matches the text: Very Low ... Very High)
levels(r_modal)<- data.frame(value = 1:5,
                             class = c("Very Low", "Low", "Moderate", "High", "Very High"))
r_pmod<- mk_raster(p_modal)
r_phigh<- mk_raster(p_high)
r_ent<- mk_raster(ent)

saveRDS(r_modal, "output/conf_modalclass.rds")
saveRDS(r_pmod,  "output/conf_pmodal.rds")
saveRDS(r_phigh, "output/conf_phigh.rds")
saveRDS(r_ent,   "output/conf_entropy.rds")

## 5. Headline figure: hazard (modal) beside confidence (P modal)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
png("output/figures/confidence_map.png", width = 1800, height = 900, res = 150)
op<- par(mfrow = c(1, 2), mar = c(2, 2, 3, 4))
plot(r_modal, main = "Posterior modal susceptibility class",
     col = hcl.colors(5, "YlOrRd", rev = TRUE))
plot(r_pmod,  main = "Decision confidence  P(modal class)",
     col = hcl.colors(20, "Blues", rev = TRUE))
par(op); dev.off()

## 6. Summary table
summ<- data.frame(
  metric = c("mean P(modal)", "median P(modal)",
             "% pixels P(modal) > 0.8", "% pixels P(modal) < 0.5",
             "mean normalised entropy"),
  value  = c(round(mean(p_modal), 3), round(median(p_modal), 3),
             round(100 * mean(p_modal > 0.8), 1),
             round(100 * mean(p_modal < 0.5), 1),
             round(mean(ent), 3))
)
write.csv(summ, "output/tables/confidence_summary.csv", row.names = FALSE)
cat("\nConfidence summary:\n"); print(summ)

# Expected qualitative check (RQ3 narrative): low-confidence pixels should
# cluster where Soil/Slope dominate and take ambiguous values -- i.e. exactly
# where the SoilType flip ambiguity lived. Inspect r_pmod against the Soil layer.

cat("\n12_hazard_posterior.R complete.\n")

