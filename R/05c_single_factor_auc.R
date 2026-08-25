# R/05c_single_factor_auc.R
# Robustness: Single-factor discrimination check.
#
# What it does: A quick reality check. How well does each factor predict landslides 
# ON ITS OWN, ignoring the expert weights entirely? 
# Continuous factors get raw AUC, categorical ones get class lift.
# Note: Need to drop the raster attribute table so terra gives me raw codes.

suppressMessages({library(terra); library(pROC)})

hazard  <- rast("output/hazard_continuous.tif")          # defines the scored frame
inventory  <- project(vect("data/LSM Tutorial/Landslide_Inventory.shp"), crs(hazard))
base <- "data/LSM Tutorial/"

prep <- function(tif) {                                # raw values on the hazard frame
  r <- rast(paste0(base, tif))
  levels(r) <- NULL                                    # drop category table -> raw codes
  mask(crop(r, hazard), hazard)                              # keep only cells where hazard is defined
}
vals <- function(r) { p  <- extract(r, inventory, ID = FALSE)[, 1]; p <- p[is.finite(p)]
                      bg <- as.numeric(values(r));     bg <- bg[is.finite(bg)]
                      list(p = p, bg = bg) }

cat("Continuous single-factor AUC (raw value, higher = scored):\n")
for (tif in c("Slope.tif", "NDVI.tif", "Distance to road (m).tif")) {
  v <- vals(prep(tif))
  a <- as.numeric(auc(roc(c(rep(1, length(v$p)), rep(0, length(v$bg))),
                          c(v$p, v$bg), quiet = TRUE, direction = "<")))
  cat(sprintf("  %-26s AUC %.3f  (n_presence = %d)\n", tif, a, length(v$p)))
}

class_lift <- function(tif, labels) {
  v <- vals(prep(tif))
  cat(sprintf("\n%s per-class landslide lift (n = %d):\n", tif, length(v$p)))
  for (k in names(labels)) {
    kk <- as.integer(k)
    ipct <- 100 * mean(v$p == kk); bpct <- 100 * mean(v$bg == kk)
    cat(sprintf("  %-30s inv%%=%5.1f  bg%%=%5.1f  lift=%.2f\n",
                labels[[k]], ipct, bpct, ifelse(bpct > 0, ipct / bpct, NA)))
  }
}
class_lift("Soil Type.tif",
           list("0"="Noncalc. alluvium","1"="Noncalc. grey floodplain",
                "2"="Noncalc. brown floodplain","3"="Acid sulphate",
                "4"="Grey piedmont","5"="Brown hill"))
class_lift("Land Cover.tif",
           list("0"="Bare land","1"="Built-up","2"="Crop","3"="Vegetation","4"="Water"))
