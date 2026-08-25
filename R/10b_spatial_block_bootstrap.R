# R/10b_spatial_block_bootstrap.R
# Robustness: Spatial block bootstrap for AUC.
#
# Why: Standard bootstrap treats points as independent, which is too optimistic 
# since landslides cluster spatially. Resampling spatial blocks (e.g. 5km) gives 
# a more honest (and wider) confidence interval.

suppressMessages({library(terra)})

hazard    <- rast("output/hazard_continuous.tif")   # baseline WLC surface, [1,3]
inventory <- project(vect("data/LSM Tutorial/Landslide_Inventory.shp"), crs(hazard))

# presence scores + all valid background scores
presence <- extract(hazard, inventory, ID = FALSE)[, 1]
keep <- is.finite(presence); presence <- presence[keep]   # 355 in-frame points
bg   <- as.numeric(values(hazard)); bg <- sort(bg[is.finite(bg)])
xy   <- crds(inventory)[keep, ]

# fast rank-based AUC of a presence set vs the fixed background:
# AUC = mean over presence of the mid-rank fraction of its score in bg
mid_rank <- function(p) {
  low  <- findInterval(p, bg, left.open = TRUE)     # count bg < p
  high <- findInterval(p, bg)                        # count bg <= p
  ((low + high) / 2) / length(bg)
}
f0 <- mid_rank(presence)
cat(sprintf("point AUC = %.4f  (n_presence = %d)\n", mean(f0), length(presence)))

B <- 3000
boot_ci <- function(v) c(mean = mean(v), low = unname(quantile(v, .025)),
                         high = unname(quantile(v, .975)))
res <- data.frame(label = character(), mean = numeric(),
                  low = numeric(), high = numeric())

# (a) ordinary (independent) bootstrap -- reproduces 10_bootstrap.R
set.seed(1)
ordinary <- replicate(B, mean(sample(f0, replace = TRUE)))
c0 <- boot_ci(ordinary)
cat(sprintf("ordinary   : mean %.3f  CI [%.3f, %.3f]  width %.3f\n",
            c0[1], c0[2], c0[3], c0[3] - c0[2]))
res <- rbind(res, data.frame(label = "Ordinary (independent)",
                             mean = c0[1], low = c0[2], high = c0[3]))

# (b) spatial block bootstrap at several block sizes (metres)
for (bs in c(3000, 5000, 10000)) {
  blk    <- paste(floor(xy[, 1] / bs), floor(xy[, 2] / bs))
  groups <- split(seq_along(blk), blk)
  keys   <- names(groups); nb <- length(keys)
  set.seed(1)
  spatial <- replicate(B, {
    pick <- sample(keys, nb, replace = TRUE)
    mean(f0[unlist(groups[pick], use.names = FALSE)])
  })
  ci <- boot_ci(spatial)
  cat(sprintf("spatial %2dkm (%d blocks): mean %.3f  CI [%.3f, %.3f]  width %.3f\n",
              bs / 1000, nb, ci[1], ci[2], ci[3], ci[3] - ci[2]))
  res <- rbind(res, data.frame(label = sprintf("Spatial %d km", bs / 1000),
                               mean = ci[1], low = ci[2], high = ci[3]))
}

# ---- forest plot of the intervals ----
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
png("output/figures/spatial_bootstrap.png", width = 1600, height = 650, res = 200)
op <- par(mar = c(4, 12, 1, 6))
yy <- rev(seq_len(nrow(res)))
plot(NA, xlim = c(0.35, 0.95), ylim = c(0.5, nrow(res) + 0.5),
     xlab = "AUC", ylab = "", yaxt = "n", bty = "n")
abline(v = 0.5, lty = 2, col = "grey40")
for (i in seq_len(nrow(res))) {
  col <- if (i == 1) "grey40" else "#2c7fb8"
  segments(res$low[i], yy[i], res$high[i], yy[i], lwd = 3, col = col)
  points(res$mean[i], yy[i], pch = 19, col = col, cex = 1.1)
  text(res$high[i] + 0.01, yy[i],
       sprintf("[%.2f, %.2f]", res$low[i], res$high[i]), pos = 4, cex = .8, xpd = NA)
}
axis(2, at = yy, labels = res$label, las = 1, tick = FALSE)
par(op); dev.off()
cat("saved output/figures/spatial_bootstrap.png\n")
