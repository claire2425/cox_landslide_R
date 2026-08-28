# R/06_make_map.R
# Step 6: GIS-style map figures.
#   (a) hazard_map.png        - 5-class hazard map
#   (b) inventory_overlay.png - hazard + inventory + Rohingya camp

cat("\n--- 06_make_map: produce GIS-style figures ---\n")

suppressPackageStartupMessages({ library(terra); library(sf); library(tmap) })

if (!file.exists("output/hazard_class.rds")) stop("Run R/04_hazard_overlay.R first.")
hazard_class<- unwrap(readRDS("output/hazard_class.rds"))
boundary<- readRDS("output/boundary.rds")
inventory<- readRDS("output/inventory.rds")
camp<- readRDS("output/camp.rds")

# Make the raster categorical so tmap shows discrete classes with labels.
labs<- c("Very Low","Low","Moderate","High","Very High")
levels(hazard_class)<- data.frame(value = 1:5, class = labs)

pal<- c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")  # blue->red
inv_pts<- tryCatch(st_centroid(inventory), error = function(e) inventory)

tmap_mode("plot")
tmap_v4<- utils::packageVersion("tmap") >= "4.0.0"
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

cat(sprintf("Detected tmap %s -> using %s syntax.\n",
            utils::packageVersion("tmap"), if (tmap_v4) "v4" else "v3"))

# ---- Figure (a): hazard map -------------------------------------------------
if (tmap_v4) {
  fig_a<- tm_shape(hazard_class) +
    tm_raster(col.scale = tm_scale_categorical(values = pal),
              col.legend = tm_legend(title = "Landslide hazard")) +
    tm_shape(boundary) + tm_borders(col = "black", lwd = 1.2) +
    tm_shape(camp)     + tm_borders(col = "blue",  lwd = 1.5) +
    tm_title("Cox's Bazar - Landslide Susceptibility (baseline AHP)") +
    tm_scalebar(position = c("left","bottom")) +
    tm_compass(position = c("right","top")) +
    tm_layout(legend.outside = TRUE, frame = TRUE)
} else {
  fig_a<- tm_shape(hazard_class) +
    tm_raster(palette = pal, title = "Landslide hazard", labels = labs, style = "cat") +
    tm_shape(boundary) + tm_borders(col = "black", lwd = 1.2) +
    tm_shape(camp)     + tm_borders(col = "blue",  lwd = 1.5) +
    tm_layout(main.title = "Cox's Bazar - Landslide Susceptibility (baseline AHP)",
              main.title.size = 1.0, legend.outside = TRUE, frame = TRUE) +
    tm_scale_bar(position = c("left","bottom")) +
    tm_compass(position = c("right","top"))
}
tmap_save(fig_a, "output/figures/hazard_map.png", width = 7, height = 6, dpi = 200)
cat("Saved: output/figures/hazard_map.png\n")

# ---- Figure (b): overlay ----------------------------------------------------
if (tmap_v4) {
  fig_b<- tm_shape(hazard_class) +
    tm_raster(col.scale = tm_scale_categorical(values = pal),
              col.legend = tm_legend(title = "Landslide hazard")) +
    tm_shape(boundary) + tm_borders(col = "black", lwd = 1.2) +
    tm_shape(camp)     + tm_borders(col = "blue",  lwd = 1.5) +
    tm_shape(inv_pts)  + tm_dots(fill = "black", size = 0.15) +
    tm_title("Hazard with landslide inventory + Rohingya camp") +
    tm_scalebar(position = c("left","bottom")) +
    tm_layout(legend.outside = TRUE, frame = TRUE)
} else {
  fig_b<- tm_shape(hazard_class) +
    tm_raster(palette = pal, title = "Landslide hazard", labels = labs, style = "cat") +
    tm_shape(boundary) + tm_borders(col = "black", lwd = 1.2) +
    tm_shape(camp)     + tm_borders(col = "blue",  lwd = 1.5) +
    tm_shape(inv_pts)  + tm_dots(col = "black", size = 0.15) +
    tm_layout(main.title = "Hazard with landslide inventory + Rohingya camp",
              main.title.size = 0.9, legend.outside = TRUE, frame = TRUE) +
    tm_scale_bar(position = c("left","bottom"))
}
tmap_save(fig_b, "output/figures/inventory_overlay.png", width = 8, height = 6, dpi = 200)
cat("Saved: output/figures/inventory_overlay.png\n\n")

cat("All figures produced. Open output/figures/ to view.\n\n")

