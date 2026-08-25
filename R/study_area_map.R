# ================================================================
# Study-area location map: Cox's Bazar district + Bangladesh inset
# Self-contained: boundaries are downloaded automatically (GADM).
#   install.packages(c("geodata","sf","terra","ggplot2","ggspatial","cowplot"))
# ================================================================
library(geodata); library(sf); library(terra)
library(ggplot2); library(ggspatial); library(cowplot)

dir.create("gadm", showWarnings = FALSE)

## 1. Download Bangladesh boundaries (GADM) ------------------------------
bgd0 <- sf::st_as_sf(geodata::gadm("BGD", level = 0, path = "gadm"))  # country
bgd2 <- sf::st_as_sf(geodata::gadm("BGD", level = 2, path = "gadm"))  # districts

## 2. Cox's Bazar district ----------------------------------------------
cox    <- bgd2[grepl("Cox", bgd2$NAME_2, ignore.case = TRUE), ]
cox_ll <- sf::st_transform(cox, 4326)

## 3. Analysis-grid extent (read straight from one of your rasters) ------
r <- terra::rast("data/LSM Tutorial/slope.tif")          # <-- any aligned raster, EPSG:32646

# safer approach: use terra to turn the extent into a polygon, then to an sf object
grid <- sf::st_as_sf(terra::as.polygons(terra::ext(r), crs = terra::crs(r)))
grid_ll <- sf::st_transform(grid, 4326)

## 4. Main panel --------------------------------------------------------
main <- ggplot() +
  geom_sf(data = cox_ll,  fill = "grey95", colour = "black", linewidth = 0.4) +
  geom_sf(data = grid_ll, fill = NA,       colour = "red",   linewidth = 0.7) +
  annotation_scale(location = "bl", width_hint = 0.35) +
  annotation_north_arrow(location = "tr", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  theme_minimal(base_size = 11) + labs(x = NULL, y = NULL)


## 5. Bangladesh inset (with an opaque white background + border) --------
inset <- ggplot() +
  geom_sf(data = bgd0, fill = "grey85", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = sf::st_geometry(sf::st_centroid(sf::st_union(cox_ll))),
          colour = "red", size = 2.5) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "black",
                                       linewidth = 0.6),
        plot.margin = margin(3, 3, 3, 3))

## 6. Combine + save (matches the \includegraphics box in 1.2) ----------
dir.create("figures", showWarnings = FALSE)
final <- ggdraw(main) +
  draw_plot(inset, x = 0.66, y = 0.66, width = 0.32, height = 0.32)
ggsave("output/figures/study_area_map.png", final, width = 7, height = 8, dpi = 300)
cat("Saved output/figures/study_area_map.png\n")
