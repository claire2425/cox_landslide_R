# ================================================================
# Figure 2.1: Landslide inventory and Cox's Bazar study area
#
# This script is located in:
#   cox_landslide_R/R/
#
# It therefore detects the project root automatically and writes
# the output to:
#   cox_landslide_R/output/figures/
#
# Figure 1.1 is NOT affected.
# ================================================================

suppressPackageStartupMessages({
  library(geodata)
  library(sf)
  library(terra)
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
})

# ---------------------------------------------------------------
# 1. Automatically determine project root
# ---------------------------------------------------------------

script_dir <- normalizePath(
  getwd(),
  winslash = "/"
)

# This script is expected to be run from cox_landslide_R/R/
# Move one level up to obtain the project root.
project_root <- normalizePath(
  file.path(script_dir, ".."),
  winslash = "/",
  mustWork = TRUE
)

data_dir <- file.path(
  project_root,
  "data",
  "LSM Tutorial"
)

output_dir <- file.path(
  project_root,
  "output",
  "figures"
)

# ---------------------------------------------------------------
# 2. Check paths
# ---------------------------------------------------------------

if (!dir.exists(data_dir)) {
  stop(
    "Could not find the data directory:\n",
    data_dir,
    "\n\nCurrent working directory:\n",
    getwd()
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Project root:\n", project_root, "\n")
cat("Data directory:\n", data_dir, "\n")
cat("Output directory:\n", output_dir, "\n\n")

# ---------------------------------------------------------------
# 3. Bangladesh boundaries
# ---------------------------------------------------------------

dir.create(
  file.path(project_root, "R", "gadm"),
  showWarnings = FALSE
)

bgd0 <- sf::st_as_sf(
  geodata::gadm(
    "BGD",
    level = 0,
    path = file.path(project_root, "R", "gadm")
  )
)

bgd2 <- sf::st_as_sf(
  geodata::gadm(
    "BGD",
    level = 2,
    path = file.path(project_root, "R", "gadm")
  )
)

# ---------------------------------------------------------------
# 4. Cox's Bazar district
# ---------------------------------------------------------------

cox <- bgd2[
  grepl("Cox", bgd2$NAME_2, ignore.case = TRUE),
]

if (nrow(cox) == 0) {
  stop("Cox's Bazar district could not be identified.")
}

cox_ll <- sf::st_transform(
  cox,
  4326
)

# ---------------------------------------------------------------
# 5. Analysis-grid extent
#
# Same raster used by your original study_area_map.R.
# ---------------------------------------------------------------

slope_path <- file.path(
  data_dir,
  "slope.tif"
)

if (!file.exists(slope_path)) {
  stop(
    "Could not find slope raster:\n",
    slope_path
  )
}

r <- terra::rast(slope_path)

grid <- sf::st_as_sf(
  terra::as.polygons(
    terra::ext(r),
    crs = terra::crs(r)
  )
)

grid_ll <- sf::st_transform(
  grid,
  4326
)

# ---------------------------------------------------------------
# 6. Landslide inventory
# ---------------------------------------------------------------

inventory_path <- file.path(
  data_dir,
  "Landslide_Inventory.shp"
)

if (!file.exists(inventory_path)) {
  stop(
    "Could not find landslide inventory:\n",
    inventory_path
  )
}

inventory <- sf::st_read(
  inventory_path,
  quiet = TRUE
)

if (nrow(inventory) == 0) {
  stop("Landslide inventory contains no features.")
}

inventory_ll <- sf::st_transform(
  inventory,
  4326
)

cat(
  "Landslide inventory features loaded:",
  nrow(inventory_ll),
  "\n"
)

# ---------------------------------------------------------------
# 7. Main map
#
# No susceptibility/hazard raster is shown.
# ---------------------------------------------------------------

main <- ggplot() +
  
  geom_sf(
    data = cox_ll,
    fill = "grey95",
    colour = "black",
    linewidth = 0.5
  ) +
  
  geom_sf(
    data = grid_ll,
    fill = NA,
    colour = "grey45",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  
  geom_sf(
    data = inventory_ll,
    colour = "black",
    size = 0.8,
    alpha = 0.75
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.32
  ) +
  
  annotation_north_arrow(
    location = "tr",
    style = north_arrow_minimal(),
    height = unit(1, "cm"),
    width = unit(1, "cm")
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  coord_sf() +
  
  theme_minimal(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

# ---------------------------------------------------------------
# 8. Bangladesh inset
# ---------------------------------------------------------------

inset <- ggplot() +
  
  geom_sf(
    data = bgd0,
    fill = "grey90",
    colour = "grey45",
    linewidth = 0.3
  ) +
  
  geom_sf(
    data = cox_ll,
    fill = NA,
    colour = "black",
    linewidth = 0.8
  ) +
  
  theme_void() +
  
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.6
    ),
    plot.margin = margin(3, 3, 3, 3)
  )

# ---------------------------------------------------------------
# 9. Combine
# ---------------------------------------------------------------

final <- cowplot::ggdraw(main) +
  cowplot::draw_plot(
    inset,
    x = 0.66,
    y = 0.66,
    width = 0.30,
    height = 0.30
  )

# ---------------------------------------------------------------
# 10. Save
# ---------------------------------------------------------------

output_path <- file.path(
  output_dir,
  "study_area_inventory_map.png"
)

ggsave(
  filename = output_path,
  plot = final,
  width = 7,
  height = 8,
  dpi = 300
)

cat(
  "\nSaved Figure 2.1 to:\n",
  output_path,
  "\n"
)
